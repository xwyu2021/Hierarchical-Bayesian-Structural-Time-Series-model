run_flex_bsts <- function(data=bsts.zoo.qua,pre.period=pre.period.qua,
                          post.period=post.period.qua,
                          alpha=0.05,model.args=list(nseasons=1,
                                                     prior.level.sd = 0.01,
                                                     standardize.data=TRUE),
                          trend, modelsize,prior.inc=NULL,
                          cycles=c(1),
                          niter=5000,random=FALSE,realy=rep(0,50)){
  pre.period = GetPeriodIndices(pre.period, time(data))
  post.period = GetPeriodIndices(post.period, time(data))
  times <- time(data)
  time(data) <- seq_len(nrow(data))
  
  # Zoom in on data in modeling range, remove original time indices.
  pre.period[1] <- max(pre.period[1], which.max(!is.na(data[, 1])))
  data.modeling <- window(data, start = pre.period[1])
  times.modeling <- window(times, start = pre.period[1])
  if (is.null(ncol(data.modeling))) {
    dim(data.modeling) <- c(length(data.modeling), 1)
  }
  
  # Standardize all variables?
  UnStandardize <- identity
  if (model.args$standardize.data) {
    fit.range <- c(1, diff(pre.period) + 1)
    sd.results <- StandardizeAllVariables(data.modeling, fit.range)
    data.modeling <- sd.results$data
    UnStandardize <- sd.results$UnStandardize
  }
  
  # Set observed response after pre-period to NA.
  window(data.modeling[, 1], start = pre.period[2] + 1) <- NA
  sd.prior <- SdPrior(sigma.guess = 0.01 * 1,
                      upper.limit = 1,
                      sample.size = 32)
  
  # Construct model and perform inference
  if(trend == 'LLT'){
    ss.train <- AddLocalLinearTrend(list(),y=data.modeling[,1],
                                    level.sigma.prior= sd.prior,
                                    slope.sigma.prior= sd.prior,
                                    initial.level.prior = NormalPrior(0,1),
                                    initial.slope.prior = NormalPrior(0,1))
  }else{
    if(trend == 'AR1'){
      ss.train <- AddAr(list(),y=data.modeling[,1],sigma.prior= sd.prior,
                        initial.state.prior=NormalPrior(0,1))
    }else{
      if(trend == 'LL'){
        ss.train <- AddLocalLevel(list(),y=data.modeling[,1],
                                  sigma.prior= sd.prior,
                                  initial.state.prior = NormalPrior(0,1))
      }else{
        if(trend == 'GLLT'){
          ss.train <-  AddSemilocalLinearTrend(list(),y=data.modeling[,1],
                                               level.sigma.prior = sd.prior,
                                               slope.sigma.prior = sd.prior,
                                               initial.level.prior = NormalPrior(0,1),
                                               initial.slope.prior = NormalPrior(0,1))
        }
      }
    }
  }
  if (model.args$nseasons > 1) {
    ss.train <- AddSeasonal(ss.train, y=data.modeling[,1],
                            nseasons = model.args$nseasons)
  }
  formula <- paste0(names(data)[1], " ~ .")
  if(!is.null(modelsize)){
    modeli <- bsts(formula,data=data.modeling,ss.train,
                   niter = niter, expected.r2 = 0.8, 
                   prior.df = 50,
                   expected.model.size = modelsize)
  }else{
    modeli <- bsts(formula,data=data.modeling,ss.train,
                   niter = niter, expected.r2 = 0.8, 
                   prior.df = 50,
                   prior.inclusion.probabilities = prior.inc)
  }
  
  time(modeli$original.series) <- time(data)
  
  y.cf <- window(data[, 1], start = pre.period[2] + 1)
  
  if(random == TRUE){
    y.cf <- realy
  }
  inferences <- CompilePosteriorInferences(modeli, y.cf,
                                           post.period - pre.period[1] + 1,
                                           alpha, UnStandardize)
  
  temporaleffect <- get.point.effect(bsts.model=modeli,y.cf,post.period, alpha,
                                     UnStandardize,cycles)
  time(inferences$series) <- times.modeling
  
  empty <- zoo(, times)
  inferences$series <- merge(inferences$series, empty, all = TRUE)
  assert_that(nrow(inferences$series) == nrow(data))
  
  # Replace <y.model> by full original response
  inferences$series[, 1] <- data[, 1]
  
  # Assign response-variable names
  names(inferences$series)[1] <- "response"
  names(inferences$series)[2] <- "cum.response"
  
  # Return 'CausalImpact' object
  model <- list(pre.period = times[pre.period],
                post.period = times[post.period],
                model.args = model.args,
                bsts.model = modeli,
                alpha = alpha,
                posterior.samples = inferences$posterior.samples)
  impact <- list(series = inferences$series,
                 summary = inferences$summary,
                 report = inferences$report,
                 model = model,
                 temporaleffect=temporaleffect,UnStandardize=UnStandardize)
  class(impact) <- "CausalImpact"
  return(impact)
  
}





effect.ssbsts <- function(modeli,y.cf,post.period,UnStandardize,alpha=0.05){
  
  inferences <- CompilePosteriorInferences(modeli, y.cf,
                                           post.period,
                                           alpha, UnStandardize)
  return(inferences$summary)
}






get.point.effect <- function(bsts.model=modeli,y.cf,post.period, alpha,
                             UnStandardize,cycles){
  y.samples <- ComputeResponseTrajectories(bsts.model)
  state.samples <- GetPosteriorStateSamples(bsts.model)
  point.pred <- ComputePointPredictions(y.samples, state.samples, alpha)
  # point pred mean coming from mean of states, ignoring noises
  # point pred lower/upper coming from posterior y.samples
  y.samples <- UnStandardize(y.samples)
  point.pred <- UnStandardize(point.pred)
  y.model <- UnStandardize(bsts.model$original.series)
  
  indices <- seq_along(y.model)
  is.cf.period <- (indices >= length(y.model) - length(y.cf) + 1)
  y.model[is.cf.period] <- y.cf
  
  is.post.period <- (indices >= post.period[1]) & (indices <= post.period[2])
  y.samples.post <- y.samples[, is.post.period, drop = FALSE]
  point.pred.mean.post <- point.pred$point.pred[is.post.period]
  y.post <- y.cf[tail(is.post.period, length(y.cf))]
  y.repmat.post <- matrix(y.post, nrow = nrow(y.samples.post),
                          ncol = length(y.post), byrow = TRUE)
  prob.lower <- alpha / 2      # e.g., 0.025 when alpha = 0.05
  prob.upper <- 1 - alpha / 2  # e.g., 0.975 when alpha = 0.05
  
  p_vals <- list()
  temp_effect <- list()
  for(k in 1:length(cycles)){
    cycle <- cycles[k]
    p_vals[[k]] <- find_p(y.post,y.samples.post,cycle)
    if(cycle > 1){
      temp_effect[[k]] <- find_effect(y.post,y.samples.post,point.pred.mean.post,y.repmat.post,cycle,prob.lower,prob.upper)
    }
  }
  return(list(p_vals=p_vals,temp_effect=temp_effect))
}



find_p <- function(y.post,y.samples.post,cycle){
  if(cycle==1){
    p <- rep(0,length(y.post))
    for(j in 1:length(y.post)){
      x = as.numeric(y.post[j])
      p[j] <- min(sum(c(y.samples.post[,j],y.post[j]) >=x),
                        sum(c(y.samples.post[,j],y.post[j]) <=x))/(length(y.samples.post[,j])+1)
    }
  }else{
    p <- rep(0,ceiling(length(y.post)/cycle))
    for(j in 1:length(p)){
      index <- c((j-1)* cycle+1,min(j*cycle,length(y.post)))
      if(index[1] == index[2]){
        y.samples.post.sum.sub <- rowSums(y.samples.post[,drop=FALSE,index[1]])
        y.post.sum.sub <- sum(y.post[index[1]])
      }else{
        y.samples.post.sum.sub <- rowSums(y.samples.post[,c(index[1]:index[2])])
        y.post.sum.sub <- sum(y.post[c(index[1]:index[2])])
      }
      
      p[j] <- min(sum(c(y.samples.post.sum.sub,y.post.sum.sub) >= y.post.sum.sub),
                       sum(c(y.samples.post.sum.sub,y.post.sum.sub) <= y.post.sum.sub))/
        (length(y.samples.post.sum.sub)+1)
    }
  }
  return(p)
}


find_effect <- function(y.post,y.samples.post,point.pred.mean.post,y.repmat.post,cycle,prob.lower,prob.upper){
  temp.effect <- temp.effect.lower <- temp.effect.upper <- rep(0,ceiling(length(y.post)/cycle))
  for(j in 1:length(temp.effect)){
    index <- c((j-1)* cycle+1,min(j*cycle,length(y.post)))
    if(index[1] == index[2]){
      y.samples.post.sub <- y.samples.post[,drop=FALSE,index[1]]
      y.post.sub <- y.post[index[1]]
      y.repmat.post.sub <- y.repmat.post[,drop=FALSE,index[1]]
      point.pred.mean.post.sub <- point.pred.mean.post[index[1]]
    }else{
      y.samples.post.sub <- y.samples.post[,c(index[1]:index[2])]
      y.post.sub <- y.post[c(index[1]:index[2])]
      y.repmat.post.sub <- y.repmat.post[,c(index[1]:index[2])]
      point.pred.mean.post.sub <- point.pred.mean.post[c(index[1]:index[2])]
    }
    temp.effect[j] = mean(y.post.sub) - mean(point.pred.mean.post.sub)
    temp.effect.lower[j] <- quantile(rowMeans(y.repmat.post.sub - y.samples.post.sub),
                                     prob.lower)
    temp.effect.upper[j] <- quantile(rowMeans(y.repmat.post.sub - y.samples.post.sub),
                                     prob.upper)
  }
  return(list(tempeffect = temp.effect,
              tempeffect.lower=temp.effect.lower,
              tempeffect.upper=temp.effect.upper))
}






