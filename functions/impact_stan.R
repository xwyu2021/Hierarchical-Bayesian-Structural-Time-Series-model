stan_mod_ll <- stan_model(file='algorithms (stan)/bstshsSAVS_img_ll.stan',verbose=T)
stan_mod_ar <- stan_model(file='algorithms (stan)/bstshsSAVS_img.stan',verbose=T)
stan_mod_llt <- stan_model(file='algorithms (stan)/bstshsSAVS_img_llt.stan',verbose=T)
stan_mod_gllt <- stan_model(file='algorithms (stan)/bstshsSAVS_img_gllt.stan',verbose=T)

run_stanbsts <- function(pre.period,post.period,data=bsts.zoo.qua,modelsize,trend,
                         nchain=1,niter=5000,alpha=0.05,cycles,adapt=0.8,burn=3000,random=FALSE,realy=rep(0,50)){
  pre.time <- pre.period[2] - pre.period[1] + 1
  post.time <- post.period[2] - post.period[1] + 1
  origin.data <- as.matrix(data)
  sd.results <- StandardizeAllVariables(origin.data, c(1,pre.time))
  data.modeling <- sd.results$data
  UnStandardize <- sd.results$UnStandardize
  
  times <- time(data)
  pre.y.scale <- as.vector(data.modeling[1:pre.time,1])
  x.scale <- data.modeling[,drop=FALSE,2:dim(data.modeling)[2]]
  dimx <- dim(x.scale)[2]
  y.cf <- as.vector(origin.data[c((pre.time+1):(pre.time+post.time)),1])
  y.model <- as.vector(origin.data[,1])
  if(random==TRUE){
    y.cf <- realy
    y.model[c((pre.time+1):(pre.time+post.time))] <- realy
  }
  tau0 <- modelsize/((dimx-modelsize)*sqrt(pre.time))
  if(modelsize==dimx){
    tau0 <-1
  }
  if(trend == 'LLT'){
    ls <- list(mu_aux=1,delta_aux=1,mu_err=rep(0,pre.time),delta_err=rep(0,pre.time))
    lls <- vector("list",nchain)
    init_fun <- lapply(lls,function(...) ls)
    modeli <- sampling(stan_mod_llt,
                       data = list( T =pre.time, T_forecast = post.time, M=dimx, 
                                    y = pre.y.scale, x =x.scale,
                                    sigma_itcpt = 0.1, scale_global=tau0,
                                    nu_global = 1, nu_local = 1, # nu_local=1 means horseshoe, nu_global=1 gives tau a half cauchy
                                    slab_scale = sqrt(0.2), 
                                    slab_df = 50,Meff=1,
                                    sdy=0.01,sdx=rep(1,dimx),c_df=50,tr_df=32),
                       iter=niter, chains=nchain, control = list(adapt_delta=adapt),
                       warmup=burn,init=init_fun)
  }else{
    if(trend == 'AR1'){
      ls <- list(mu_aux=1,mu_err=rep(0,pre.time),phi=0)
      lls <- vector("list",nchain)
      init_fun <- lapply(lls,function(...) ls)
      modeli <- sampling(stan_mod_ar,
                         data = list( T =pre.time, T_forecast = post.time, M=dimx, 
                                      y = pre.y.scale, x =x.scale,
                                      sigma_itcpt = 0.1, scale_global=tau0,
                                      nu_global = 1, nu_local = 1, # nu_local=1 means horseshoe, nu_global=1 gives tau a half cauchy
                                      slab_scale = sqrt(0.2), 
                                      slab_df = 50,Meff=1,
                                      sdy=0.01,sdx=rep(1,dimx),c_df=50,tr_df=32),
                         iter=niter, chains=nchain, control = list(adapt_delta=adapt),
                         warmup=burn,init=init_fun)
    }else{
      if(trend == 'LL'){
        ls <- list(mu_aux=1,mu_err=rep(0,pre.time))
        lls <- vector("list",nchain)
        init_fun <- lapply(lls,function(...) ls)
        #init_fun <- function(...) list(sigma_mu=0.01)
        modeli <- sampling(stan_mod_ll,
                       data = list( T =pre.time, T_forecast = post.time, M=dimx, 
                                    y = pre.y.scale, x =x.scale,
                                    sigma_itcpt = 0.1, scale_global=tau0,
                                    nu_global = 1, nu_local = 1, # nu_local=1 means horseshoe, nu_global=1 gives tau a half cauchy
                                    slab_scale = sqrt(0.2), 
                                    slab_df = 50,Meff=1,
                                    sdy=0.01,sdx=rep(1,dimx),c_df=50,tr_df=32),
                       iter=niter, chains=nchain, control = list(adapt_delta=adapt),
                       warmup=burn,init=init_fun)
      }else{
        if(trend == 'GLLT'){
          ls <- list(mu_aux=1,delta_aux=1,phi=0,D=0,mu_err=rep(0,pre.time),delta_err=rep(0,pre.time))
          lls <- vector("list",nchain)
          init_fun <- lapply(lls,function(...) ls)
          modeli <- sampling(stan_mod_gllt,
                             data = list( T =pre.time, T_forecast = post.time, M=dimx, 
                                          y = pre.y.scale, x =x.scale,
                                          sigma_itcpt = 0.1, scale_global=tau0,
                                          nu_global = 1, nu_local = 1, # nu_local=1 means horseshoe, nu_global=1 gives tau a half cauchy
                                          slab_scale = sqrt(0.2), 
                                          slab_df = 50,Meff=1,
                                          sdy=0.01,sdx=rep(1,dimx),c_df=50,tr_df=32),
                             iter=niter, chains=nchain, control = list(adapt_delta=adapt),
                             warmup=burn,init=init_fun) 
          
        }
      }
    }
  }
  fitted <- rstan::extract(modeli)
  y.samples <- cbind(fitted$y_fitted,fitted$y_forecast)
  state.samples <- cbind(fitted$mu+fitted$f, fitted$mu_forecast+fitted$f_forecast) 
  prob.lower <- alpha / 2      # e.g., 0.025 when alpha = 0.05
  prob.upper <- 1 - alpha / 2 
  point.pred.mean <-colMeans(state.samples)
  point.pred.lower <- as.numeric(t(apply(y.samples, 2, quantile, prob.lower)))
  point.pred.upper <- as.numeric(t(apply(y.samples, 2, quantile, prob.upper)))
  point.pred <- data.frame(point.pred = point.pred.mean,
                           point.pred.lower, point.pred.upper)
  y.samples <- UnStandardize(y.samples)
  point.pred <- UnStandardize(point.pred)
  
  cum.pred <- ComputeCumulativePredictions(y.samples, point.pred, y.model,
                                           pre.time+1, alpha)
  
  
  y.samples.post <- y.samples[, c((pre.time+1):(pre.time+post.time)), drop = FALSE]
  point.pred.mean.post <- point.pred$point.pred[c((pre.time+1):(pre.time+post.time))]
  y.post <- tail(y.cf, post.time)
  
  summary <- CompileSummaryTable(y.post, y.samples.post, point.pred.mean.post,
                                 alpha)
  temporaleffect <- get.point.effect.stan(y.samples.post,point.pred.mean.post,y.post,
                                          y.model,prob.lower,prob.upper,cycles)
  
  cum.y.model <- cumsum.na.rm(y.model)
  series <- zoo(data.frame(y.model, cum.y.model, point.pred, cum.pred),
                time(y.model))
  series$point.effect <- series$y.model - series$point.pred
  series$point.effect.lower <- series$y.model - series$point.pred.upper
  series$point.effect.upper <- series$y.model - series$point.pred.lower
  series$cum.effect <- series$cum.y.model - series$cum.pred
  series$cum.effect.lower <- series$cum.y.model - series$cum.pred.upper
  series$cum.effect.upper <- series$cum.y.model - series$cum.pred.lower
  
  time(series) <- times
  empty <- zoo(, times)
  series <- merge(series, empty, all = TRUE)
  # Replace <y.model> by full original response
  series[, 1] <- origin.data[, 1]
  # Assign response-variable names
  names(series)[1] <- "response"
  names(series)[2] <- "cum.response"
  model <- list(pre.period = times[c(1,pre.time)],
                post.period = c(tail(times,post.time)[1],last(times)),
                model.args = NULL,
                bsts.model = modeli,
                alpha = alpha,
                posterior.samples = y.samples)
  impact <- list(series = series,
                 summary = summary,
                 report=NULL,
                 model = model,
                 temporaleffect=temporaleffect,
                 UnStandardize=UnStandardize)
  class(impact) <- "CausalImpact"
  return(impact)
}

get.point.effect.stan <- function(y.samples.post,point.pred.mean.post,y.post,
                                  y.model,prob.lower,prob.upper,cycles){
  y.repmat.post <- matrix(y.post, nrow = nrow(y.samples.post),
                          ncol = length(y.post), byrow = TRUE)
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


effect.hsbsts <- function(modeli,alpha=0.05,UnStandardize,y.model,y.cf,pre.time,post.time){
  fitted <- rstan::extract(modeli)
  y.samples <- cbind(fitted$y_fitted,fitted$y_forecast)
  state.samples <- cbind(fitted$mu+fitted$f, fitted$mu_forecast+fitted$f_forecast) 
  prob.lower <- alpha / 2      # e.g., 0.025 when alpha = 0.05
  prob.upper <- 1 - alpha / 2 
  point.pred.mean <-colMeans(state.samples)
  point.pred.lower <- as.numeric(t(apply(y.samples, 2, quantile, prob.lower)))
  point.pred.upper <- as.numeric(t(apply(y.samples, 2, quantile, prob.upper)))
  point.pred <- data.frame(point.pred = point.pred.mean,
                           point.pred.lower, point.pred.upper)
  y.samples <- UnStandardize(y.samples)
  point.pred <- UnStandardize(point.pred)
  
  cum.pred <- ComputeCumulativePredictions(y.samples, point.pred, y.model,
                                           pre.time+1, alpha)
  
  
  y.samples.post <- y.samples[, c((pre.time+1):(pre.time+post.time)), drop = FALSE]
  point.pred.mean.post <- point.pred$point.pred[c((pre.time+1):(pre.time+post.time))]
  y.post <- tail(y.cf, post.time)
  
  summary <- CompileSummaryTable(y.post, y.samples.post, point.pred.mean.post,
                                 alpha)
  return(summary)
}





