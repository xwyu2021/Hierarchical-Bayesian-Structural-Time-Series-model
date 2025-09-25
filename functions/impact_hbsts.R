stan_mod_hbstsLL <- stan_model(file='algorithm (stan)/hbsts_ll.stan',verbose=T)
stan_mod_hbstsAR <- stan_model(file='algorithm (stan)/hbsts_ar.stan',verbose=T)
stan_mod_hbstsLLT <- stan_model(file='algorithm (stan)/hbsts_llt.stan',verbose=T)
stan_mod_hbstsGLLT <- stan_model(file='algorithm (stan)/hbsts_gllt.stan',verbose=T)

run_hbsts <- function(pre.time=pre.time.qua,post.time=post.time.qua,
                      pre.obs=pre.obs.qua,
                      pre.ys=pre.ys.qua,dimx=14,
                      post.obs =post.obs.qua,index_pre=index_pre,
                      K=K,x.obs=x.obs, x.cum.count=x.cum.count,index_post=index_post,
                      niter=4000,nchains=1,time_index=time_index,alpha=0.05,cycles=c(1),
                      emp.mean.y,trend,modelsize=3,burn=1000,adapt=0.92,ind=FALSE,hdi=FALSE,presence_matrix=presence_matrix){
  
  expect.size <- modelsize
  total.size <- dimx
  if(expect.size == dimx){
    tau0=1
  }else{
    tau0 <- expect.size/((total.size-expect.size)*sqrt(pre.time))
  }
  
  if(trend == 'LLT'){
    ls <- list(sigma2_mu=1,mu_err=rep(0,pre.time),sigma2_delta=1,delta_err=rep(0,pre.time))
    lls <- vector("list",1)
    init_fun <- lapply(lls,function(...) ls)
    modeli <- sampling(stan_mod_hbstsLLT,
                       data = list(T=pre.time,T_forecast=post.time,
                                   N=pre.obs,N_forecast=post.obs,
                                   M=dimx,y=pre.ys,
                                   index_obs_t = index_pre,
                                   K=K,x_obs=x.obs,
                                   index_x_t = x.cum.count,index_obs_forecast = index_post,
                                   sigma_itcpt=0.1,
                                   scale_global=tau0, 
                                   nu_global=1,
                                   nu_local=1, 
                                   slab_scale = sqrt(0.2), 
                                   slab_df = 50,c_df=50,tr_df=32,
                                   sdy = 0.01,presence = presence_matrix),
                       iter=niter, chains=nchains, control = list(adapt_delta=adapt),
                       warmup=burn,refresh=burn,init=init_fun)
  }else{
    if(trend == 'AR1'){
      ls <- list(sigma_mu=1,mu_err=rep(0,pre.time),phi=0)
      lls <- vector("list",nchains)
      init_fun <- lapply(lls,function(...) ls)
      modeli <- sampling(stan_mod_hbstsAR,
                         data = list(T=pre.time,T_forecast=post.time,
                                     N=pre.obs,N_forecast=post.obs,
                                     M=dimx,y=pre.ys,
                                     index_obs_t = index_pre,
                                     K=K,x_obs=x.obs,
                                     index_x_t = x.cum.count,index_obs_forecast = index_post,
                                     sigma_itcpt=0.1,
                                     scale_global=tau0, 
                                     nu_global=1,
                                     nu_local=1, 
                                     slab_scale = sqrt(0.2), 
                                     slab_df = 50,c_df=50,tr_df=32,
                                     sdy = 0.01,presence = presence_matrix),
                         iter=niter, chains=nchains, control = list(adapt_delta=adapt),
                         warmup=burn,refresh=burn,init=init_fun)
    }else{
      if(trend == 'LL'){
        ls <- list(sigma2_mu=1,mu_err=rep(0,pre.time))
        lls <- vector("list",1)
        init_fun <- lapply(lls,function(...) ls)
        modeli <- sampling(stan_mod_hbstsLL,
                       data = list(T=pre.time,T_forecast=post.time,
                                   N=pre.obs,N_forecast=post.obs,
                                   M=dimx,y=pre.ys,
                                   index_obs_t = index_pre,
                                   K=K,x_obs=x.obs,
                                   index_x_t = x.cum.count,index_obs_forecast = index_post,
                                   sigma_itcpt=0.1,
                                   scale_global=tau0, 
                                   nu_global=1,
                                   nu_local=1, 
                                   slab_scale = sqrt(0.2), 
                                   slab_df = 50,c_df=50,tr_df=32,
                                   sdy = 0.01,presence = presence_matrix),
                       iter=niter, chains=nchains, control = list(adapt_delta=adapt),
                       warmup=burn,refresh=burn,init=init_fun)
      }else{
        if(trend == 'GLLT'){
          ls <- list(sigma2_mu=1,mu_err=rep(0,pre.time),sigma2_delta=1,delta_err=rep(0,pre.time))
          lls <- vector("list",1)
          init_fun <- lapply(lls,function(...) ls)
          modeli <- sampling(stan_mod_hbstsGLLT,
                             data = list(T=pre.time,T_forecast=post.time,
                                         N=pre.obs,N_forecast=post.obs,
                                         M=dimx,y=pre.ys,
                                         index_obs_t = index_pre,
                                         K=K,x_obs=x.obs,
                                         index_x_t = x.cum.count,index_obs_forecast = index_post,
                                         sigma_itcpt=0.1,
                                         scale_global=tau0, 
                                         nu_global=1,
                                         nu_local=1, 
                                         slab_scale = sqrt(0.2), 
                                         slab_df = 50,c_df=50,tr_df=32,
                                         sdy = 0.01,presence = presence_matrix),
                             iter=niter, chains=nchains, control = list(adapt_delta=adapt),
                             warmup=burn,refresh=burn,init=init_fun)
        }
      }
    }
  }
  
  
  fitted <- rstan::extract(modeli)
  if(ind==TRUE){
    yhat <- fitted$ys
    yhatest <- NULL
    for(i in 1:post.time){
      yhatest <- cbind(yhatest,rowMeans(yhat[,index_post[i]:(index_post[i+1]-1)]))
    }
    y.samples <- cbind(fitted$alpha[,1:pre.time],yhatest)
  }else{
    #y.samples <- cbind(fitted$alpha[,1:pre.time],fitted$alpha_forecast_unscale)
    y.samples <- cbind(fitted$alpha_pred_unscale[,1:pre.time],fitted$alpha_forecast_unscale)
  }
  
  #y.samples <- cbind(exp(fitted$alpha[,1:pre.time] + matrix(fitted$sigma2_obs,nrow=niter-burn,ncol=pre.time)/2),
  #                   exp(fitted$alpha_forecast_unscale + matrix(fitted$sigma2_obs,nrow=niter-burn,ncol=15)/2))
  
  state.samples <- cbind(fitted$mu+fitted$f, fitted$mu_forecast+fitted$f_forecast) 
  state.samples <- matrix(fitted$alpha_sd,nrow=(niter-burn)*nchains,ncol=pre.time+post.time) * state.samples + matrix(fitted$alpha_mean,nrow=(niter-burn)*nchains,ncol=pre.time+post.time)
  prob.lower <- alpha / 2      # e.g., 0.025 when alpha = 0.05
  prob.upper <- 1 - alpha / 2 
  point.pred.mean <- colMeans(y.samples)
  if(hdi==TRUE){
    bound <- hdi(y.samples, credMass = 1-alpha)
    point.pred.lower <- bound[1,]
    point.pred.upper <- bound[2,]
  }else{
    point.pred.lower <- as.numeric(t(apply(y.samples, 2, quantile, prob.lower)))
    point.pred.upper <- as.numeric(t(apply(y.samples, 2, quantile, prob.upper)))
  }
  
  point.pred <- data.frame(point.pred = point.pred.mean,
                           point.pred.lower, point.pred.upper)
  y.model <-emp.mean.y#colMeans(fitted$alpha_bar) ### y.model=origin.mean.y???
  y.cf <- y.model[c((pre.time+1):(pre.time+post.time))]
  if(hdi == TRUE){
    cum.pred <- ComputeCumulativePredictions.hdi(y.samples, point.pred, y.model,
                                             pre.time+1, alpha)
  }else{
    cum.pred <- ComputeCumulativePredictions(y.samples, point.pred, y.model,
                                             pre.time+1, alpha)
  }
  
  y.samples.post <- y.samples[, c((pre.time+1):(pre.time+post.time)), drop = FALSE]
  point.pred.mean.post <- point.pred$point.pred[c((pre.time+1):(pre.time+post.time))]
  y.post <- tail(y.cf, post.time)
  
  if(hdi == TRUE){
    summary <- CompileSummaryTable.hdi(y.post, y.samples.post, point.pred.mean.post,
                                   alpha)
    temporaleffect <- get.point.effect.stan.hdi(y.samples.post,point.pred.mean.post,y.post,
                                            y.model,prob.lower,prob.upper,cycles)
  }else{
    summary <- CompileSummaryTable(y.post, y.samples.post, point.pred.mean.post,
                                   alpha)
    temporaleffect <- get.point.effect.stan(y.samples.post,point.pred.mean.post,y.post,
                                            y.model,prob.lower,prob.upper,cycles)
  }
  
  
  
  cum.y.model <- cumsum.na.rm(y.model)
  series <- zoo(data.frame(y.model, cum.y.model, point.pred, cum.pred),
                time(y.model))
  series$point.effect <- series$y.model - series$point.pred
  series$point.effect.lower <- series$y.model - series$point.pred.upper
  series$point.effect.upper <- series$y.model - series$point.pred.lower
  series$cum.effect <- series$cum.y.model - series$cum.pred
  series$cum.effect.lower <- series$cum.y.model - series$cum.pred.upper
  series$cum.effect.upper <- series$cum.y.model - series$cum.pred.lower
  time(series) <- time_index
  empty <- zoo(, time_index)
  series <- merge(series, empty, all = TRUE)
  # Replace <y.model> by full original response
  series[, 1] <- y.model
  # Assign response-variable names
  names(series)[1] <- "response"
  names(series)[2] <- "cum.response"
  model <- list(pre.period = time_index[c(1,pre.time)],
                post.period = c(tail(time_index,post.time)[1],last(time_index)),
                model.args = NULL,
                bsts.model = modeli,
                alpha = alpha,
                posterior.samples = y.samples)
  impact <- list(series = series,
                 summary = summary,
                 report=NULL,
                 model = model,temporaleffect=temporaleffect)
  class(impact) <- "CausalImpact"
  
  return(impact)
  
  
}  


## -- ADDING RANDOM EFFECT -- ##
stan_mod_hbstsLL_re <- stan_model(file='algorithm (stan)/hbsts_ll_re.stan',verbose=T)
stan_mod_hbstsAR_re <- stan_model(file='algorithm (stan)/hbsts_ar_re.stan',verbose=T)
stan_mod_hbstsLLT_re <- stan_model(file='algorithm (stan)/hbsts_llt_re.stan',verbose=T)
stan_mod_hbstsGLLT_re <- stan_model(file='algorithm (stan)/hbsts_gllt_re.stan',verbose=T)

run_hbsts_re <- function(pre.time=pre.time.qua,post.time=post.time.qua,
                      pre.obs=pre.obs.qua,
                      pre.ys=pre.ys.qua,dimx=14,
                      post.obs =post.obs.qua,index_pre=index_pre,
                      K=K,x.obs=x.obs, x.cum.count=x.cum.count,index_post=index_post,
                      x_strata=x.strata,NS=NS,pre_strata = pre.strata.mon,
                      x_psu=x.psu,NC=NC,pre_psu = pre.psu.mon,
                      niter=4000,nchains=1,time_index=time_index,alpha=0.05,cycles=c(1),
                      emp.mean.y,trend,modelsize=3,burn=1000,adapt=0.92,ind=FALSE,hdi=FALSE,presence_matrix = presence_matrix){
  
  expect.size <- modelsize
  total.size <- dimx
  if(expect.size == dimx){
    tau0=1
  }else{
    tau0 <- expect.size/((total.size-expect.size)*sqrt(pre.time))
  }
  
  if(trend == 'LLT'){
    ls <- list(sigma2_mu=1,mu_err=rep(0,pre.time),sigma2_delta=1,delta_err=rep(0,pre.time))
    lls <- vector("list",1)
    init_fun <- lapply(lls,function(...) ls)
    modeli <- sampling(stan_mod_hbstsLLT_re,
                   data = list(T=pre.time,T_forecast=post.time,
                               N=pre.obs,N_forecast=post.obs,
                               M=dimx,y=pre.ys,
                               index_obs_t = index_pre,
                               K=K,x_obs=x.obs,
                               index_x_t = x.cum.count,index_obs_forecast = index_post,
                               sigma_itcpt=0.1,
                               scale_global=tau0, 
                               nu_global=1,
                               nu_local=1, 
                               slab_scale = 0.1, 
                               slab_df = 50,c_df=50,tr_df=32,
                               sdy = 0.01,x_strata,NS,pre_strata,
                               x_psu,NC,pre_psu,presence = presence_matrix),
                   iter=4000, chains=1, control = list(adapt_delta=adapt),
                   warmup=3000,refresh=1000,init=init_fun)
  }else{
    if(trend == 'AR1'){
      ls <- list(mu_err=rep(0,pre.time),phi=0)
      lls <- vector("list",nchains)
      init_fun <- lapply(lls,function(...) ls)
      modeli <- sampling(stan_mod_hbstsAR_re,
                     data = list(T=pre.time,T_forecast=post.time,
                                 N=pre.obs,N_forecast=post.obs,
                                 M=dimx,y=pre.ys,
                                 index_obs_t = index_pre,
                                 K=K,x_obs=x.obs,
                                 index_x_t = x.cum.count,index_obs_forecast = index_post,
                                 sigma_itcpt=0.1,
                                 scale_global=tau0, 
                                 nu_global=1,
                                 nu_local=1, 
                                 slab_scale = sqrt(0.2), 
                                 slab_df = 50,c_df=50,tr_df=32,
                                 sdy = 0.01,x_strata,NS,pre_strata,
                                 x_psu,NC,pre_psu,presence = presence_matrix),
                     iter=niter, chains=nchains, control = list(adapt_delta=adapt),
                     warmup=burn,refresh=1000,init=init_fun)
    }else{
      if(trend == 'LL'){
        ls <- list(sigma2_mu=1,mu_err=rep(0,pre.time))
        lls <- vector("list",1)
        init_fun <- lapply(lls,function(...) ls)
        modeli <- sampling(stan_mod_hbstsLL_re,
                       data = list(T=pre.time,T_forecast=post.time,
                                   N=pre.obs,N_forecast=post.obs,
                                   M=dimx,y=pre.ys,
                                   index_obs_t = index_pre,
                                   K=K,x_obs=x.obs,
                                   index_x_t = x.cum.count,index_obs_forecast = index_post,
                                   sigma_itcpt=0.1,
                                   scale_global=tau0, 
                                   nu_global=1,
                                   nu_local=1, 
                                   slab_scale = 0.1, 
                                   slab_df = 50,c_df=50,tr_df=32,
                                   sdy = 0.01,x_strata,NS,pre_strata,
                                   x_psu,NC,pre_psu,presence = presence_matrix),
                       iter=4000, chains=1, control = list(adapt_delta=adapt),
                       warmup=3000,refresh=1000,init=init_fun)
      }else{
        if(trend == 'GLLT'){
          ls <- list(sigma2_mu=1,mu_err=rep(0,pre.time),sigma2_delta=1,delta_err=rep(0,pre.time))
          lls <- vector("list",1)
          init_fun <- lapply(lls,function(...) ls)
          modeli <- sampling(stan_mod_hbstsGLLT_re,
                         data = list(T=pre.time,T_forecast=post.time,
                                     N=pre.obs,N_forecast=post.obs,
                                     M=dimx,y=pre.ys,
                                     index_obs_t = index_pre,
                                     K=K,x_obs=x.obs,
                                     index_x_t = x.cum.count,index_obs_forecast = index_post,
                                     sigma_itcpt=0.1,
                                     scale_global=tau0, 
                                     nu_global=1,
                                     nu_local=1, 
                                     slab_scale = 0.1, 
                                     slab_df = 50,c_df=50,tr_df=32,
                                     sdy = 0.01,x_strata,NS,pre_strata,
                                     x_psu,NC,pre_psu,presence = presence_matrix),
                         iter=4000, chains=1, control = list(adapt_delta=adapt),
                         warmup=3000,refresh=1000,init=init_fun)
        }
      }
    }
  }
  
  
  fitted <- rstan::extract(modeli)
  if(ind==TRUE){
    yhat <- fitted$ys
    yhatest <- NULL
    for(i in 1:post.time){
      yhatest <- cbind(yhatest,rowMeans(yhat[,index_post[i]:(index_post[i+1]-1)]))
    }
    y.samples <- cbind(fitted$alpha[,1:pre.time],yhatest)
  }else{
    y.samples <- cbind(fitted$alpha_pred_unscale[,1:pre.time],fitted$alpha_forecast_unscale)
  }
  state.samples <- cbind(fitted$mu+fitted$f, fitted$mu_forecast+fitted$f_forecast) 
  state.samples <- matrix(fitted$alpha_sd,nrow=(niter-burn)*nchains,ncol=pre.time+post.time) * state.samples + matrix(fitted$alpha_mean,nrow=(niter-burn)*nchains,ncol=pre.time+post.time)
  prob.lower <- alpha / 2      # e.g., 0.025 when alpha = 0.05
  prob.upper <- 1 - alpha / 2 
  point.pred.mean <- colMeans(state.samples)
  if(hdi==TRUE){
    bound <- hdi(y.samples, credMass = 1-alpha)
    point.pred.lower <- bound[1,]
    point.pred.upper <- bound[2,]
  }else{
    point.pred.lower <- as.numeric(t(apply(y.samples, 2, quantile, prob.lower)))
    point.pred.upper <- as.numeric(t(apply(y.samples, 2, quantile, prob.upper)))
  }
  
  point.pred <- data.frame(point.pred = point.pred.mean,
                           point.pred.lower, point.pred.upper)
  y.model <-emp.mean.y#colMeans(fitted$alpha_bar) ### y.model=origin.mean.y???
  y.cf <- y.model[c((pre.time+1):(pre.time+post.time))]
  if(hdi == TRUE){
    cum.pred <- ComputeCumulativePredictions.hdi(y.samples, point.pred, y.model,
                                                 pre.time+1, alpha)
  }else{
    cum.pred <- ComputeCumulativePredictions(y.samples, point.pred, y.model,
                                             pre.time+1, alpha)
  }
  
  y.samples.post <- y.samples[, c((pre.time+1):(pre.time+post.time)), drop = FALSE]
  point.pred.mean.post <- point.pred$point.pred[c((pre.time+1):(pre.time+post.time))]
  y.post <- tail(y.cf, post.time)
  
  if(hdi == TRUE){
    summary <- CompileSummaryTable.hdi(y.post, y.samples.post, point.pred.mean.post,
                                       alpha)
    temporaleffect <- get.point.effect.stan.hdi(y.samples.post,point.pred.mean.post,y.post,
                                                y.model,prob.lower,prob.upper,cycles)
  }else{
    summary <- CompileSummaryTable(y.post, y.samples.post, point.pred.mean.post,
                                   alpha)
    temporaleffect <- get.point.effect.stan(y.samples.post,point.pred.mean.post,y.post,
                                            y.model,prob.lower,prob.upper,cycles)
  }
  
  cum.y.model <- cumsum.na.rm(y.model)
  series <- zoo(data.frame(y.model, cum.y.model, point.pred, cum.pred),
                time(y.model))
  series$point.effect <- series$y.model - series$point.pred
  series$point.effect.lower <- series$y.model - series$point.pred.upper
  series$point.effect.upper <- series$y.model - series$point.pred.lower
  series$cum.effect <- series$cum.y.model - series$cum.pred
  series$cum.effect.lower <- series$cum.y.model - series$cum.pred.upper
  series$cum.effect.upper <- series$cum.y.model - series$cum.pred.lower
  time(series) <- time_index
  empty <- zoo(, time_index)
  series <- merge(series, empty, all = TRUE)
  # Replace <y.model> by full original response
  series[, 1] <- y.model
  # Assign response-variable names
  names(series)[1] <- "response"
  names(series)[2] <- "cum.response"
  model <- list(pre.period = time_index[c(1,pre.time)],
                post.period = c(tail(time_index,post.time)[1],last(time_index)),
                model.args = NULL,
                bsts.model = modeli,
                alpha = alpha,
                posterior.samples = y.samples)
  impact <- list(series = series,
                 summary = summary,
                 report=NULL,
                 model = model,temporaleffect=temporaleffect)
  class(impact) <- "CausalImpact"
  
  return(impact)
  
  
} 



effect.hbsts <- function(modeli,pre.time,post.time,niter,burn,nchains=1,alpha=0.05,emp.mean.y,hdi=FALSE){
  fitted <- rstan::extract(modeli)
  y.samples <- cbind(fitted$alpha[,1:pre.time],fitted$alpha_forecast_unscale)
  state.samples <- cbind(fitted$mu+fitted$f, fitted$mu_forecast+fitted$f_forecast) 
  state.samples <- matrix(fitted$alpha_sd,nrow=(niter-burn)*nchains,ncol=pre.time+post.time) * state.samples + matrix(fitted$alpha_mean,nrow=(niter-burn)*nchains,ncol=pre.time+post.time)
  prob.lower <- alpha / 2      # e.g., 0.025 when alpha = 0.05
  prob.upper <- 1 - alpha / 2 
  point.pred.mean <- colMeans(state.samples)
  if(hdi==TRUE){
    bound <- hdi(y.samples, credMass = 1-alpha)
    point.pred.lower <- bound[1,]
    point.pred.upper <- bound[2,]
  }else{
    point.pred.lower <- as.numeric(t(apply(y.samples, 2, quantile, prob.lower)))
    point.pred.upper <- as.numeric(t(apply(y.samples, 2, quantile, prob.upper)))
  }
  
  point.pred <- data.frame(point.pred = point.pred.mean,
                           point.pred.lower, point.pred.upper)
  y.model <-emp.mean.y#colMeans(fitted$alpha_bar) ### y.model=origin.mean.y???
  y.cf <- y.model[c((pre.time+1):(pre.time+post.time))]
  if(hdi == TRUE){
    cum.pred <- ComputeCumulativePredictions.hdi(y.samples, point.pred, y.model,
                                                 pre.time+1, alpha)
  }else{
    cum.pred <- ComputeCumulativePredictions(y.samples, point.pred, y.model,
                                             pre.time+1, alpha)
  }
  
  y.samples.post <- y.samples[, c((pre.time+1):(pre.time+post.time)), drop = FALSE]
  point.pred.mean.post <- point.pred$point.pred[c((pre.time+1):(pre.time+post.time))]
  y.post <- tail(y.cf, post.time)
  
  if(hdi == TRUE){
    summary <- CompileSummaryTable.hdi(y.post, y.samples.post, point.pred.mean.post,
                                       alpha)
  }else{
    summary <- CompileSummaryTable(y.post, y.samples.post, point.pred.mean.post,
                                   alpha)
  }
  return(summary)
}

                             
                             
                             
                             
ComputeCumulativePredictions.hdi <- function(y.samples, point.pred, y,
                                             post.period.begin, alpha = 0.05) {
  
  # Compute posterior mean
  is.post.period <- seq_along(y) >= post.period.begin
  cum.pred.mean.pre <- cumsum.na.rm(as.vector(y)[!is.post.period])
  non.na.indices <- which(!is.na(cum.pred.mean.pre))
  assert_that(length(non.na.indices) > 0)
  last.non.na.index <- max(non.na.indices)
  cum.pred.mean.post <- cumsum(point.pred$point.pred[is.post.period]) +
    cum.pred.mean.pre[last.non.na.index]
  cum.pred.mean <- c(cum.pred.mean.pre, cum.pred.mean.post)
  
  # Check for overflow
  assert_that(identical(which(is.na(cum.pred.mean)),
                        which(is.na(y[!is.post.period]))),
              msg = "unexpected NA found in cum.pred.mean")
  
  # Compute posterior interval
  cum.pred.lower.pre <- cum.pred.mean.pre
  cum.pred.upper.pre <- cum.pred.mean.pre
  y.samples.cum.post <- t(apply(y.samples[, is.post.period, drop = FALSE], 1,
                                cumsum)) +
    cum.pred.mean.pre[last.non.na.index]
  if (sum(is.post.period) == 1) {
    y.samples.cum.post <- t(y.samples.cum.post)
  }
  assert_that(is.scalar(alpha), alpha > 0, alpha < 1)
  prob.lower <- alpha / 2      # e.g., 0.025 when alpha = 0.05
  prob.upper <- 1 - alpha / 2  # e.g., 0.975 when alpha = 0.05
  
  cum.pred.hdi <- apply(y.samples.cum.post,2,hdi,credMass=1-alpha)
  cum.pred.lower.post <- cum.pred.hdi[1,]
  cum.pred.upper.post <- cum.pred.hdi[2,]
  cum.pred.lower <- c(cum.pred.lower.pre, cum.pred.lower.post)
  cum.pred.upper <- c(cum.pred.upper.pre, cum.pred.upper.post)
  
  # Put cumulative prediction together
  cum.pred <- data.frame(cum.pred = cum.pred.mean,
                         cum.pred.lower, cum.pred.upper)
  return(cum.pred)
}

# Tell R CMD check to treat columns of data frames used in `dplyr::mutate` as
# global variables; this avoids false positives of "no visible binding for
# global variable ..." during the check.
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("AbsEffect", "AbsEffect.lower", "AbsEffect.upper",
                           "AbsEffect.sd", "Pred"))
}



CompileSummaryTable.hdi <- function(y.post, y.samples.post,
                                    point.pred.mean.post, alpha = 0.05) {
  # Creates a table of statistics that summarise the post-intervention period.
  # This will later be accessible through \code{impact$model$summary}.
  #
  # Args:
  #   y.post:               Actual observed response during the post-period.
  #   y.samples.post:       Matrix of sampled response trajectories for the
  #                         post-period.
  #   point.pred.mean.post: Posterior predictive mean for the post-period. Note
  #                         that colMeans(y.samples.post) = point.pred.mean.post
  #                         in expectation (i.e., in the limit of an infinite
  #                         number of MCMC iterations); but for any given finite
  #                         simulation, y.samples.post contains sampled
  #                         observation noise. Therefore, to obtain a summary of
  #                         the posterior mean series, we consider the mean of
  #                         the posterior predictive level, without additional
  #                         simulated (centered) observation noise.
  #   alpha:                The resulting coverage of the posterior intervals
  #                         will be \code{1 - alpha}.
  #
  # Returns:
  #   data frame of post-period summary statistics
  
  # Check input
  assert_that(ncol(y.samples.post) == length(y.post),
              msg = "inconsistent y.post")
  assert_that(length(point.pred.mean.post) == length(y.post),
              msg = "inconsistent y.post")
  
  # We will compare the matrix of predicted trajectories (e.g., 900 x 201)
  # with a matrix of replicated observations (e.g., 900 x 201)
  y.repmat.post <- matrix(y.post, nrow = nrow(y.samples.post),
                          ncol = length(y.post), byrow = TRUE)
  assert_that(all(dim(y.repmat.post) == dim(y.samples.post)))
  
  # Define quantiles
  assert_that(is.scalar(alpha), alpha > 0, alpha < 1)
  prob.lower <- alpha / 2      # e.g., 0.025 when alpha = 0.05
  prob.upper <- 1 - alpha / 2  # e.g., 0.975 when alpha = 0.05
  
  # Obtain matrix of posterior samples for relative effect at each time point
  # for post period
  if (length(y.post) == 1) {
    RelEffect.post <- (1 / y.samples.post) * y.post - 1
  } else {
    RelEffect.post <- (1 / y.samples.post) %*% diag(y.post) - 1
  }
  
  # Obtain posterior sample for cumulative relative effect for post period
  RelEffect.cum.post <- sum(y.post)/rowSums(y.samples.post) - 1
  # Obtain posterior sample for average relative effect for post period
  # make average relative effect and cumulative relative effect identical
  RelEffect.average.post <- RelEffect.cum.post
  avg = hdi(rowMeans(y.samples.post),credMass = 1-alpha)
  cum = hdi(rowSums(y.samples.post),credMass = 1-alpha)
  absavg = hdi(rowMeans(y.repmat.post - y.samples.post),credMass = 1-alpha)
  abscum = hdi(rowSums(y.repmat.post - y.samples.post),credMass = 1-alpha)
  relavg = hdi(RelEffect.average.post,credMass = 1-alpha)
  relcum = hdi(RelEffect.cum.post,credMass = 1-alpha)
  # Compile summary statistics
  summary <- data.frame(
    Actual = c(mean(y.post), sum(y.post)),
    Pred = c(mean(point.pred.mean.post), sum(point.pred.mean.post)),
    Pred.lower = c(avg[1],cum[1]),
    Pred.upper = c(avg[2],cum[2]),
    Pred.sd = c(sd(rowMeans(y.samples.post)),
                sd(rowSums(y.samples.post))),
    AbsEffect = c(mean(y.post) - mean(point.pred.mean.post),
                  sum(y.post) - sum(point.pred.mean.post)),
    AbsEffect.lower = c(absavg[1],abscum[1]),
    AbsEffect.upper = c(absavg[2],abscum[2]),
    AbsEffect.sd = c(sd(rowMeans(y.repmat.post - y.samples.post)),
                     sd(rowSums(y.repmat.post - y.samples.post))),
    RelEffect = c(mean(RelEffect.average.post), mean(RelEffect.cum.post)),
    RelEffect.lower = c(relavg[1],relcum[1]),
    RelEffect.upper = c(relavg[2],relcum[2]),
    RelEffect.sd = c(sd(RelEffect.average.post),
                     sd(RelEffect.cum.post))
  )
  rownames(summary) <- c("Average", "Cumulative")
  
  # Add interval coverage, defined by alpha
  summary$alpha <- alpha
  
  # Add one-sided tail-area probability of overall impact, p
  y.samples.post.sum <- rowSums(y.samples.post)
  y.post.sum <- sum(y.post)
  p <- min(sum(c(y.samples.post.sum, y.post.sum) >= y.post.sum),
           sum(c(y.samples.post.sum, y.post.sum) <= y.post.sum)) /
      (length(y.samples.post.sum) + 1)
  assert_that(p > 0, p < 1)
  summary$p <- p
  return(summary)
}



get.point.effect.stan.hdi <- function(y.samples.post,point.pred.mean.post,y.post,
                                      y.model,prob.lower,prob.upper,cycles){
  y.repmat.post <- matrix(y.post, nrow = nrow(y.samples.post),
                          ncol = length(y.post), byrow = TRUE)
  p_vals <- list()
  temp_effect <- list()
  for(k in 1:length(cycles)){
    cycle <- cycles[k]
    p_vals[[k]] <- find_p(y.post,y.samples.post,cycle)
    if(cycle > 1){
      temp_effect[[k]] <- find_effect.hdi(y.post,y.samples.post,point.pred.mean.post,y.repmat.post,cycle,prob.lower,prob.upper)
    }
  }
  return(list(p_vals=p_vals,temp_effect=temp_effect))
}

find_effect.hdi <- function(y.post,y.samples.post,point.pred.mean.post,y.repmat.post,cycle,prob.lower,prob.upper){
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
    temphdi <- hdi(rowMeans(y.repmat.post.sub - y.samples.post.sub),credMass = 2*(1-prob.upper))
    temp.effect.lower[j] <- temphdi[1]
    temp.effect.upper[j] <- temphdi[2]
  }
  return(list(tempeffect = temp.effect,
              tempeffect.lower=temp.effect.lower,
              tempeffect.upper=temp.effect.upper))
}
