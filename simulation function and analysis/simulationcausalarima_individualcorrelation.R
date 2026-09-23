totT <- 100 #100 time points
sd_id <- 4 #individual variation from population mean
min_new <- 50 # min number of new participants at each month
max_new <- 100 # max number of new participants at each month
max_obs_per_id <- 7 # maximum followup times assumed for each person along 100 month
gap_min <- 15 # min followup interval for participants
gap_max <- 24 # max followup interval for participants


sim <- function(M=20,totT=100,sd_mu=0.01,realbeta=c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),correlation=FALSE,seasontrend=FALSE,
                sd_id = 4,min_new = 50,max_new = 100,max_obs_per_id =7,gap_min=15,gap_max=24){
  pred.error <- beta.error  <- inc.error <-
    effect <- cumeffect <- p <- spec <- specavg <- 0
  effect.ci <- cumeffect.ci <- c(0,0)
  sen <- p.sen <- apee <- cumeffectall<-avgeffectall<- senavg <- avgeffectall.ci.lower <- avgeffectall.ci.upper<-cumeffectall.ci.lower <- cumeffectall.ci.upper <- rep(0,7)
  realeffect <- matrix(0,nrow=2,ncol=7)
  newdata <-datagen.arima.ar1(T=totT,dimx=M,sig.x=sig.x,
                              mu.x=rnorm(M,0,0.01),                           beta=realbeta,#c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),#c(rep(1/3,4),0,0),#c(0.2,0.3,0.8,-0.03),#c(0.1,-0.5,1,0.5,rep(0,2)),
                              ar = TRUE,sd_trend = sd_mu,sd_y = 0.5,
                              r = rep(0.9,M),corr=correlation,phi=1)
  
  if(seasontrend == TRUE){
    newdata <- datagen.arima.season(T=totT,dimx = M,sig.x,mu.x=rnorm(M,0,0.1),
                                    beta=realbeta,sd_trend = sd_mu,r = rep(0.99, M),
                                    corr=correlation,phi=1,sd_y=0.5)
  }
  realbeta <- newdata$beta
  
  
  
  
  simemp.mean.y <-  newdata$y
  t.star <- 0.6*totT
  t.post <- t.star+1
  simemp.mean.y1 <- c(simemp.mean.y[1:t.star],simemp.mean.y[t.post:totT] + abs(simemp.mean.y[t.post:totT])*0.01)
  simemp.mean.y10 <- c(simemp.mean.y[1:t.star],simemp.mean.y[t.post:totT] + abs(simemp.mean.y[t.post:totT])*0.03)
  simemp.mean.y30 <- c(simemp.mean.y[1:t.star],simemp.mean.y[t.post:totT]+ abs(simemp.mean.y[t.post:totT])*0.05)
  simemp.mean.y50 <- c(simemp.mean.y[1:t.star],simemp.mean.y[t.post:totT]+ abs(simemp.mean.y[t.post:totT])*0.07)
  simemp.mean.y100 <- c(simemp.mean.y[1:t.star],simemp.mean.y[t.post:totT]+ abs(simemp.mean.y[t.post:totT]*0.1))
  simemp.mean.y200 <- c(simemp.mean.y[1:t.star],simemp.mean.y[t.post:totT]+ abs(simemp.mean.y[t.post:totT])*0.3)
  simemp.mean.y300 <- c(simemp.mean.y[1:t.star],simemp.mean.y[t.post:totT]+ abs(simemp.mean.y[t.post:totT])*0.5)
  
  ci1 <- mean(simemp.mean.y1[t.post:totT] - simemp.mean.y[t.post:totT])
  ci10 <- mean(simemp.mean.y10[t.post:totT] - simemp.mean.y[t.post:totT])
  ci30 <- mean(simemp.mean.y30[t.post:totT] - simemp.mean.y[t.post:totT])
  ci50 <- mean(simemp.mean.y50[t.post:totT] - simemp.mean.y[t.post:totT])
  ci100 <- mean(simemp.mean.y100[t.post:totT] - simemp.mean.y[t.post:totT])
  ci200 <- mean(simemp.mean.y200[t.post:totT] - simemp.mean.y[t.post:totT])
  ci300 <- mean(simemp.mean.y300[t.post:totT] - simemp.mean.y[t.post:totT])
  
  
  realeffect[1,] <- c(ci1,ci10,ci30,ci50,ci100,ci200,ci300)
  realeffect[2,] <- c(sum(simemp.mean.y1[t.post:totT] - simemp.mean.y[t.post:totT]),
                      sum(simemp.mean.y10[t.post:totT] - simemp.mean.y[t.post:totT]),
                      sum(simemp.mean.y30[t.post:totT] - simemp.mean.y[t.post:totT]),
                      sum(simemp.mean.y50[t.post:totT] - simemp.mean.y[t.post:totT]),
                      sum(simemp.mean.y100[t.post:totT] - simemp.mean.y[t.post:totT]),
                      sum(simemp.mean.y200[t.post:totT] - simemp.mean.y[t.post:totT]),
                      sum(simemp.mean.y300[t.post:totT] - simemp.mean.y[t.post:totT]))
  
  
  
  ys <- ty <- id <- c()
  current_id <- 0
  for(entry_t in 1:totT){
    # number of new participants entering this month
    n_new <- sample(min_new:max_new, 1)
    for(k in 1:n_new){
      current_id <- current_id + 1
      
      # individual random intercept
      b_i <- rnorm(1, mean = 0, sd = sd_id)
      
      # intended maximum number of observations for this person
      m_i <- sample(1:max_obs_per_id, 1)
      
      # generate observation times
      obs_times <- entry_t
      
      if (m_i > 1) {
        current_t <- entry_t
        
        for (jj in 2:m_i) {
          next_gap <- sample(gap_min:gap_max, 1)
          next_t <- current_t + next_gap
          
          if (next_t > totT) break
          
          obs_times <- c(obs_times, next_t)
          current_t <- next_t
        }
      }
      
      # generate y only at observed months for this person
      y_i <- sapply(obs_times, function(tt) {
        rtruncnorm(1,a = 0,b = 36,mean = simemp.mean.y[tt] + b_i,sd = 6)})
      
      ys <- c(ys, y_i)
      ty <- c(ty, obs_times)
      id <- c(id, rep(current_id, length(obs_times)))
    }
  }
  
  simy <- data.frame(y=ys,t=ty)
  
  simy_agg <- aggregate(y ~t,data=simy, mean)
  simx <- newdata$X[1:totT,-1]
  ind <- xs <- t <- c()
  sdx <- rtruncnorm(M,a=0,b=10,6.5,3)
  
  for(i in 1:M){
    for(entry_t in 1:totT){
      n_new <- sample(min_new:max_new, 1)
      for(k in 1:n_new){
        current_id <- current_id + 1
        
        b_i <- rnorm(1, mean = 0, sd = sd_id)
        
        m_i <- sample(1:max_obs_per_id, 1)
        
        obs_times <- entry_t
        
        if(m_i > 1){
          current_t <- entry_t
          
          for(jj in 2:m_i){
            next_gap <- sample(gap_min:gap_max, 1)
            next_t <- current_t + next_gap
            
            if(next_t > totT){break}
            
            obs_times <- c(obs_times, next_t)
            current_t <- next_t
          }
          
        }
        
        # generate y only at observed months
        x_i <- sapply(obs_times, function(tt) {
          rtruncnorm(1,a = 0,b = 36,mean = simx[tt,i] + b_i,sd = sdx[i])})
        xs <- c(xs,x_i)
        t <- c(t,obs_times)
        ind <- c(ind,rep(i,length(obs_times)))
      }
    }
  }
  
  
  simxobs <- data.frame(x=xs,t=t,ind=ind)
  
  simxobs <- simxobs %>% dplyr::arrange(t,ind)
  simx_agg <- aggregate(x ~ t+ind,data=simxobs,mean)
  simx_agg_df <- dcast(simx_agg,t~ind)
  simx_agg_df <- simx_agg_df[-1]
  
  
  simpre.time <- t.star
  simpost.time <- totT-t.star
  simpre.ys <- as.vector(simy$y[simy$t <=t.star])
  simpost.ys <- as.vector(simy$y[simy$t >t.star])
  
  simpre.obs <- sum(simy$t <=t.star)
  simpost.obs <- sum(simy$t >t.star)
  
  simobs.count <- as.vector(table(simy$t))
  simindex_pre <- cumsum(c(1,simobs.count[1:simpre.time]))
  
  simindex_post <- cumsum(c(1,simobs.count[(simpre.time+1):(simpre.time+simpost.time)]))
  
  checksim <- as.data.frame(simxobs[,2:3])
  simx.count <- table(checksim)
  simx.cum.count <- t(matrix(cumsum(t(simx.count)),nrow=M,ncol=totT))
  simxobs <- as.vector(simxobs[,1])
  #simxobs <- as.vector(stand.simxobs[,1])
  simK <- length(simxobs)
  
  
  calendar.time <- as.Date(1:totT)
  calendar.intervene <- calendar.time[t.post]
  
  cycle <- ifelse(seasontrend == TRUE,4,1)
  
  safe_impact <- function(model, digits = 3) {
    assign("ce", model, envir = .GlobalEnv)
    impact(ce, digits = digits)
  }
  camod <- CausalArima(y=ts(simy_agg$y,frequency = cycle),dates = calendar.time, int.date = calendar.intervene,xreg =as.matrix(simx_agg_df), nboot = 1000)
  camod_summary <- safe_impact(camod)
  
  
  pred.error <- RMSE(camod$forecast,simemp.mean.y[t.post:totT])
  if(length(camod_summary$arima$param[,1])!=M){
    beta.error <- NA
  }else{
    beta.error <- RMSE(camod_summary$arima$param[,1],realbeta[2:(M+1)]) 
  }
  
  effect <- camod_summary$impact_boot$average$estimates[3]
  effect.ci <- c(camod_summary$impact_boot$average$inf[3],camod_summary$impact_boot$average$sup[3])
  cumeffect <- camod_summary$impact_boot$effect_cum$estimates[3]
  cumeffect.ci <- c(camod_summary$impact_boot$effect_cum$inf[3],camod_summary$impact_boot$effect_cum$sup[3])
  p <- camod_summary$impact_boot$p_values[2]
  spec <- as.numeric((0>= c(camod_summary$impact_boot$effect_cum$inf[3]) & (0 <= camod_summary$impact_boot$effect_cum$sup[3])))
  specavg <- as.numeric((0>= c(camod_summary$impact_boot$average$inf[3]) & (0 <= camod_summary$impact_boot$average$sup[3])))
  
  ys <- ty <- c()
  current_id <- 0
  for(entry_t in 1:totT){
    
    n_new <- sample(min_new:max_new, 1)
    for(k in 1:n_new){
      current_id <- current_id + 1
      
      b_i <- rnorm(1, mean = 0, sd = sd_id)
      
      m_i <- sample(1:max_obs_per_id, 1)
      
      obs_times <- entry_t
      
      if (m_i > 1) {
        current_t <- entry_t
        
        for (jj in 2:m_i) {
          next_gap <- sample(gap_min:gap_max, 1)
          next_t <- current_t + next_gap
          
          if (next_t > totT) break
          
          obs_times <- c(obs_times, next_t)
          current_t <- next_t
        }
      }
      
      # generate y only at observed months for this person
      y_i <- sapply(obs_times, function(tt) {
        rtruncnorm(1,a = 0,b = 36,mean = simemp.mean.y1[tt] + b_i,sd = 6)})
      
      ys <- c(ys, y_i)
      ty <- c(ty, obs_times)
    }
  }
  
  simy <- data.frame(y=ys,t=ty)
  simy_agg <- aggregate(y ~t,data=simy, mean)
  simy1 <- simy_agg$y
  
  ys <- ty <- c()
  current_id <- 0
  for(entry_t in 1:totT){
    
    n_new <- sample(min_new:max_new, 1)
    for(k in 1:n_new){
      current_id <- current_id + 1
      
      b_i <- rnorm(1, mean = 0, sd = sd_id)
      
      m_i <- sample(1:max_obs_per_id, 1)
      
      obs_times <- entry_t
      
      if (m_i > 1) {
        current_t <- entry_t
        
        for (jj in 2:m_i) {
          next_gap <- sample(gap_min:gap_max, 1)
          next_t <- current_t + next_gap
          
          if (next_t > totT) break
          
          obs_times <- c(obs_times, next_t)
          current_t <- next_t
        }
      }
      
      # generate y only at observed months for this person
      y_i <- sapply(obs_times, function(tt) {
        rtruncnorm(1,a = 0,b = 36,mean = simemp.mean.y10[tt] + b_i,sd = 6)})
      
      ys <- c(ys, y_i)
      ty <- c(ty, obs_times)
    }
  }
  
  simy <- data.frame(y=ys,t=ty)
  simy_agg <- aggregate(y ~t,data=simy, mean)
  simy10 <- simy_agg$y
  
  
  ys <- ty <- c()
  current_id <- 0
  for(entry_t in 1:totT){
    
    n_new <- sample(min_new:max_new, 1)
    for(k in 1:n_new){
      current_id <- current_id + 1
      
      b_i <- rnorm(1, mean = 0, sd = sd_id)
      
      m_i <- sample(1:max_obs_per_id, 1)
      
      obs_times <- entry_t
      
      if (m_i > 1) {
        current_t <- entry_t
        
        for (jj in 2:m_i) {
          next_gap <- sample(gap_min:gap_max, 1)
          next_t <- current_t + next_gap
          
          if (next_t > totT) break
          
          obs_times <- c(obs_times, next_t)
          current_t <- next_t
        }
      }
      
      # generate y only at observed months for this person
      y_i <- sapply(obs_times, function(tt) {
        rtruncnorm(1,a = 0,b = 36,mean = simemp.mean.y30[tt] + b_i,sd = 6)})
      
      ys <- c(ys, y_i)
      ty <- c(ty, obs_times)
      
    }
  }
  
  simy <- data.frame(y=ys,t=ty)
  simy_agg <- aggregate(y ~t,data=simy, mean)
  simy30 <- simy_agg$y
  
  
  ys <- ty <- c()
  current_id <- 0
  for(entry_t in 1:totT){
    
    n_new <- sample(min_new:max_new, 1)
    for(k in 1:n_new){
      current_id <- current_id + 1
      
      b_i <- rnorm(1, mean = 0, sd = sd_id)
      
      m_i <- sample(1:max_obs_per_id, 1)
      
      obs_times <- entry_t
      
      if (m_i > 1) {
        current_t <- entry_t
        
        for (jj in 2:m_i) {
          next_gap <- sample(gap_min:gap_max, 1)
          next_t <- current_t + next_gap
          
          if (next_t > totT) break
          
          obs_times <- c(obs_times, next_t)
          current_t <- next_t
        }
      }
      
      # generate y only at observed months for this person
      y_i <- sapply(obs_times, function(tt) {
        rtruncnorm(1,a = 0,b = 36,mean = simemp.mean.y50[tt] + b_i,sd = 6)})
      
      ys <- c(ys, y_i)
      ty <- c(ty, obs_times)
      
    }
  }
  
  simy <- data.frame(y=ys,t=ty)
  simy_agg <- aggregate(y ~t,data=simy, mean)
  simy50 <- simy_agg$y
  
  ys <- ty <- c()
  current_id <- 0
  for(entry_t in 1:totT){
    
    n_new <- sample(min_new:max_new, 1)
    for(k in 1:n_new){
      current_id <- current_id + 1
      
      b_i <- rnorm(1, mean = 0, sd = sd_id)
      
      m_i <- sample(1:max_obs_per_id, 1)
      
      obs_times <- entry_t
      
      if (m_i > 1) {
        current_t <- entry_t
        
        for (jj in 2:m_i) {
          next_gap <- sample(gap_min:gap_max, 1)
          next_t <- current_t + next_gap
          
          if (next_t > totT) break
          
          obs_times <- c(obs_times, next_t)
          current_t <- next_t
        }
      }
      
      # generate y only at observed months for this person
      y_i <- sapply(obs_times, function(tt) {
        rtruncnorm(1,a = 0,b = 36,mean = simemp.mean.y100[tt] + b_i,sd = 6)})
      
      ys <- c(ys, y_i)
      ty <- c(ty, obs_times)
      
    }
  }
  
  simy <- data.frame(y=ys,t=ty)
  simy_agg <- aggregate(y ~t,data=simy, mean)
  simy100 <- simy_agg$y
  
  
  
  ys <- ty <- c()
  current_id <- 0
  for(entry_t in 1:totT){
    
    n_new <- sample(min_new:max_new, 1)
    for(k in 1:n_new){
      current_id <- current_id + 1
      
      b_i <- rnorm(1, mean = 0, sd = sd_id)
      
      m_i <- sample(1:max_obs_per_id, 1)
      
      obs_times <- entry_t
      
      if (m_i > 1) {
        current_t <- entry_t
        
        for (jj in 2:m_i) {
          next_gap <- sample(gap_min:gap_max, 1)
          next_t <- current_t + next_gap
          
          if (next_t > totT) break
          
          obs_times <- c(obs_times, next_t)
          current_t <- next_t
        }
      }
      
      # generate y only at observed months for this person
      y_i <- sapply(obs_times, function(tt) {
        rtruncnorm(1,a = 0,b = 36,mean = simemp.mean.y200[tt] + b_i,sd = 6)})
      
      ys <- c(ys, y_i)
      ty <- c(ty, obs_times)
      
    }
  }
  
  simy <- data.frame(y=ys,t=ty)
  simy_agg <- aggregate(y ~t,data=simy, mean)
  simy200 <- simy_agg$y
  
  
  ys <- ty <- c()
  current_id <- 0
  for(entry_t in 1:totT){
    
    n_new <- sample(min_new:max_new, 1)
    for(k in 1:n_new){
      current_id <- current_id + 1
      
      b_i <- rnorm(1, mean = 0, sd = sd_id)
      
      m_i <- sample(1:max_obs_per_id, 1)
      
      obs_times <- entry_t
      
      if (m_i > 1) {
        current_t <- entry_t
        
        for (jj in 2:m_i) {
          next_gap <- sample(gap_min:gap_max, 1)
          next_t <- current_t + next_gap
          
          if (next_t > totT) break
          
          obs_times <- c(obs_times, next_t)
          current_t <- next_t
        }
      }
      
      # generate y only at observed months for this person
      y_i <- sapply(obs_times, function(tt) {
        rtruncnorm(1,a = 0,b = 36,mean = simemp.mean.y300[tt] + b_i,sd = 6)})
      
      ys <- c(ys, y_i)
      ty <- c(ty, obs_times)
      
    }
  }
  
  simy <- data.frame(y=ys,t=ty)
  simy_agg <- aggregate(y ~t,data=simy, mean)
  simy300 <- simy_agg$y
  
  
  camod1 <- CausalArima(y=ts(c(simy_agg$y[1:t.star],simy1[t.post:totT]),frequency = cycle),dates = calendar.time, int.date = calendar.intervene,xreg =as.matrix(simx_agg_df), nboot = 1000)
  camod_summary1 <- safe_impact(camod1)
  camod10 <- CausalArima(y=ts(c(simy_agg$y[1:t.star],simy10[t.post:totT]),frequency = cycle),dates = calendar.time, int.date = calendar.intervene,xreg =as.matrix(simx_agg_df), nboot = 1000)
  camod_summary10 <- safe_impact(camod10)
  camod30 <- CausalArima(y=ts(c(simy_agg$y[1:t.star],simy30[t.post:totT]),frequency = cycle),dates = calendar.time, int.date = calendar.intervene,xreg =as.matrix(simx_agg_df), nboot = 1000)
  camod_summary30 <- safe_impact(camod30)
  camod50 <- CausalArima(y=ts(c(simy_agg$y[1:t.star],simy50[t.post:totT]),frequency = cycle),dates = calendar.time, int.date = calendar.intervene,xreg =as.matrix(simx_agg_df), nboot = 1000)
  camod_summary50 <- safe_impact(camod50)
  camod100 <- CausalArima(y=ts(c(simy_agg$y[1:t.star],simy100[t.post:totT]),frequency = cycle),dates = calendar.time, int.date = calendar.intervene,xreg =as.matrix(simx_agg_df), nboot = 1000)
  camod_summary100 <- safe_impact(camod100)
  camod200 <- CausalArima(y=ts(c(simy_agg$y[1:t.star],simy200[t.post:totT]),frequency = cycle),dates = calendar.time, int.date = calendar.intervene,xreg =as.matrix(simx_agg_df), nboot = 1000)
  camod_summary200 <- safe_impact(camod200)
  camod300 <- CausalArima(y=ts(c(simy_agg$y[1:t.star],simy300[t.post:totT]),frequency = cycle),dates = calendar.time, int.date = calendar.intervene,xreg =as.matrix(simx_agg_df), nboot = 1000)
  camod_summary300 <- safe_impact(camod300)
  
  
  avgeffectall[1] <- camod_summary1$impact_boot$average$estimates[3]
  avgeffectall[2] <- camod_summary10$impact_boot$average$estimates[3]
  avgeffectall[3] <- camod_summary30$impact_boot$average$estimates[3]
  avgeffectall[4] <- camod_summary50$impact_boot$average$estimates[3]
  avgeffectall[5] <- camod_summary100$impact_boot$average$estimates[3]
  avgeffectall[6] <- camod_summary200$impact_boot$average$estimates[3]
  avgeffectall[7] <- camod_summary300$impact_boot$average$estimates[3]
  
  cumeffectall[1] <- camod_summary1$impact_boot$effect_cum$estimates[3]
  cumeffectall[2] <- camod_summary10$impact_boot$effect_cum$estimates[3]
  cumeffectall[3] <- camod_summary30$impact_boot$effect_cum$estimates[3]
  cumeffectall[4] <- camod_summary50$impact_boot$effect_cum$estimates[3]
  cumeffectall[5] <- camod_summary100$impact_boot$effect_cum$estimates[3]
  cumeffectall[6] <- camod_summary200$impact_boot$effect_cum$estimates[3]
  cumeffectall[7] <- camod_summary300$impact_boot$effect_cum$estimates[3]
  
  
  avgeffectall.ci.lower[1] <- camod_summary1$impact_boot$average$inf[3]
  avgeffectall.ci.lower[2] <- camod_summary10$impact_boot$average$inf[3]
  avgeffectall.ci.lower[3] <- camod_summary30$impact_boot$average$inf[3]
  avgeffectall.ci.lower[4] <- camod_summary50$impact_boot$average$inf[3]
  avgeffectall.ci.lower[5] <- camod_summary100$impact_boot$average$inf[3]
  avgeffectall.ci.lower[6] <- camod_summary200$impact_boot$average$inf[3]
  avgeffectall.ci.lower[7] <- camod_summary300$impact_boot$average$inf[3]
  
  avgeffectall.ci.upper[1] <- camod_summary1$impact_boot$average$sup[3]
  avgeffectall.ci.upper[2] <- camod_summary10$impact_boot$average$sup[3]
  avgeffectall.ci.upper[3] <- camod_summary30$impact_boot$average$sup[3]
  avgeffectall.ci.upper[4] <- camod_summary50$impact_boot$average$sup[3]
  avgeffectall.ci.upper[5] <- camod_summary100$impact_boot$average$sup[3]
  avgeffectall.ci.upper[6] <- camod_summary200$impact_boot$average$sup[3]
  avgeffectall.ci.upper[7] <- camod_summary300$impact_boot$average$sup[3]
  
  cumeffectall.ci.lower[1] <- camod_summary1$impact_boot$effect_cum$inf[3]
  cumeffectall.ci.lower[2] <- camod_summary10$impact_boot$effect_cum$inf[3]
  cumeffectall.ci.lower[3] <- camod_summary30$impact_boot$effect_cum$inf[3]
  cumeffectall.ci.lower[4] <- camod_summary50$impact_boot$effect_cum$inf[3]
  cumeffectall.ci.lower[5] <- camod_summary100$impact_boot$effect_cum$inf[3]
  cumeffectall.ci.lower[6] <- camod_summary200$impact_boot$effect_cum$inf[3]
  cumeffectall.ci.lower[7] <- camod_summary300$impact_boot$effect_cum$inf[3]
  
  cumeffectall.ci.upper[1] <- camod_summary1$impact_boot$effect_cum$sup[3]
  cumeffectall.ci.upper[2] <- camod_summary10$impact_boot$effect_cum$sup[3]
  cumeffectall.ci.upper[3] <- camod_summary30$impact_boot$effect_cum$sup[3]
  cumeffectall.ci.upper[4] <- camod_summary50$impact_boot$effect_cum$sup[3]
  cumeffectall.ci.upper[5] <- camod_summary100$impact_boot$effect_cum$sup[3]
  cumeffectall.ci.upper[6] <- camod_summary200$impact_boot$effect_cum$sup[3]
  cumeffectall.ci.upper[7] <- camod_summary300$impact_boot$effect_cum$sup[3]
  
  
  p.sen[1] <- camod_summary1$impact_boot$p_values[2]
  p.sen[2] <- camod_summary10$impact_boot$p_values[2]
  p.sen[3] <- camod_summary30$impact_boot$p_values[2]
  p.sen[4] <- camod_summary50$impact_boot$p_values[2]
  p.sen[5] <- camod_summary100$impact_boot$p_values[2]
  p.sen[6] <- camod_summary200$impact_boot$p_values[2]
  p.sen[7] <- camod_summary300$impact_boot$p_values[2]
  
  apee[1] <- apee_est(camod_summary1$impact_boot$average$estimates[3],ci1)
  apee[2] <- apee_est(camod_summary10$impact_boot$average$estimates[3],ci10)
  apee[3] <- apee_est(camod_summary30$impact_boot$average$estimates[3],ci30)
  apee[4] <- apee_est(camod_summary50$impact_boot$average$estimates[3],ci50)
  apee[5] <- apee_est(camod_summary100$impact_boot$average$estimates[3],ci100)
  apee[6] <- apee_est(camod_summary200$impact_boot$average$estimates[3],ci200)
  apee[7] <- apee_est(camod_summary300$impact_boot$average$estimates[3],ci300)
  
  
  senavg[1] <- as.numeric((0 <= camod_summary1$impact_boot$average$inf[3])|(0>=camod_summary1$impact_boot$average$sup[3]))
  senavg[2] <- as.numeric((0 <= camod_summary10$impact_boot$average$inf[3])|(0>=camod_summary10$impact_boot$average$sup[3]))
  senavg[3] <- as.numeric((0 <= camod_summary30$impact_boot$average$inf[3])|(0>=camod_summary30$impact_boot$average$sup[3]))
  senavg[4] <- as.numeric((0 <= camod_summary50$impact_boot$average$inf[3])|(0>=camod_summary50$impact_boot$average$sup[3]))
  senavg[5] <- as.numeric((0 <= camod_summary100$impact_boot$average$inf[3])|(0>=camod_summary100$impact_boot$average$sup[3]))
  senavg[6] <- as.numeric((0 <= camod_summary200$impact_boot$average$inf[3])|(0>=camod_summary200$impact_boot$average$sup[3]))
  senavg[7] <- as.numeric((0 <= camod_summary300$impact_boot$average$inf[3])|(0>=camod_summary300$impact_boot$average$sup[3]))
  
  sen[1] <- as.numeric((0 <= c(camod_summary1$impact_boot$effect_cum$inf[3])|(0 >= camod_summary1$impact_boot$effect_cum$sup[3])))
  sen[2] <- as.numeric((0 <= c(camod_summary10$impact_boot$effect_cum$inf[3])|(0 >= camod_summary10$impact_boot$effect_cum$sup[3])))
  sen[3] <- as.numeric((0 <= c(camod_summary30$impact_boot$effect_cum$inf[3])|(0 >= camod_summary30$impact_boot$effect_cum$sup[3])))
  sen[4] <- as.numeric((0 <= c(camod_summary50$impact_boot$effect_cum$inf[3])|(0 >= camod_summary50$impact_boot$effect_cum$sup[3])))
  sen[5] <- as.numeric((0 <= c(camod_summary100$impact_boot$effect_cum$inf[3])|(0 >= camod_summary100$impact_boot$effect_cum$sup[3])))
  sen[6] <- as.numeric((0 <= c(camod_summary200$impact_boot$effect_cum$inf[3])|(0 >= camod_summary200$impact_boot$effect_cum$sup[3])))
  sen[7] <- as.numeric((0 <= c(camod_summary300$impact_boot$effect_cum$inf[3])|(0 >= camod_summary300$impact_boot$effect_cum$sup[3])))
  
  
  
  return(list(pred.error = pred.error,beta.error=beta.error,inc.error=inc.error,p = p,p.sen=p.sen,apee=apee,
              effect=effect, spec=spec, sen=sen,cumeffect=cumeffect,
              specavg=specavg,cumeffectall=cumeffectall,avgeffectall=avgeffectall,senavg=senavg,
              effect.ci=effect.ci,cumeffect.ci=cumeffect.ci,
              avgeffectall.ci.lower=avgeffectall.ci.lower,avgeffectall.ci.upper=avgeffectall.ci.upper,
              cumeffectall.ci.lower=cumeffectall.ci.lower,cumeffectall.ci.upper=cumeffectall.ci.upper,
              realeffect=realeffect))
  
}  



## this function merge the paralleled outputs

clean <- function(res){
  out <- list(
    pred.error = sapply(res, `[[`, "pred.error"),
    beta.error = sapply(res, `[[`, "beta.error"),
    inc.error  = sapply(res, `[[`, "inc.error"),
    p          = sapply(res, `[[`, "p"),
    effect     = sapply(res, `[[`, "effect"),
    spec       = sapply(res, `[[`, "spec"),
    cumeffect  = sapply(res, `[[`, "cumeffect"),
    specavg    = sapply(res, `[[`, "specavg")
  )  
  out$p.sen <- do.call(rbind, lapply(res, `[[`, "p.sen"))
  out$apee <- do.call(rbind, lapply(res, `[[`, "apee"))
  out$sen <- do.call(rbind, lapply(res, `[[`, "sen"))
  out$cumeffectall <- do.call(rbind, lapply(res, `[[`, "cumeffectall"))
  out$avgeffectall <- do.call(rbind, lapply(res, `[[`, "avgeffectall"))
  out$senavg <- do.call(rbind, lapply(res, `[[`, "senavg"))
  
  out$avgeffectall.ci.lower <- do.call(rbind, lapply(res, `[[`, "avgeffectall.ci.lower"))
  out$avgeffectall.ci.upper <- do.call(rbind, lapply(res, `[[`, "avgeffectall.ci.upper"))
  out$cumeffectall.ci.lower <- do.call(rbind, lapply(res, `[[`, "cumeffectall.ci.lower"))
  out$cumeffectall.ci.upper <- do.call(rbind, lapply(res, `[[`, "cumeffectall.ci.upper"))
  
  out$realeffect <- array(
    unlist(lapply(res, `[[`, "realeffect")),
    dim = c(2, 7, length(res))
  )
  return(out)
}

B <- 200
ncores <- detectCores() - 1
set.seed(123)


res <- mclapply(1:10,function(b) {sim(M=20,totT=100,sd_mu=0.01,realbeta=c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),correlation=FALSE,seasontrend=FALSE,
                                     sd_id = 4,min_new = 50,max_new = 100,max_obs_per_id =7,gap_min=15,gap_max=24)},
                mc.cores = ncores
)

 ca_individualcorr <- clean(res)
saveRDS(ca_individualcorr,"ca_individualcorr.RData")


res1 <- mclapply(1:200,function(b) {sim(M=5,totT=100,sd_mu=0.01,realbeta=c(rep(1/3,4),0,0),correlation=FALSE,seasontrend=FALSE,
                                     sd_id = 4,min_new = 50,max_new = 100,max_obs_per_id =7,gap_min=15,gap_max=24)},
                mc.cores = ncores
)

ca_individualcorr_5 <- clean(res1)
saveRDS(ca_individualcorr_5,"ca_individualcorr_5.RData")

