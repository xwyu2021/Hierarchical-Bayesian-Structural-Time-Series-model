library(parallel)
library(lubridate)
library(MASS)
library(dplyr)
library(gplots)
library(reshape2)
library(bsts)
library(CausalImpact)
library(assertthat)
library(truncnorm)
library(tmvtnorm)
library(tidyr)
library('rstan')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(HDInterval)


#source("/root/autodl-tmp/lse project/impact_additional.R")
#source("/root/autodl-tmp/lse project/impact_stan.R")
#source("/root/autodl-tmp/lse project/impact_hbsts.R")
#source("/root/autodl-tmp/lse project/metrics.R")
#source("/root/autodl-tmp/lse project/impact_hbsts_stand.R")
#source("/root/autodl-tmp/lse project/impact_model.R")
#source("/root/autodl-tmp/lse project/impact_misc.R")
#source("/root/autodl-tmp/lse project/impact_analysis.R")
#source("/root/autodl-tmp/lse project/impact_inference.R")

####################################################
### absolute percentatge error estimation (APEE) ###
####################################################
apee_est <- function(est_ef,real_ef){
  return(abs(est_ef - real_ef)/real_ef)
}

###################################
### construct covariance matrix ###
###################################
qij <-function(corrij=0.7,ri,rj){
  # corrij: correlation between ith and jth times series
  # ri: ar(1) autocorrelation coefficient
  return(corrij * (1-ri*rj)/sqrt((1-ri^2)*(1-rj^2)))
}

#################################################
### data generating function with seasonality ###
#################################################
datagen.arima.season <- function(T=100,dimx = 5,sig.x,mu.x=rep(1,5),beta,sd_trend = 0.1,
                              phi=0.5,sd_y=0.1,corr =FALSE,r){
  trend=rep(0,T) #Trend component
  X <- matrix(0,nrow=T,ncol=dimx)
  if(corr == TRUE){
    #Q <- matrix(0,dimx,dimx)
    #for(i in 1:dimx){
    #  for(j in 1:dimx){
    #    if(i != j){
    #      Q[i,j] <- qij(0.7,r[i],r[j])
    #    }else{
    #      Q[i,j] <- qij(1,r[i],r[j])
    #    }
    #  }
    #}
    eps <- mvrnorm(T,mu = rep(0,dimx),Sigma=sig.x)
    eps[1,] <-  mvrnorm(1,mu = mu.x,Sigma=sig.x)   # desired initial value
  }else{
    eps <- matrix(0,nrow=T,ncol=dimx)
    for(i in 1:T){
      for(j in 1:dimx){
        eps[i,j] <- rnorm(1,0,1)
      }
    }
    eps[1,] <- rnorm(5,0,1)
  }
  # (1): simulate control series,inc. intercept
  
  #for(i in 1:(dimx)){
  #  X[,i] <- arima.sim(list(order = c(1,0,0),ar=r[i]),T,innov=eps[,i])
  #  X[,i] <- (X[,i] - min(X[,i]))+0.01
  #}
  xs <- runif(dimx,0.1,4)#c(7,2,0.1,5,8) #c(2.5,2,0.1,1,0.5)  
     #xs <- c(0.2,1.5,0.8) 
    xm <- rnorm(dimx,11,0.05)  
    for(i in 1:dimx){
        X[,i] <- rnorm(T,11,xs[i])
    } 
  Xtilde <- X  
  for(i in 1:dimx){
      Xtilde[,i] <- (Xtilde[,i] - mean(Xtilde[,i]))/sd(Xtilde[,i])
  }
  Xtilde <- cbind(rep(1,T),Xtilde) 
  X <- cbind(rep(1,T),X)
  
  # (2): simulate trend
  trend=rep(0,T) #Trend component
  trend[1] <- rnorm(1,mean=0,sd=sd_trend)
  
    # SET AR(1) model FOR TREND
    for(i in 2:T){
      trend[i] = phi * trend[i-1] + rnorm(n=1,mean=0,sd=sd_trend)
    }
  
  # (3): simulate seasonal component
  
  gamma <- trend
  gamma[1] <- rnorm(1,mean=0,sd=0.5)
  gamma[2] <- rnorm(1,mean=0,sd=0.5)
  gamma[3] <- rnorm(1,mean=0,sd=0.5)
  for(i in 4:T){
    gamma[i] = -(gamma[i-1] + gamma[i-2] + gamma[i-3]) + rnorm(n=1,mean=0,sd = 0.5)
  }
  
  #y <- X%*%beta + rnorm(T,0,0.1)
  y <- Xtilde%*%beta + trend + rnorm(T,0,sd_y)
  y <- y*1.5 + 11
  out <- list(y=y,X=X,beta=beta)
  return(out)
}

sig.x <- matrix(c(1,0.8,0.85,0.75,0.82,
                  0.8,1,0.78,0.83,0.79,
                  0.85,0.78,1,0.81,0.77,
                  0.75,0.83,0.81,1,0.84,
                  0.82,0.79,0.77,0.84,1),nrow=5,ncol=5)

###################################
### data generating function ######
###################################
datagen.arima.ar1 <- function(T=200,dimx = 5,sig.x=sig.x,mu.x=mu.x,beta,
                              ar = TRUE,rho = 0,Dtilde = 1,
                              sd_trend = 0.1,
                              phi=0.5,sd_y=0.1,r,corr =FALSE){
  trend=rep(0,T) #Trend component
  X <- matrix(0,nrow=T,ncol=dimx)
  if(corr == TRUE){
    R <- matrix(0.85, nrow = dimx, ncol = dimx)
    diag(R) <- 1
   xs <- runif(dimx,0.1,4)  
  Sgm <- diag(xs) %*% R %*% diag(xs)
  X <- mvrnorm(n = T, mu = rnorm(dimx,11,0.05), Sigma = Sgm)  
    
  }else{
    #eps <- matrix(0,nrow=T,ncol=dimx)
    #for(i in 1:T){
    #  for(j in 1:dimx){
    #    eps[i,j] <- rnorm(1,0,1)
    #  }
    #} 
    #X[1,] <- rnorm(dimx,11,1)
    #for(t in 2:T){
    #  for(i in 1:dimx){
    #    X[t,i] <- 11*(1-r[i]) + r[i] * X[t-1,i] + rnorm(1,0,1)
    #  }
    #}
      #xs <- rep(c(3,2,0.1,4,5),each=4)#c(5,3,0.1,5,6) # c(7,2,0.1,5,8)
      xs <- runif(dimx,0.1,4)
     for(i in 1:dimx){
        X[,i] <- rnorm(T,11,xs[i])
    } 
  }
  # (1): simulate control series,inc. intercept
  #col_means <- colMeans(X)
  #col_sds <- apply(X, 2, sd)
  #beta_dir <- rnorm(dimx)
  #beta_dir <- beta_dir / sqrt(sum(beta_dir^2))
  #X_centered <- scale(X, center = TRUE, scale = FALSE)
  #Xbeta_temp <- X_centered %*% beta_dir
  #target_var <- 1.2^2 - sd_trend^2-sd_y^2  
  #beta_scale <- as.numeric(sqrt(target_var / var(Xbeta_temp)))
  #newbeta <- beta_dir * beta_scale  # now sd(X %*% beta) ≈ 1.41

  # Compute the mean of X %*% beta
  #mean_X <- colMeans(X)
  #mean_scale <- sum(mean_X * newbeta)

  # Adjust beta to also make mean(X %*% beta) ≈ 10
  #newbeta <- newbeta * (11 / mean_scale)  
  #beta <- c(beta[1],newbeta)  
  Xtilde <- X  
  for(i in 1:dimx){
      Xtilde[,i] <- (Xtilde[,i] - mean(Xtilde[,i]))/sd(Xtilde[,i])
  }
  Xtilde <- cbind(rep(1,T),Xtilde)  
  X <- cbind(rep(1,T),X)
   
  # (3): simulate y
  trend=rep(0,T) #Trend component
  delta=rep(0,T) #Slope
  trend[1] <- rnorm(1,mean=0,sd=sd_trend)
  delta[1] <- rnorm(1,mean=0,sd=sd_trend)
  if(ar == FALSE){
    # SET A LOCAL LEVEL MODEL FOR TREND
    for (i in 2:T){
      trend[i]<-trend[i-1]+delta[i-1]+rnorm(n=1,mean=mean_trend,sd=sd_trend)
      if(rho!=0){
        delta[i]<-Dtilde+rho*(delta[i-1]-Dtilde)+rnorm(n=1,mean=mean_trend,sd=sd_trend)
      }else{
        delta[i] <- delta[i-1] + rnorm(n=1,mean=0,sd=0.1)
      }
    }
  }else{
    # SET AR(1) model FOR TREND
    for(i in 2:T){
      trend[i] = phi * trend[i-1] + rnorm(n=1,mean=0,sd=sd_trend)
    }
  }
  y <- Xtilde%*%beta + trend + rnorm(T,0,sd_y)
  y <- y*1.5 + 11
  out <- list(y=y,X=X,beta=beta)
  return(out)
}
##############################
### simulating function ######
##############################

sim <- function(miss=1,M=20,totT=100,sd_mu=0.01,realbeta=c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),correlation=FALSE,seasontrend=FALSE){

  pred.error <- beta.error  <- inc.error <-
    effect <- cumeffect <- p <- spec <- specavg <- 0
  effect.ci <- cumeffect.ci <- c(0,0)
  sen <- p.sen <- apee <- cumeffectall<-avgeffectall<- senavg <- avgeffectall.ci.lower <- avgeffectall.ci.upper<-cumeffectall.ci.lower <- cumeffectall.ci.upper <- rep(0,7)
  realeffect <- matrix(0,nrow=2,ncol=7)
    newdata <-datagen.arima.ar1(T=totT,dimx=M,sig.x=sig.x,
                            mu.x=rnorm(M,0,0.01),                           beta=realbeta,#c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),#c(rep(1/3,4),0,0),#c(0.2,0.3,0.8,-0.03),#c(0.1,-0.5,1,0.5,rep(0,2)),
                            ar = TRUE,sd_trend = sd_mu,sd_y = 0.5,
                           r = rep(0.9,M),corr=correlation,phi=1)
    #newdata <-datagen.arima.ar1(T=100,dimx=20,sig.x=sig.x,
    #                        mu.x=rnorm(20,0,0.01),                           beta=c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),#c(rep(1/3,4),0,0),#c(0.2,0.3,0.8,-0.03),#c(0.1,-0.5,1,0.5,rep(0,2)),
  #                          ar = TRUE,sd_trend = 0.01,sd_y = 0.5,
  #                         r = rep(0.9,20),corr=TRUE,phi=1)
    if(seasontrend == TRUE){
         newdata <- datagen.arima.season(T=totT,dimx = M,sig.x,mu.x=rnorm(M,0,0.1),
                                         beta=realbeta,sd_trend = sd_mu,r = rep(0.99, M),
                                        corr=correlation,phi=1,sd_y=0.5)
    }
    realbeta <- newdata$beta
 
 
    ys <- ty <- c()
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
    
    
  dup <- 100
  # non missing
  
  if(miss == 1){
    for(j in 1:totT){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y[j],6))
        ty <- c(ty,rep(j,dup))
    }
  }else{
    if(miss == 2){
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y[j],6)
        newsample <- c()
        for(s in sample){
          if(rbinom(1,n=1,p=0.4) == 0){
            newsample <- c(newsample,s)
          }
        }
        if(length(newsample) == 0){
          newsample <- sample
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
    }else{
      # MNAR
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y[j],6)
        newsample <- c()
        for(s in sample){
          if(s >=13){
      if(rbinom(1,n=1,p=0.5) == 0){
        newsample <- c(newsample,s)
      }
    }else{
      if(s < 4){
        if(rbinom(1,n=1,p=0.7) == 0){
          newsample <- c(newsample,s)
        }
      }else{
        newsample <- c(newsample,s)
      }
    }
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
    }
  }
  
  simy <- data.frame(y=ys,t=ty)

  simy_agg <- aggregate(y ~t,data=simy, mean)
  simx <- newdata$X[1:totT,-1]
  ind <- xs <- t <- c()
  sdx <- rtruncnorm(M,a=0,b=10,6.5,3)#c(7,6,3,4,5)#rnorm(M,6,1)
  #sdx <- c(7.4,8.2,7.3)  
  dupx <- rep(100,M)
  if(miss == 1){
    # NON MISSING
    for(i in 1:M){
      for(j in 1:totT){
        xs <-c(xs,rtruncnorm(dupx[i],a=0,b=36,simx[j,i],sdx[i]))
        t <- c(t,rep(j,dupx[i]))
        ind <- c(ind,rep(i,dupx[i]))
      }
    }
  }else{
    if(miss == 2){
      # MCAR
      for(i in 1:M){
        for(j in 1:totT){
          sample <- rtruncnorm(dup,a=0,b=36,simx[j,i],sdx[i])
          newsample <- c()
          for(s in sample){
            if(rbinom(1,n=1,p=0.4) == 0){
              newsample <- c(newsample,s)
            }
          }
          if(length(newsample) == 0){
            newsample <- sample
          }
          xs <- c(xs,newsample)
          t <- c(t,rep(j,length(newsample)))
          ind <- c(ind,rep(i,length(newsample)))
        }
      }
    }else{
      for(i in 1:M){
        #p1 <- runif(1,0,1)
        #p2 <- runif(1,0,1)
        for(j in 1:totT){
          sample <- rtruncnorm(dup,a=0,b=36,simx[j,i],sdx[i])
          newsample <- c()
          for(s in sample){
            if(s >=round(12+0.08*(i-1))){
      if(rbinom(1,n=1,p=0.4+0.01*(i-1)) == 0){
        newsample <- c(newsample,s)
      }
    }else{
      if(s < round(3+0.08*(i-1))){
        if(rbinom(1,n=1,p=0.6+0.01*(i-1)) == 0){
          newsample <- c(newsample,s)
        }
      }else{
        newsample <- c(newsample,s)
      }
    }
          }
          xs <- c(xs,newsample)
          t <- c(t,rep(j,length(newsample)))
          ind <- c(ind,rep(i,length(newsample)))
        }
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
  
  

  if(miss == 1){
    ys <- ty <- c()
    # non missing
    for(j in 1:totT){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y1[j],6))
      ty <- c(ty,rep(j,dup))
    }
    simy <- data.frame(y=ys,t=ty)
    simy_agg <- aggregate(y ~t,data=simy, mean)
    simy1 <- simy_agg$y
    
    ys <- ty <- c()
    for(j in 1:totT){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y10[j],6))
      ty <- c(ty,rep(j,dup))
    }
    simy <- data.frame(y=ys,t=ty)
    simy_agg <- aggregate(y ~t,data=simy, mean)
    simy10 <- simy_agg$y
    
    ys <- ty <- c()
    for(j in 1:totT){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y30[j],6))
      ty <- c(ty,rep(j,dup))
    }
    simy <- data.frame(y=ys,t=ty)
    simy_agg <- aggregate(y ~t,data=simy, mean)
    simy30 <- simy_agg$y
    
    ys <- ty <- c()
    for(j in 1:totT){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y50[j],6))
      ty <- c(ty,rep(j,dup))
    }
    simy <- data.frame(y=ys,t=ty)
    simy_agg <- aggregate(y ~t,data=simy, mean)
    simy50 <- simy_agg$y
    
    ys <- ty <- c()
    for(j in 1:totT){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y100[j],6))
      ty <- c(ty,rep(j,dup))
    }
    simy <- data.frame(y=ys,t=ty)
    simy_agg <- aggregate(y ~t,data=simy, mean)
    simy100 <- simy_agg$y
    
    ys <- ty <- c()
    for(j in 1:totT){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y200[j],6))
      ty <- c(ty,rep(j,dup))
    }
    simy <- data.frame(y=ys,t=ty)
    simy_agg <- aggregate(y ~t,data=simy, mean)
    simy200 <- simy_agg$y
    
    ys <- ty <- c()
    for(j in 1:totT){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y300[j],6))
      ty <- c(ty,rep(j,dup))
    }
    simy <- data.frame(y=ys,t=ty)
    simy_agg <- aggregate(y ~t,data=simy, mean)
    simy300 <- simy_agg$y
  }else{
    if(miss == 2){
      #### ----  MCAR ----
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y1[j],6)
        newsample <- c()
        for(s in sample){
          if(rbinom(1,n=1,p=0.4) == 0){
            newsample <- c(newsample,s)
          }
        }
        if(length(newsample) == 0){
          newsample <- sample
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy1 <- simy_agg$y
      
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y10[j],6)
        newsample <- c()
        for(s in sample){
          if(rbinom(1,n=1,p=0.4) == 0){
            newsample <- c(newsample,s)
          }
        }
        if(length(newsample) == 0){
          newsample <- sample
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy10 <- simy_agg$y
      
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y30[j],6)
        newsample <- c()
        for(s in sample){
          if(rbinom(1,n=1,p=0.4) == 0){
            newsample <- c(newsample,s)
          }
        }
        if(length(newsample) == 0){
          newsample <- sample
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy30 <- simy_agg$y
      
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y50[j],6)
        newsample <- c()
        for(s in sample){
          if(rbinom(1,n=1,p=0.4) == 0){
            newsample <- c(newsample,s)
          }
        }
        if(length(newsample) == 0){
          newsample <- sample
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy50 <- simy_agg$y
      
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y100[j],6)
        newsample <- c()
        for(s in sample){
          if(rbinom(1,n=1,p=0.4) == 0){
            newsample <- c(newsample,s)
          }
        }
        if(length(newsample) == 0){
          newsample <- sample
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy100 <- simy_agg$y
      
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y200[j],6)
        newsample <- c()
        for(s in sample){
          if(rbinom(1,n=1,p=0.4) == 0){
            newsample <- c(newsample,s)
          }
        }
        if(length(newsample) == 0){
          newsample <- sample
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy200 <- simy_agg$y
      
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y300[j],6)
        newsample <- c()
        for(s in sample){
          if(rbinom(1,n=1,p=0.4) == 0){
            newsample <- c(newsample,s)
          }
        }
        if(length(newsample) == 0){
          newsample <- sample
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy300 <- simy_agg$y
    }else{
      ## ----- MNAR ---
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y1[j],6)
        newsample <- c()
        for(s in sample){
          if(s >=13){
      if(rbinom(1,n=1,p=0.5) == 0){
        newsample <- c(newsample,s)
      }
    }else{
      if(s < 4){
        if(rbinom(1,n=1,p=0.7) == 0){
          newsample <- c(newsample,s)
        }
      }else{
        newsample <- c(newsample,s)
      }
    }
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy1 <- simy_agg$y
      
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y10[j],6)
        newsample <- c()
        for(s in sample){
          if(s >=13){
      if(rbinom(1,n=1,p=0.5) == 0){
        newsample <- c(newsample,s)
      }
    }else{
      if(s < 4){
        if(rbinom(1,n=1,p=0.7) == 0){
          newsample <- c(newsample,s)
        }
      }else{
        newsample <- c(newsample,s)
      }
    }
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy10 <- simy_agg$y
      
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y30[j],6)
        newsample <- c()
        for(s in sample){
          if(s >=13){
      if(rbinom(1,n=1,p=0.5) == 0){
        newsample <- c(newsample,s)
      }
    }else{
      if(s < 4){
        if(rbinom(1,n=1,p=0.7) == 0){
          newsample <- c(newsample,s)
        }
      }else{
        newsample <- c(newsample,s)
      }
    }
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy30 <- simy_agg$y
      
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y50[j],6)
        newsample <- c()
        for(s in sample){
          if(s >=13){
      if(rbinom(1,n=1,p=0.5) == 0){
        newsample <- c(newsample,s)
      }
    }else{
      if(s < 4){
        if(rbinom(1,n=1,p=0.7) == 0){
          newsample <- c(newsample,s)
        }
      }else{
        newsample <- c(newsample,s)
      }
    }
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy50 <- simy_agg$y
      
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y100[j],6)
        newsample <- c()
        for(s in sample){
          if(s >=13){
      if(rbinom(1,n=1,p=0.5) == 0){
        newsample <- c(newsample,s)
      }
    }else{
      if(s < 4){
        if(rbinom(1,n=1,p=0.7) == 0){
          newsample <- c(newsample,s)
        }
      }else{
        newsample <- c(newsample,s)
      }
    }
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy100 <- simy_agg$y
      
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y200[j],6)
        newsample <- c()
        for(s in sample){
          if(s >=13){
      if(rbinom(1,n=1,p=0.5) == 0){
        newsample <- c(newsample,s)
      }
    }else{
      if(s < 4){
        if(rbinom(1,n=1,p=0.7) == 0){
          newsample <- c(newsample,s)
        }
      }else{
        newsample <- c(newsample,s)
      }
    }
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy200 <- simy_agg$y
      
      ys <- ty <- c()
      for(j in 1:totT){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y300[j],6)
        newsample <- c()
        for(s in sample){
          if(s >=13){
      if(rbinom(1,n=1,p=0.5) == 0){
        newsample <- c(newsample,s)
      }
    }else{
      if(s < 4){
        if(rbinom(1,n=1,p=0.7) == 0){
          newsample <- c(newsample,s)
        }
      }else{
        newsample <- c(newsample,s)
      }
    }
        }
        ys <-c(ys,newsample)
        ty <- c(ty,rep(j,length(newsample)))
      }
      simy <- data.frame(y=ys,t=ty)
      simy_agg <- aggregate(y ~t,data=simy, mean)
      simy300 <- simy_agg$y
      
    }
  }
 

  
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

# ca_20_100_miss1
res <- mclapply(1:B,function(b) {sim(miss = 1,M = 20,totT = 100,sd_mu = 0.01,
                                     realbeta = c(rep(1/5, 6), rep(-1/5, 5), rep(1/10, 10)),correlation = FALSE)},
                mc.cores = ncores
)
ca_20_100_miss1 <- clean(res)
saveRDS(ca_20_100_miss1,"ca_20_100_miss1.RData")
# ca_20_100_miss2
res <- mclapply(1:B,function(b) {sim(miss = 2,M = 20,totT = 100,sd_mu = 0.01,
                                     realbeta = c(rep(1/5, 6), rep(-1/5, 5), rep(1/10, 10)),correlation = FALSE)},
                mc.cores = ncores
)
ca_20_100_miss2 <- clean(res)
saveRDS(ca_20_100_miss2,"ca_20_100_miss2.RData")
# ca_20_100_miss3
set.seed(123)
res <- mclapply(1:B,function(b) {sim(miss = 3,M = 20,totT = 100,sd_mu = 0.01,
                                     realbeta = c(rep(1/5, 6), rep(-1/5, 5), rep(1/10, 10)),correlation = FALSE)},
                mc.cores = ncores
)
ca_20_100_miss3 <- clean(res)
saveRDS(ca_20_100_miss3,"ca_20_100_miss3.RData")
# ca_20_100_corr
res <- mclapply(1:B,function(b) {sim(miss=1,M=20,totT=100,sd_mu=0.01,realbeta=c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),correlation=TRUE)},
                mc.cores = ncores)
ca_20_100_corr <- clean(res)
saveRDS(ca_20_100_corr,"ca_20_100_corr.RData")

# ca_20_100_season
set.seed(123)
res <- mclapply(1:B,function(b) {sim(miss=1,M=20,totT=100,sd_mu=0.01,realbeta=c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),correlation=FALSE,seasontrend=TRUE)},
                mc.cores = ncores)
ca_20_100_season <- clean(res)
saveRDS(ca_20_100_season,"ca_20_100_season.RData")
# ca_20_100_mu1
set.seed(123)
res <- mclapply(1:B,function(b) {sim(miss=1,M=20,totT=100,sd_mu=1,realbeta=c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),correlation=FALSE)},
                mc.cores = ncores)
ca_20_100_mu1 <- clean(res)
saveRDS(ca_20_100_mu1,"ca_20_100_mu1.RData")

# ca_20_100_mu05
set.seed(123)
res <- mclapply(1:B,function(b) {sim(miss=1,M=20,totT=100,sd_mu=0.5,realbeta=c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),correlation=FALSE)},
                mc.cores = ncores)
ca_20_100_mu05 <- clean(res)
saveRDS(ca_20_100_mu05,"ca_20_100_mu05.RData")

# ca_20_100_mu15
set.seed(123)
res <- mclapply(1:B,function(b) {sim(miss=1,M=20,totT=100,sd_mu=1.5,realbeta=c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),correlation=FALSE)},
                mc.cores = ncores)
ca_20_100_mu15 <- clean(res)
saveRDS(ca_20_100_mu15,"ca_20_100_mu15.RData")

# ca_20_100_mu2
set.seed(123)
res <- mclapply(1:B,function(b) {sim(miss=1,M=20,totT=100,sd_mu=2,realbeta=c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),correlation=FALSE)},
                mc.cores = ncores)
ca_20_100_mu2 <- clean(res)
saveRDS(ca_20_100_mu2,"ca_20_100_mu2.RData")

# ca_20_50
res <- mclapply(1:B,function(b) {sim(miss = 1,M = 20,totT = 50,sd_mu = 0.01,
                                     realbeta = c(rep(1/5, 6), rep(-1/5, 5), rep(1/10, 10)),correlation = FALSE)},
                mc.cores = ncores
)
ca_20_50 <- clean(res)
saveRDS(ca_20_50,"ca_20_50.RData")
# ca_20_200
res <- mclapply(1:50,function(b) {sim(miss = 1,M = 20,totT = 200,sd_mu = 0.01,
                                     realbeta = c(rep(1/5, 6), rep(-1/5, 5), rep(1/10, 10)),correlation = FALSE)},
                mc.cores = ncores
)

res2 <- mclapply(1:50,function(b) {sim(miss = 1,M = 20,totT = 200,sd_mu = 0.01,
                                      realbeta = c(rep(1/5, 6), rep(-1/5, 5), rep(1/10, 10)),correlation = FALSE)},
                mc.cores = ncores
)

res3 <- mclapply(1:50,function(b) {sim(miss = 1,M = 20,totT = 200,sd_mu = 0.01,
                                       realbeta = c(rep(1/5, 6), rep(-1/5, 5), rep(1/10, 10)),correlation = FALSE)},
                 mc.cores = ncores
)

res4 <- mclapply(1:50,function(b) {sim(miss = 1,M = 20,totT = 200,sd_mu = 0.01,
                                       realbeta = c(rep(1/5, 6), rep(-1/5, 5), rep(1/10, 10)),correlation = FALSE)},
                 mc.cores = ncores
)
res_all <- c(res,res2,res3,res4)
ca_20_200 <- clean(res_all)
saveRDS(ca_20_200,"ca_20_200.RData")





#---------------------------------M=5 ------------------------------
# ca_5_100_miss1

res <- mclapply(1:B,function(b) {sim(miss = 1,M = 5,totT = 100,sd_mu = 0.01,
                                     realbeta = c(rep(1/3,4),0,0),correlation = FALSE)},
                mc.cores = ncores, mc.set.seed = TRUE
)
ca_5_100_miss1 <- clean(res)
saveRDS(ca_5_100_miss1,"ca_5_100_miss1.RData")
# ca_5_100_miss2
set.seed(123)
res <- mclapply(1:B,function(b) {sim(miss = 2,M = 5,totT = 100,sd_mu = 0.01,
                                     realbeta = c(rep(1/3,4),0,0),correlation = FALSE)},
                mc.cores = ncores, mc.set.seed = TRUE
)
ca_5_100_miss2 <- clean(res)
saveRDS(ca_5_100_miss2,"ca_5_100_miss2.RData")

# ca_5_100_miss3

res <- mclapply(1:B,function(b) {sim(miss = 3,M = 5,totT = 100,sd_mu = 0.01,
                                     realbeta = c(rep(1/3,4),0,0),correlation = FALSE)},
                mc.cores = ncores
)
ca_5_100_miss3 <- clean(res)
saveRDS(ca_5_100_miss3,"ca_5_100_miss3.RData")


# ca_5_100_corr
res <- mclapply(1:B,function(b) {sim(miss = 1,M = 5,totT = 100,sd_mu = 0.01,
                                     realbeta = c(rep(1/3,4),0,0),correlation = TRUE)},
                mc.cores = ncores
)
ca_5_100_corr <- clean(res)
saveRDS(ca_5_100_corr,"ca_5_100_corr.RData")

# ca_5_100_season
set.seed(123)
res <- mclapply(1:B,function(b) {sim(miss=1,M=5,totT=100,sd_mu=0.01,realbeta=c(rep(1/3,4),0,0),correlation=FALSE,seasontrend=TRUE)},
                mc.cores = ncores)
ca_5_100_season <- clean(res)
saveRDS(ca_5_100_season,"ca_5_100_season.RData")

# ca_5_100_mu05
res <- mclapply(1:B,function(b) {sim(miss = 1,M = 5,totT = 100,sd_mu = 0.5,
                                     realbeta = c(rep(1/3,4),0,0),correlation = FALSE)},
                mc.cores = ncores
)
ca_5_100_mu05 <- clean(res)
saveRDS(ca_5_100_mu05,"ca_5_100_mu05.RData")

# ca_5_100_mu1
res <- mclapply(1:B,function(b) {sim(miss = 1,M = 5,totT = 100,sd_mu = 1,
                                     realbeta = c(rep(1/3,4),0,0),correlation = FALSE)},
                mc.cores = ncores
)
ca_5_100_mu1 <- clean(res)
saveRDS(ca_5_100_mu1,"ca_5_100_mu1.RData")

# ca_5_100_mu2
res <- mclapply(1:B,function(b) {sim(miss = 1,M = 5,totT = 100,sd_mu = 2,
                                     realbeta = c(rep(1/3,4),0,0),correlation = FALSE)},
                mc.cores = ncores
)
ca_5_100_mu2 <- clean(res)
saveRDS(ca_5_100_mu2,"ca_5_100_mu2.RData")

# ca_5_100_mu15
res <- mclapply(1:B,function(b) {sim(miss = 1,M = 5,totT = 100,sd_mu = 1.5,
                                     realbeta = c(rep(1/3,4),0,0),correlation = FALSE)},
                mc.cores = ncores
)
ca_5_100_mu15 <- clean(res)
saveRDS(ca_5_100_mu15,"ca_5_100_mu15.RData")

# ca_5_50

res <- mclapply(1:B,function(b) {sim(miss = 1,M = 5,totT = 50,sd_mu = 0.01,
                                     realbeta = c(rep(1/3,4),0,0),correlation = FALSE)},
                mc.cores = ncores
)
ca_5_50 <- clean(res)
saveRDS(ca_5_50,"ca_5_50.RData")

# ca_5_200
res <- lapply(1:200,function(b) {sim(miss = 1,M = 5,totT = 200,sd_mu = 0.01,
                                     realbeta = c(rep(1/3,4),0,0),correlation = FALSE)})
res <- mclapply(1:200,function(b) {sim(miss = 1,M = 5,totT = 200,sd_mu = 0.01,
                                     realbeta = c(rep(1/3,4),0,0),correlation = FALSE)},
                mc.cores = ncores
)
ca_5_200 <- clean(res)
saveRDS(ca_5_200,"ca_5_200.RData")






