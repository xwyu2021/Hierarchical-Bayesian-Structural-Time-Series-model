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


source("functions/impact_additional.R")
source("functions/impact_stan.R")
source("functions/impact_hbsts.R")
source("functions/metrics.R")

source("https://raw.githubusercontent.com/google/CausalImpact/master/R/impact_model.R")
source("https://raw.githubusercontent.com/google/CausalImpact/master/R/impact_misc.R")
source("https://raw.githubusercontent.com/google/CausalImpact/master/R/impact_analysis.R")
source("https://raw.githubusercontent.com/google/CausalImpact/master/R/impact_inference.R")

####################################################
### absolute percentatge error estimation (APEE) ###
####################################################
apee_est <- function(est_ef,real_ef){
  return(abs(est_ef - real_ef)/real_ef)
}

#################################################
### data generating function with seasonality ###
#################################################
datagen.arima.season <- function(T=100,dimx = 5,sig.x,mu.x=rep(1,5),beta,sd_trend = 0.1,
                                 phi=0.5,sd_y=0.1){
  trend=rep(0,T) #Trend component
  X <- matrix(0,nrow=T,ncol=dimx)
  # (1): simulate control series,inc. intercept
  
  xs <- runif(dimx,0.1,4)#c(7,2,0.1,5,8) 
  #xm <- rnorm(dimx,11,0.05)  
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
                              phi=0.5,sd_y=0.1,corr =FALSE){
  trend=rep(0,T) #Trend component
  X <- matrix(0,nrow=T,ncol=dimx)
  if(corr == TRUE){
    R <- matrix(0.85, nrow = dimx, ncol = dimx)
    diag(R) <- 1
    xs <- runif(dimx,0.1,4)  
    Sgm <- diag(xs) %*% R %*% diag(xs)
    X <- mvrnorm(n = T, mu = rnorm(dimx,11,0.05), Sigma = Sgm)  
    
  }else{
    xs <- runif(dimx,0.1,4)
    for(i in 1:dimx){
      X[,i] <- rnorm(T,11,xs[i])
    } 
  }
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

####################################################################
### simulating function (part 1: expected model size) ##############
####################################################################
sim_modelsize <- function(miss){
  newdata <-datagen.arima.ar1(T=100,dimx=20,sig.x=sig.x,
                              mu.x=rnorm(20,0,0.01),                           
                              beta=c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),#c(rep(1/3,4),0,0)
                              ar = TRUE,sd_trend = 0.01,sd_y = 0.5,
                              corr=FALSE,phi=1)
  realbeta <- newdata$beta
  

  
  ys <- ty <- c()
  simemp.mean.y <-  newdata$y
  
  
  dup <- 100
  # non missing
  
  if(miss == 1){
    for(j in 1:100){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y[j],6))
      ty <- c(ty,rep(j,dup))
    }
  }else{
    if(miss == 2){
      for(j in 1:100){
        sample <- rtruncnorm(dup,a=0,b=36,simemp.mean.y[j],6)
        newsample <- c()
        for(s in sample){
          if(rbinom(1,n=1,p=0.6) == 0){
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
      for(j in 1:100){
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
  simx <- newdata$X[1:100,-1]
  ind <- xs <- t <- c()
  sdx <- rtruncnorm(M,a=0,b=10,6.5,3)#c(7,6,3,4,5)
  #sdx <- c(7.4,8.2,7.3)  
  dupx <- rep(100,M)
  if(miss == 1){
    # NON MISSING
    for(i in 1:M){
      for(j in 1:100){
        xs <-c(xs,rtruncnorm(dupx[i],a=0,b=36,simx[j,i],sdx[i]))
        t <- c(t,rep(j,dupx[i]))
        ind <- c(ind,rep(i,dupx[i]))
      }
    }
  }else{
    if(miss == 2){
      # MCAR
      for(i in 1:M){
        for(j in 1:100){
          sample <- rtruncnorm(dupx[i],a=0,b=36,simx[j,i],sdx[i])
          newsample <- c()
          for(s in sample){
            if(rbinom(1,n=1,p=0.6) == 0){
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
        for(j in 1:100){
          sample <- rtruncnorm(dupx[i],a=0,b=36,simx[j,i],sdx[i])
          newsample <- c()
          for(s in sample){
            if(s >=round(12+0.4*(i-1))){
              if(rbinom(1,n=1,p=0.4+0.05*(i-1)) == 0){
                newsample <- c(newsample,s)
              }
            }else{
              if(s < round(3+0.4*(i-1))){
                if(rbinom(1,n=1,p=0.6+0.05*(i-1)) == 0){
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
  simx_agg_df <- dcast(as.data.table(simx_agg),t~ind)
  simx_agg_df <- simx_agg_df[-1]
  
  #stand.simxobs <- simxobs
  #for(id in unique(simxobs$ind)){
  #  pick <- which(simxobs$ind == id)
  #  stand.simxobs$x[pick] <- (simxobs$x[pick]- mean(simxobs$x[pick]))/sd(simxobs$x[pick])
  #}
  simpre.time <- 60
  simpost.time <- 40
  simpre.ys <- as.vector(simy$y[simy$t <=60])
  simpost.ys <- as.vector(simy$y[simy$t >60])
  #muy <- mean(simy$y)
  #sdy <- sd(simy$y)
  #standy <- (simy$y - muy)/sdy
  #simpre.ys <- as.vector(standy[simy$t <=150])
  #simpost.ys <- as.vector(standy[simy$t >150])
  simpre.obs <- sum(simy$t <=60)
  simpost.obs <- sum(simy$t >60)
  
  simobs.count <- as.vector(table(simy$t))
  simindex_pre <- cumsum(c(1,simobs.count[1:simpre.time]))
  
  simindex_post <- cumsum(c(1,simobs.count[(simpre.time+1):(simpre.time+simpost.time)]))
  
  checksim <- as.data.frame(simxobs[,2:3])
  simx.count <- table(checksim)
  simx.cum.count <- t(matrix(cumsum(t(simx.count)),nrow=M,ncol=100))
  simxobs <- as.vector(simxobs[,1])
  #simxobs <- as.vector(stand.simxobs[,1])
  simK <- length(simxobs)
  p.mat <- as.matrix(simx.count)
  p.mat[p.mat!=0] <- 1
  hsim <- run_hbsts(pre.time=simpre.time,post.time=simpost.time,
                    pre.obs=simpre.obs,
                    pre.ys=simpre.ys,dimx=M,
                    post.obs =simpost.obs,index_pre=simindex_pre,
                    K=simK,x.obs=simxobs, x.cum.count=simx.cum.count,index_post=simindex_post,
                    niter=8000,nchains=1,time_index=1:100,presence_matrix = p.mat,
                    alpha=0.05,cycles=c(1),emp.mean.y=simy_agg$y,#emp.mean.y=simemp.mean.y,
                    trend = 'LL',burn=3000,modelsize=meff,adapt=0.96,hdi=FALSE)
  
  
  newD <- cbind(simy_agg$y,as.matrix(simx_agg_df))
  hsmod <- run_stanbsts(c(1,60),c(61,100),data=newD,
                        modelsize=meff,trend='LL',
                        nchain=1,niter=8000,alpha=0.05,cycles=1,adapt=0.9,burn=3000,random = FALSE,realy=simemp.mean.y[61:100])
  
  newD <- as.data.frame(cbind(1:100,newD))
  D <- read.zoo(newD,index.column = 1)
  ssmod <-  run_flex_bsts(data=D,pre.period=c(1,60),
                          post.period=c(61,100),
                          alpha=0.05,
                          trend = 'LL',modelsize=meff,cycles = c(1),random=FALSE,realy=simemp.mean.y[61:100],niter = 8000)
  
  pred.error[1] <- RMSE(hsim$series$point.pred[61:100],simy_agg$y[61:100])
  pred.error[2] <- RMSE(hsmod$series$point.pred[61:100],simy_agg$y[61:100])
  pred.error[3] <- RMSE(ssmod$series$point.pred[61:100],simy_agg$y[61:100])
  
  #effect[1] <- RMSE(hsim$series$point.pred[61:100],simemp.mean.y[61:100])
  #effect[2] <- RMSE(hsmod$series$point.pred[61:100],simemp.mean.y[61:100])
  #effect[3] <- RMSE(ssmod$series$point.pred[61:100],simemp.mean.y[61:100])
  effect[1] <- hsim$summary$AbsEffect[1]
  effect[2] <- hsmod$summary$AbsEffect[1]
  effect[3] <- ssmod$summary$AbsEffect[1]
  
  fit <- rstan::extract(hsim$model$bsts.model)
  fit2 <- rstan::extract(hsmod$model$bsts.model)
  fit3 <- ssmod$model$bsts.model
  
  beta.error[1] <- RMSE(colMeans(cbind(fit$beta)),realbeta[2:(M+1)])
  beta.error[2] <- RMSE(colMeans(cbind(fit2$beta)),realbeta[2:(M+1)])
  beta.error[3] <- RMSE(colMeans(cbind(fit3$coefficients[,2:(M+1)])),realbeta[2:(M+1)])
  
  inc.error[1] <- RMSE(colMeans(cbind(fit$rho!=0)),realbeta[2:(M+1)]!=0)
  inc.error[2] <- RMSE(colMeans(cbind(fit2$rho!=0)),realbeta[2:(M+1)]!=0)
  inc.error[3] <- RMSE(colMeans(cbind(fit3$coefficients[,2:(M+1)]!=0)),realbeta[2:(M+1)]!=0)
  
  MAPE[1] <- mean(abs((hsim$series$point.pred[61:100]-simy_agg$y[61:100])/simy_agg$y[61:100])) 
  MAPE[2] <- mean(abs((hsmod$series$point.pred[61:100]-simy_agg$y[61:100])/simy_agg$y[61:100])) 
  MAPE[3] <- mean(abs((ssmod$series$point.pred[61:100]-simy_agg$y[61:100])/simy_agg$y[61:100])) 
  
  SMAPE[1] <-  mean(abs(hsim$series$point.pred[61:100]-simy_agg$y[61:100])/(abs(simy_agg$y[61:100]) + abs(hsim$series$point.pred[61:100])))
  SMAPE[2] <-  mean(abs(hsmod$series$point.pred[61:100]-simy_agg$y[61:100])/(abs(simy_agg$y[61:100]) + abs(hsmod$series$point.pred[61:100])))
  SMAPE[3] <-  mean(abs(ssmod$series$point.pred[61:100]-simy_agg$y[61:100])/(abs(simy_agg$y[61:100]) + abs(ssmod$series$point.pred[61:100])))
  
  return(list(pred.error = pred.error,beta.error=beta.error,inc.error=inc.error,effect=effect,
              abs.error=abs.error, MAPE=MAPE, SMAPE=SMAPE))
}




M=20
realbeta <- c(rep(1/5,6),rep(-1/5,5),rep(1/10,10))#c(rep(1/3,4),0,0)

## --- RUN THIS FUNCTION ---
for(i in c(1:20)){
  
  pred.error <- beta.error  <- inc.error <- abs.error <- 
    effect <- MAPE <- SMAPE <- rep(0,3)
  
  
  
  #source('/root/autodl-tmp/lse project/simulation100_short_modsize.R')
  meff <- i
  miss <- 1
  results <- mclapply(rep(1,200), sim_modelsize,mc.cores = numCores)
  
  emp <- c()
  for(l in 1:length(results)){
    if(is.null(results[[l]])){
      emp <- c(emp,l)
    }
  }
  if(length(emp) > 0){
    results <- results[-emp]
  }
  U <- length(results)
  
  pred.error <- t(sapply(results,function(l) l$pred.error))
  beta.error <- t(sapply(results,function(l) l$beta.error))
  inc.error <- t(sapply(results,function(l) l$inc.error))
  abs.error <- t(sapply(results,function(l) l$abs.error))
  MAPE <- t(sapply(results,function(l) l$MAPE))
  SMAPE <- t(sapply(results,function(l) l$SMAPE))
  effect <- t(sapply(results,function(l) l$effect))
  simNew_result <- list(pred.error = pred.error,beta.error=beta.error,
                        inc.error=inc.error,effect=effect,abs.error=abs.error,
                        MAPE=MAPE,SMAPE=SMAPE)             
  
  saveRDS(simNew_result,paste0('smallsim100_20N',i,'.RData'))
}





####################################################################
### simulating function (part 2: sensitivity and specificity) ######
####################################################################

sim <- function(miss){
  newdata <-datagen.arima.ar1(T=100,dimx=20,sig.x=sig.x,
                              mu.x=rnorm(20,0,0.01),           
                              beta=c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),#c(rep(1/3,4),0,0),
                              ar = TRUE,sd_trend = 2,sd_y = 0.5,
                              corr=FALSE,phi=1)
  #if(seasontrend == TRUE){
  #     newdata <- datagen.arima.season(T=100,dimx = 20,sig.x,mu.x=rnorm(20,0,0.1),
  #                                     beta=c(rep(1/5,6),rep(-1/5,5),rep(1/10,10)),
  #                                     sd_trend = 0.01,
  #                                     corr=FALSE,phi=1,sd_y=0.5)
  #}
  realbeta <- newdata$beta
  ys <- ty <- c()
  simemp.mean.y <-  newdata$y
  simemp.mean.y1 <- c(simemp.mean.y[1:60],simemp.mean.y[61:100] + abs(simemp.mean.y[61:100])*0.01)
  simemp.mean.y10 <- c(simemp.mean.y[1:60],simemp.mean.y[61:100] + abs(simemp.mean.y[61:100])*0.03)
  simemp.mean.y30 <- c(simemp.mean.y[1:60],simemp.mean.y[61:100]+ abs(simemp.mean.y[61:100])*0.05)
  simemp.mean.y50 <- c(simemp.mean.y[1:60],simemp.mean.y[61:100]+ abs(simemp.mean.y[61:100])*0.07)
  simemp.mean.y100 <- c(simemp.mean.y[1:60],simemp.mean.y[61:100]+ abs(simemp.mean.y[61:100]*0.1))
  simemp.mean.y200 <- c(simemp.mean.y[1:60],simemp.mean.y[61:100]+ abs(simemp.mean.y[61:100])*0.3)
  simemp.mean.y300 <- c(simemp.mean.y[1:60],simemp.mean.y[61:100]+ abs(simemp.mean.y[61:100])*0.5)
  
  ci1 <- mean(simemp.mean.y1[61:100] - simemp.mean.y[61:100])
  ci10 <- mean(simemp.mean.y10[61:100] - simemp.mean.y[61:100])
  ci30 <- mean(simemp.mean.y30[61:100] - simemp.mean.y[61:100])
  ci50 <- mean(simemp.mean.y50[61:100] - simemp.mean.y[61:100])
  ci100 <- mean(simemp.mean.y100[61:100] - simemp.mean.y[61:100])
  ci200 <- mean(simemp.mean.y200[61:100] - simemp.mean.y[61:100])
  ci300 <- mean(simemp.mean.y300[61:100] - simemp.mean.y[61:100])
    
  
  realeffect[1,] <- c(ci1,ci10,ci30,ci50,ci100,ci200,ci300)
  realeffect[2,] <- c(sum(simemp.mean.y1[61:100] - simemp.mean.y[61:100]),
                      sum(simemp.mean.y10[61:100] - simemp.mean.y[61:100]),
                      sum(simemp.mean.y30[61:100] - simemp.mean.y[61:100]),
                      sum(simemp.mean.y50[61:100] - simemp.mean.y[61:100]),
                      sum(simemp.mean.y100[61:100] - simemp.mean.y[61:100]),
                      sum(simemp.mean.y200[61:100] - simemp.mean.y[61:100]),
                      sum(simemp.mean.y300[61:100] - simemp.mean.y[61:100]))
    
  dup <- 100
  # non missing
  
  if(miss == 1){
    for(j in 1:100){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y[j],6))
      ty <- c(ty,rep(j,dup))
    }
  }else{
    if(miss == 2){
      for(j in 1:100){
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
      for(j in 1:100){
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
  simx <- newdata$X[1:100,-1]
  ind <- xs <- t <- c()
  sdx <- rtruncnorm(20,a=0,b=10,6.5,3)#c(7,6,3,4,5)#rnorm(M,6,1)
  #sdx <- c(7.4,8.2,7.3)  
  dupx <- rep(100,20)
  if(miss == 1){
    # NON MISSING
    for(i in 1:M){
      for(j in 1:100){
        xs <-c(xs,rtruncnorm(dupx[i],a=0,b=36,simx[j,i],sdx[i]))
        t <- c(t,rep(j,dupx[i]))
        ind <- c(ind,rep(i,dupx[i]))
      }
    }
  }else{
    if(miss == 2){
      # MCAR
      for(i in 1:M){
        for(j in 1:100){
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
        for(j in 1:100){
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
  simx_agg_df <- dcast(as.data.table(simx_agg),t~ind)
  simx_agg_df <- simx_agg_df[-1]
  
  #stand.simxobs <- simxobs
  #for(id in unique(simxobs$ind)){
  #  pick <- which(simxobs$ind == id)
  #  stand.simxobs$x[pick] <- (simxobs$x[pick]- mean(simxobs$x[pick]))/sd(simxobs$x[pick])
  #}
  simpre.time <- 60
  simpost.time <- 40
  simpre.ys <- as.vector(simy$y[simy$t <=60])
  simpost.ys <- as.vector(simy$y[simy$t >60])
  #muy <- mean(simy$y)
  #sdy <- sd(simy$y)
  #standy <- (simy$y - muy)/sdy
  #simpre.ys <- as.vector(standy[simy$t <=150])
  #simpost.ys <- as.vector(standy[simy$t >150])
  simpre.obs <- sum(simy$t <=60)
  simpost.obs <- sum(simy$t >60)
  
  simobs.count <- as.vector(table(simy$t))
  simindex_pre <- cumsum(c(1,simobs.count[1:simpre.time]))
  
  simindex_post <- cumsum(c(1,simobs.count[(simpre.time+1):(simpre.time+simpost.time)]))
  
  checksim <- as.data.frame(simxobs[,2:3])
  simx.count <- table(checksim)
  simx.cum.count <- t(matrix(cumsum(t(simx.count)),nrow=M,ncol=100))
  simxobs <- as.vector(simxobs[,1])
  #simxobs <- as.vector(stand.simxobs[,1])
  simK <- length(simxobs)
  p.mat <- as.matrix(simx.count)
  p.mat[p.mat!=0] <- 1
  hsim <- run_hbsts(pre.time=simpre.time,post.time=simpost.time,
                    pre.obs=simpre.obs,
                    pre.ys=simpre.ys,dimx=M,
                    post.obs =simpost.obs,index_pre=simindex_pre,
                    K=simK,x.obs=simxobs, x.cum.count=simx.cum.count,index_post=simindex_post,
                    niter=8000,nchains=1,time_index=1:100,presence_matrix = p.mat,
                    alpha=0.05,cycles=c(1),emp.mean.y=simy_agg$y,#emp.mean.y=simemp.mean.y,
                    trend = 'LL',burn=3000,modelsize=meff,adapt=0.96,hdi=FALSE)
  
  
  newD <- cbind(simy_agg$y,as.matrix(simx_agg_df))
  hsmod <- run_stanbsts(c(1,60),c(61,100),data=newD,
                        modelsize=meff,trend='LL',
                        nchain=1,niter=8000,alpha=0.05,cycles=1,adapt=0.9,burn=3000,random = FALSE,realy=simemp.mean.y[61:100])
  
  newD <- as.data.frame(cbind(1:100,newD))
  D <- read.zoo(newD,index.column = 1)
  ssmod <-  run_flex_bsts(data=D,pre.period=c(1,60),
                          post.period=c(61,100),
                          alpha=0.05,
                          trend = 'LL',modelsize=meff,cycles = c(1),random=FALSE,realy=simemp.mean.y[61:100],niter = 8000)
  
  pred.error[1] <- RMSE(hsim$series$point.pred[61:100],simemp.mean.y[61:100])
  pred.error[2] <- RMSE(hsmod$series$point.pred[61:100],simemp.mean.y[61:100])
  pred.error[3] <- RMSE(ssmod$series$point.pred[61:100],simemp.mean.y[61:100])
  
  fit <- rstan::extract(hsim$model$bsts.model)
  fit2 <- rstan::extract(hsmod$model$bsts.model)
  fit3 <- ssmod$model$bsts.model
  
  beta.error[1] <- RMSE(colMeans(cbind(fit$rho)),realbeta[2:(M+1)])
  beta.error[2] <- RMSE(colMeans(cbind(fit2$rho)),realbeta[2:(M+1)])
  beta.error[3] <- RMSE(colMeans(cbind(fit3$coefficients[,2:(M+1)])),realbeta[2:(M+1)])
  
  inc.error[1] <- RMSE(colMeans(cbind(fit$rho!=0)),realbeta[2:(M+1)]!=0)
  inc.error[2] <- RMSE(colMeans(cbind(fit2$rho!=0)),realbeta[2:(M+1)]!=0)
  inc.error[3] <- RMSE(colMeans(cbind(fit3$coefficients[,2:(M+1)]!=0)),realbeta[2:(M+1)]!=0)
  
  effect[1] <- hsim$summary$AbsEffect[1]
  effect[2] <- hsmod$summary$AbsEffect[1]
  effect[3] <- ssmod$summary$AbsEffect[1]
  
  effect.ci[1,] <- c(hsim$summary$AbsEffect.lower[1],hsim$summary$AbsEffect.upper[1])
  effect.ci[2,] <- c(hsmod$summary$AbsEffect.lower[1],hsmod$summary$AbsEffect.upper[1])
  effect.ci[3,] <- c(ssmod$summary$AbsEffect.lower[1],ssmod$summary$AbsEffect.upper[1])
  
  cumeffect[1] <- hsim$summary$AbsEffect[2]
  cumeffect[2] <- hsmod$summary$AbsEffect[2]
  cumeffect[3] <- ssmod$summary$AbsEffect[2]
  
  cumeffect.ci[1,] <- c(hsim$summary$AbsEffect.lower[2],hsim$summary$AbsEffect.upper[2])
  cumeffect.ci[2,] <- c(hsmod$summary$AbsEffect.lower[2],hsmod$summary$AbsEffect.upper[2])
  cumeffect.ci[3,] <- c(ssmod$summary$AbsEffect.lower[2],ssmod$summary$AbsEffect.upper[2])
  
  p[1] <- hsim$summary$p[1]
  p[2] <- hsmod$summary$p[1]
  p[3] <- ssmod$summary$p[1]
  
  spec[1] <- as.numeric((0>= c(hsim$summary$AbsEffect.lower[2]) & (0 <= hsim$summary$AbsEffect.upper[2])))
  spec[2] <- as.numeric((0>= c(hsmod$summary$AbsEffect.lower[2]) & (0 <= hsmod$summary$AbsEffect.upper[2])))
  spec[3] <- as.numeric((0>= c(ssmod$summary$AbsEffect.lower[2]) & (0 <= ssmod$summary$AbsEffect.upper[2])))
  
  specavg[1] <- as.numeric((0>= c(hsim$summary$AbsEffect.lower[1]) & (0 <= hsim$summary$AbsEffect.upper[1])))
  specavg[2] <- as.numeric((0>= c(hsmod$summary$AbsEffect.lower[1]) & (0 <= hsmod$summary$AbsEffect.upper[1])))
  specavg[3] <- as.numeric((0>= c(ssmod$summary$AbsEffect.lower[1]) & (0 <= ssmod$summary$AbsEffect.upper[1])))
  
  
  if(miss == 1){
    ys <- ty <- c()
    # non missing
    for(j in 1:100){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y1[j],6))
      ty <- c(ty,rep(j,dup))
    }
    simy <- data.frame(y=ys,t=ty)
    simy_agg <- aggregate(y ~t,data=simy, mean)
    simy1 <- simy_agg$y
    
    ys <- ty <- c()
    for(j in 1:100){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y10[j],6))
      ty <- c(ty,rep(j,dup))
    }
    simy <- data.frame(y=ys,t=ty)
    simy_agg <- aggregate(y ~t,data=simy, mean)
    simy10 <- simy_agg$y
    
    ys <- ty <- c()
    for(j in 1:100){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y30[j],6))
      ty <- c(ty,rep(j,dup))
    }
    simy <- data.frame(y=ys,t=ty)
    simy_agg <- aggregate(y ~t,data=simy, mean)
    simy30 <- simy_agg$y
    
    ys <- ty <- c()
    for(j in 1:100){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y50[j],6))
      ty <- c(ty,rep(j,dup))
    }
    simy <- data.frame(y=ys,t=ty)
    simy_agg <- aggregate(y ~t,data=simy, mean)
    simy50 <- simy_agg$y
    
    ys <- ty <- c()
    for(j in 1:100){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y100[j],6))
      ty <- c(ty,rep(j,dup))
    }
    simy <- data.frame(y=ys,t=ty)
    simy_agg <- aggregate(y ~t,data=simy, mean)
    simy100 <- simy_agg$y
    
    ys <- ty <- c()
    for(j in 1:100){
      ys <-c(ys,rtruncnorm(dup,a=0,b=36,simemp.mean.y200[j],6))
      ty <- c(ty,rep(j,dup))
    }
    simy <- data.frame(y=ys,t=ty)
    simy_agg <- aggregate(y ~t,data=simy, mean)
    simy200 <- simy_agg$y
    
    ys <- ty <- c()
    for(j in 1:100){
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
      for(j in 1:100){
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
      for(j in 1:100){
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
      for(j in 1:100){
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
      for(j in 1:100){
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
      for(j in 1:100){
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
      for(j in 1:100){
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
      for(j in 1:100){
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
      for(j in 1:100){
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
      for(j in 1:100){
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
      for(j in 1:100){
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
      for(j in 1:100){
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
      for(j in 1:100){
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
      for(j in 1:100){
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
      for(j in 1:100){
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
  
  #simy1[1:150] <- simy10[1:150] <- simy30[1:150] <- simy50[1:150] <- simy100[1:150] <- simy200[1:150] <-simy300[1:150] <-simy_agg0[1:150]
  ef1 <- effect.hbsts(hsim$model$bsts.model,pre.time=simpre.time,post.time=simpost.time,
                      niter=8000,burn=3000,nchains=1,alpha=0.05,emp.mean.y=simy1,hdi=FALSE)
  ef10 <- effect.hbsts(hsim$model$bsts.model,pre.time=simpre.time,post.time=simpost.time,
                       niter=8000,burn=3000,nchains=1,alpha=0.05,emp.mean.y=simy10,hdi=FALSE)
  ef30 <- effect.hbsts(hsim$model$bsts.model,pre.time=simpre.time,post.time=simpost.time,
                       niter=8000,burn=3000,nchains=1,alpha=0.05,emp.mean.y=simy30,hdi=FALSE)
  ef50 <- effect.hbsts(hsim$model$bsts.model,pre.time=simpre.time,post.time=simpost.time,
                       niter=8000,burn=3000,nchains=1,alpha=0.05,emp.mean.y=simy50,hdi=FALSE)
  ef100 <- effect.hbsts(hsim$model$bsts.model,pre.time=simpre.time,post.time=simpost.time,
                        niter=8000,burn=3000,nchains=1,alpha=0.05,emp.mean.y=simy100,hdi=FALSE)
  ef200 <- effect.hbsts(hsim$model$bsts.model,pre.time=simpre.time,post.time=simpost.time,
                        niter=8000,burn=3000,nchains=1,alpha=0.05,emp.mean.y=simy200,hdi=FALSE)
  ef300 <- effect.hbsts(hsim$model$bsts.model,pre.time=simpre.time,post.time=simpost.time,
                        niter=8000,burn=3000,nchains=1,alpha=0.05,emp.mean.y=simy300,hdi=FALSE)
  
  sen[1,1] <- as.numeric((0 <= ef1$AbsEffect.lower[2])|(0>=ef1$AbsEffect.upper[2]))
  sen[1,2] <- as.numeric((0 <= ef10$AbsEffect.lower[2])|(0>=ef10$AbsEffect.upper[2]))
  sen[1,3] <- as.numeric((0 <= ef30$AbsEffect.lower[2])|(0>=ef30$AbsEffect.upper[2]))
  sen[1,4] <- as.numeric((0 <= ef50$AbsEffect.lower[2])|(0>=ef50$AbsEffect.upper[2]))
  sen[1,5] <- as.numeric((0 <= ef100$AbsEffect.lower[2])|(0>=ef100$AbsEffect.upper[2]))
  sen[1,6] <- as.numeric((0 <= ef200$AbsEffect.lower[2])|(0>=ef200$AbsEffect.upper[2]))
  sen[1,7] <- as.numeric((0 <= ef300$AbsEffect.lower[2])|(0>=ef300$AbsEffect.upper[2]))
  
  senavg[1,1] <- as.numeric((0 <= ef1$AbsEffect.lower[1])|(0>=ef1$AbsEffect.upper[1]))
  senavg[1,2] <- as.numeric((0 <= ef10$AbsEffect.lower[1])|(0>=ef10$AbsEffect.upper[1]))
  senavg[1,3] <- as.numeric((0 <= ef30$AbsEffect.lower[1])|(0>=ef30$AbsEffect.upper[1]))
  senavg[1,4] <- as.numeric((0 <= ef50$AbsEffect.lower[1])|(0>=ef50$AbsEffect.upper[1]))
  senavg[1,5] <- as.numeric((0 <= ef100$AbsEffect.lower[1])|(0>=ef100$AbsEffect.upper[1]))
  senavg[1,6] <- as.numeric((0 <= ef200$AbsEffect.lower[1])|(0>=ef200$AbsEffect.upper[1]))
  senavg[1,7] <- as.numeric((0 <= ef300$AbsEffect.lower[1])|(0>=ef300$AbsEffect.upper[1]))
  
  avgeffectall[1,1] <- ef1$AbsEffect[1]
  avgeffectall[1,2] <- ef10$AbsEffect[1]
  avgeffectall[1,3] <- ef30$AbsEffect[1]
  avgeffectall[1,4] <- ef50$AbsEffect[1]
  avgeffectall[1,5] <- ef100$AbsEffect[1]
  avgeffectall[1,6] <- ef200$AbsEffect[1]
  avgeffectall[1,7] <- ef300$AbsEffect[1]
  
  cumeffectall[1,1] <- ef1$AbsEffect[2]
  cumeffectall[1,2] <- ef10$AbsEffect[2]
  cumeffectall[1,3] <- ef30$AbsEffect[2]
  cumeffectall[1,4] <- ef50$AbsEffect[2]
  cumeffectall[1,5] <- ef100$AbsEffect[2]
  cumeffectall[1,6] <- ef200$AbsEffect[2]
  cumeffectall[1,7] <- ef300$AbsEffect[2]
  
  avgeffectall.ci.lower[1,1] <- ef1$AbsEffect.lower[1]
  avgeffectall.ci.lower[1,2] <- ef10$AbsEffect.lower[1]
  avgeffectall.ci.lower[1,3] <- ef30$AbsEffect.lower[1]
  avgeffectall.ci.lower[1,4] <- ef50$AbsEffect.lower[1]
  avgeffectall.ci.lower[1,5] <- ef100$AbsEffect.lower[1]
  avgeffectall.ci.lower[1,6] <- ef200$AbsEffect.lower[1]
  avgeffectall.ci.lower[1,7] <- ef300$AbsEffect.lower[1]
  
  avgeffectall.ci.upper[1,1] <- ef1$AbsEffect.upper[1]
  avgeffectall.ci.upper[1,2] <- ef10$AbsEffect.upper[1]
  avgeffectall.ci.upper[1,3] <- ef30$AbsEffect.upper[1]
  avgeffectall.ci.upper[1,4] <- ef50$AbsEffect.upper[1]
  avgeffectall.ci.upper[1,5] <- ef100$AbsEffect.upper[1]
  avgeffectall.ci.upper[1,6] <- ef200$AbsEffect.upper[1]
  avgeffectall.ci.upper[1,7] <- ef300$AbsEffect.upper[1]
  
  cumeffectall.ci.lower[1,1] <- ef1$AbsEffect.lower[2]
  cumeffectall.ci.lower[1,2] <- ef10$AbsEffect.lower[2]
  cumeffectall.ci.lower[1,3] <- ef30$AbsEffect.lower[2]
  cumeffectall.ci.lower[1,4] <- ef50$AbsEffect.lower[2]
  cumeffectall.ci.lower[1,5] <- ef100$AbsEffect.lower[2]
  cumeffectall.ci.lower[1,6] <- ef200$AbsEffect.lower[2]
  cumeffectall.ci.lower[1,7] <- ef300$AbsEffect.lower[2]
  
  cumeffectall.ci.upper[1,1] <- ef1$AbsEffect.upper[2]
  cumeffectall.ci.upper[1,2] <- ef10$AbsEffect.upper[2]
  cumeffectall.ci.upper[1,3] <- ef30$AbsEffect.upper[2]
  cumeffectall.ci.upper[1,4] <- ef50$AbsEffect.upper[2]
  cumeffectall.ci.upper[1,5] <- ef100$AbsEffect.upper[2]
  cumeffectall.ci.upper[1,6] <- ef200$AbsEffect.upper[2]
  cumeffectall.ci.upper[1,7] <- ef300$AbsEffect.upper[2]
  
  p.sen[1,1] <- ef1$p[1]
  p.sen[1,2] <- ef10$p[1]
  p.sen[1,3] <- ef30$p[1]
  p.sen[1,4] <- ef50$p[1]
  p.sen[1,5] <- ef100$p[1]
  p.sen[1,6] <- ef200$p[1]
  p.sen[1,7] <- ef300$p[1]
  
  apee[1,1] <- apee_est(ef1$AbsEffect[1],ci1)
  apee[1,2] <- apee_est(ef10$AbsEffect[1],ci10)
  apee[1,3] <- apee_est(ef30$AbsEffect[1],ci30)
  apee[1,4] <- apee_est(ef50$AbsEffect[1],ci50)
  apee[1,5] <- apee_est(ef100$AbsEffect[1],ci100)
  apee[1,6] <- apee_est(ef200$AbsEffect[1],ci200)
  apee[1,7] <- apee_est(ef300$AbsEffect[1],ci300)
  
  
  
  ef1 <- effect.hsbsts(hsmod$model$bsts.model,UnStandardize=hsmod$UnStandardize,
                       y.model=simemp.mean.y1[1:60],y.cf=simy1[61:100],pre.time=60,post.time=40)
  ef10 <- effect.hsbsts(hsmod$model$bsts.model,UnStandardize=hsmod$UnStandardize,
                        y.model=simemp.mean.y10[1:60],y.cf=simy10[61:100],pre.time=60,post.time=40)    
  ef30 <- effect.hsbsts(hsmod$model$bsts.model,UnStandardize=hsmod$UnStandardize,
                        y.model=simemp.mean.y30[1:60],y.cf=simy30[61:100],pre.time=60,post.time=40)
  ef50 <- effect.hsbsts(hsmod$model$bsts.model,UnStandardize=hsmod$UnStandardize,
                        y.model=simemp.mean.y50[1:60],y.cf=simy50[61:100],pre.time=60,post.time=40)
  ef100 <- effect.hsbsts(hsmod$model$bsts.model,UnStandardize=hsmod$UnStandardize,
                         y.model=simemp.mean.y100[1:60],y.cf=simy100[61:100],pre.time=60,post.time=40)
  ef200 <- effect.hsbsts(hsmod$model$bsts.model,UnStandardize=hsmod$UnStandardize,
                         y.model=simemp.mean.y200[1:60],y.cf=simy200[61:100],pre.time=60,post.time=40)
  ef300 <- effect.hsbsts(hsmod$model$bsts.model,UnStandardize=hsmod$UnStandardize,
                         y.model=simemp.mean.y300[1:60],y.cf=simy300[61:100],pre.time=60,post.time=40)
  
  
  avgeffectall[2,1] <- ef1$AbsEffect[1]
  avgeffectall[2,2] <- ef10$AbsEffect[1]
  avgeffectall[2,3] <- ef30$AbsEffect[1]
  avgeffectall[2,4] <- ef50$AbsEffect[1]
  avgeffectall[2,5] <- ef100$AbsEffect[1]
  avgeffectall[2,6] <- ef200$AbsEffect[1]
  avgeffectall[2,7] <- ef300$AbsEffect[1]
  
  cumeffectall[2,1] <- ef1$AbsEffect[2]
  cumeffectall[2,2] <- ef10$AbsEffect[2]
  cumeffectall[2,3] <- ef30$AbsEffect[2]
  cumeffectall[2,4] <- ef50$AbsEffect[2]
  cumeffectall[2,5] <- ef100$AbsEffect[2]
  cumeffectall[2,6] <- ef200$AbsEffect[2]
  cumeffectall[2,7] <- ef300$AbsEffect[2]
  
  avgeffectall.ci.lower[2,1] <- ef1$AbsEffect.lower[1]
  avgeffectall.ci.lower[2,2] <- ef10$AbsEffect.lower[1]
  avgeffectall.ci.lower[2,3] <- ef30$AbsEffect.lower[1]
  avgeffectall.ci.lower[2,4] <- ef50$AbsEffect.lower[1]
  avgeffectall.ci.lower[2,5] <- ef100$AbsEffect.lower[1]
  avgeffectall.ci.lower[2,6] <- ef200$AbsEffect.lower[1]
  avgeffectall.ci.lower[2,7] <- ef300$AbsEffect.lower[1]
  
  avgeffectall.ci.upper[2,1] <- ef1$AbsEffect.upper[1]
  avgeffectall.ci.upper[2,2] <- ef10$AbsEffect.upper[1]
  avgeffectall.ci.upper[2,3] <- ef30$AbsEffect.upper[1]
  avgeffectall.ci.upper[2,4] <- ef50$AbsEffect.upper[1]
  avgeffectall.ci.upper[2,5] <- ef100$AbsEffect.upper[1]
  avgeffectall.ci.upper[2,6] <- ef200$AbsEffect.upper[1]
  avgeffectall.ci.upper[2,7] <- ef300$AbsEffect.upper[1]
  
  cumeffectall.ci.lower[2,1] <- ef1$AbsEffect.lower[2]
  cumeffectall.ci.lower[2,2] <- ef10$AbsEffect.lower[2]
  cumeffectall.ci.lower[2,3] <- ef30$AbsEffect.lower[2]
  cumeffectall.ci.lower[2,4] <- ef50$AbsEffect.lower[2]
  cumeffectall.ci.lower[2,5] <- ef100$AbsEffect.lower[2]
  cumeffectall.ci.lower[2,6] <- ef200$AbsEffect.lower[2]
  cumeffectall.ci.lower[2,7] <- ef300$AbsEffect.lower[2]
  
  cumeffectall.ci.upper[2,1] <- ef1$AbsEffect.upper[2]
  cumeffectall.ci.upper[2,2] <- ef10$AbsEffect.upper[2]
  cumeffectall.ci.upper[2,3] <- ef30$AbsEffect.upper[2]
  cumeffectall.ci.upper[2,4] <- ef50$AbsEffect.upper[2]
  cumeffectall.ci.upper[2,5] <- ef100$AbsEffect.upper[2]
  cumeffectall.ci.upper[2,6] <- ef200$AbsEffect.upper[2]
  cumeffectall.ci.upper[2,7] <- ef300$AbsEffect.upper[2]
  
  senavg[2,1] <- as.numeric((0 <= ef1$AbsEffect.lower[1])|(0>=ef1$AbsEffect.upper[1]))
  senavg[2,2] <- as.numeric((0 <= ef10$AbsEffect.lower[1])|(0>=ef10$AbsEffect.upper[1]))
  senavg[2,3] <- as.numeric((0 <= ef30$AbsEffect.lower[1])|(0>=ef30$AbsEffect.upper[1]))
  senavg[2,4] <- as.numeric((0 <= ef50$AbsEffect.lower[1])|(0>=ef50$AbsEffect.upper[1]))
  senavg[2,5] <- as.numeric((0 <= ef100$AbsEffect.lower[1])|(0>=ef100$AbsEffect.upper[1]))
  senavg[2,6] <- as.numeric((0 <= ef200$AbsEffect.lower[1])|(0>=ef200$AbsEffect.upper[1]))
  senavg[2,7] <- as.numeric((0 <= ef300$AbsEffect.lower[1])|(0>=ef300$AbsEffect.upper[1]))
  
  sen[2,1] <- as.numeric((0 <= ef1$AbsEffect.lower[2])|(0>=ef1$AbsEffect.upper[2]))
  sen[2,2] <- as.numeric((0 <= ef10$AbsEffect.lower[2])|(0>=ef10$AbsEffect.upper[2]))
  sen[2,3] <- as.numeric((0 <= ef30$AbsEffect.lower[2])|(0>=ef30$AbsEffect.upper[2]))
  sen[2,4] <- as.numeric((0 <= ef50$AbsEffect.lower[2])|(0>=ef50$AbsEffect.upper[2]))
  sen[2,5] <- as.numeric((0 <= ef100$AbsEffect.lower[2])|(0>=ef100$AbsEffect.upper[2]))
  sen[2,6] <- as.numeric((0 <= ef200$AbsEffect.lower[2])|(0>=ef200$AbsEffect.upper[2]))
  sen[2,7] <- as.numeric((0 <= ef300$AbsEffect.lower[2])|(0>=ef300$AbsEffect.upper[2]))
  
  p.sen[2,1] <- ef1$p[1]
  p.sen[2,2] <- ef10$p[1]
  p.sen[2,3] <- ef30$p[1]
  p.sen[2,4] <- ef50$p[1]
  p.sen[2,5] <- ef100$p[1]
  p.sen[2,6] <- ef200$p[1]
  p.sen[2,7] <- ef300$p[1]
  
  apee[2,1] <- apee_est(ef1$AbsEffect[1],ci1)
  apee[2,2] <- apee_est(ef10$AbsEffect[1],ci10)
  apee[2,3] <- apee_est(ef30$AbsEffect[1],ci30)
  apee[2,4] <- apee_est(ef50$AbsEffect[1],ci50)
  apee[2,5] <- apee_est(ef100$AbsEffect[1],ci100)
  apee[2,6] <- apee_est(ef200$AbsEffect[1],ci200)
  apee[2,7] <- apee_est(ef300$AbsEffect[1],ci300)
  
  ef1 <- effect.ssbsts(ssmod$model$bsts.model,y.cf=simy1[61:100],post.period=c(61,100),UnStandardize = ssmod$UnStandardize)
  ef10 <- effect.ssbsts(ssmod$model$bsts.model,y.cf=simy10[61:100],post.period=c(61,100),UnStandardize = ssmod$UnStandardize)
  ef30 <- effect.ssbsts(ssmod$model$bsts.model,y.cf=simy30[61:100],post.period=c(61,100),UnStandardize = ssmod$UnStandardize)
  ef50 <- effect.ssbsts(ssmod$model$bsts.model,y.cf=simy50[61:100],post.period=c(61,100),UnStandardize = ssmod$UnStandardize)
  ef100 <- effect.ssbsts(ssmod$model$bsts.model,y.cf=simy100[61:100],post.period=c(61,100),UnStandardize = ssmod$UnStandardize)
  ef200 <- effect.ssbsts(ssmod$model$bsts.model,y.cf=simy200[61:100],post.period=c(61,100),UnStandardize = ssmod$UnStandardize)
  ef300 <- effect.ssbsts(ssmod$model$bsts.model,y.cf=simy300[61:100],post.period=c(61,100),UnStandardize = ssmod$UnStandardize)
  
  
  
  avgeffectall[3,1] <- ef1$AbsEffect[1]
  avgeffectall[3,2] <- ef10$AbsEffect[1]
  avgeffectall[3,3] <- ef30$AbsEffect[1]
  avgeffectall[3,4] <- ef50$AbsEffect[1]
  avgeffectall[3,5] <- ef100$AbsEffect[1]
  avgeffectall[3,6] <- ef200$AbsEffect[1]
  avgeffectall[3,7] <- ef300$AbsEffect[1]
  
  cumeffectall[3,1] <- ef1$AbsEffect[2]
  cumeffectall[3,2] <- ef10$AbsEffect[2]
  cumeffectall[3,3] <- ef30$AbsEffect[2]
  cumeffectall[3,4] <- ef50$AbsEffect[2]
  cumeffectall[3,5] <- ef100$AbsEffect[2]
  cumeffectall[3,6] <- ef200$AbsEffect[2]
  cumeffectall[3,7] <- ef300$AbsEffect[2]
  
  avgeffectall.ci.lower[3,1] <- ef1$AbsEffect.lower[1]
  avgeffectall.ci.lower[3,2] <- ef10$AbsEffect.lower[1]
  avgeffectall.ci.lower[3,3] <- ef30$AbsEffect.lower[1]
  avgeffectall.ci.lower[3,4] <- ef50$AbsEffect.lower[1]
  avgeffectall.ci.lower[3,5] <- ef100$AbsEffect.lower[1]
  avgeffectall.ci.lower[3,6] <- ef200$AbsEffect.lower[1]
  avgeffectall.ci.lower[3,7] <- ef300$AbsEffect.lower[1]
  
  avgeffectall.ci.upper[3,1] <- ef1$AbsEffect.upper[1]
  avgeffectall.ci.upper[3,2] <- ef10$AbsEffect.upper[1]
  avgeffectall.ci.upper[3,3] <- ef30$AbsEffect.upper[1]
  avgeffectall.ci.upper[3,4] <- ef50$AbsEffect.upper[1]
  avgeffectall.ci.upper[3,5] <- ef100$AbsEffect.upper[1]
  avgeffectall.ci.upper[3,6] <- ef200$AbsEffect.upper[1]
  avgeffectall.ci.upper[3,7] <- ef300$AbsEffect.upper[1]
  
  cumeffectall.ci.lower[3,1] <- ef1$AbsEffect.lower[2]
  cumeffectall.ci.lower[3,2] <- ef10$AbsEffect.lower[2]
  cumeffectall.ci.lower[3,3] <- ef30$AbsEffect.lower[2]
  cumeffectall.ci.lower[3,4] <- ef50$AbsEffect.lower[2]
  cumeffectall.ci.lower[3,5] <- ef100$AbsEffect.lower[2]
  cumeffectall.ci.lower[3,6] <- ef200$AbsEffect.lower[2]
  cumeffectall.ci.lower[3,7] <- ef300$AbsEffect.lower[2]
  
  cumeffectall.ci.upper[3,1] <- ef1$AbsEffect.upper[2]
  cumeffectall.ci.upper[3,2] <- ef10$AbsEffect.upper[2]
  cumeffectall.ci.upper[3,3] <- ef30$AbsEffect.upper[2]
  cumeffectall.ci.upper[3,4] <- ef50$AbsEffect.upper[2]
  cumeffectall.ci.upper[3,5] <- ef100$AbsEffect.upper[2]
  cumeffectall.ci.upper[3,6] <- ef200$AbsEffect.upper[2]
  cumeffectall.ci.upper[3,7] <- ef300$AbsEffect.upper[2]
  
  senavg[3,1] <- as.numeric((0 <= ef1$AbsEffect.lower[1])|(0>=ef1$AbsEffect.upper[1]))
  senavg[3,2] <- as.numeric((0 <= ef10$AbsEffect.lower[1])|(0>=ef10$AbsEffect.upper[1]))
  senavg[3,3] <- as.numeric((0 <= ef30$AbsEffect.lower[1])|(0>=ef30$AbsEffect.upper[1]))
  senavg[3,4] <- as.numeric((0 <= ef50$AbsEffect.lower[1])|(0>=ef50$AbsEffect.upper[1]))
  senavg[3,5] <- as.numeric((0 <= ef100$AbsEffect.lower[1])|(0>=ef100$AbsEffect.upper[1]))
  senavg[3,6] <- as.numeric((0 <= ef200$AbsEffect.lower[1])|(0>=ef200$AbsEffect.upper[1]))
  senavg[3,7] <- as.numeric((0 <= ef300$AbsEffect.lower[1])|(0>=ef300$AbsEffect.upper[1]))
  
  sen[3,1] <- as.numeric((0 <= ef1$AbsEffect.lower[2])|(0>=ef1$AbsEffect.upper[2]))
  sen[3,2] <- as.numeric((0 <= ef10$AbsEffect.lower[2])|(0>=ef10$AbsEffect.upper[2]))
  sen[3,3] <- as.numeric((0 <= ef30$AbsEffect.lower[2])|(0>=ef30$AbsEffect.upper[2]))
  sen[3,4] <- as.numeric((0 <= ef50$AbsEffect.lower[2])|(0>=ef50$AbsEffect.upper[2]))
  sen[3,5] <- as.numeric((0 <= ef100$AbsEffect.lower[2])|(0>=ef100$AbsEffect.upper[2]))
  sen[3,6] <- as.numeric((0 <= ef200$AbsEffect.lower[2])|(0>=ef200$AbsEffect.upper[2]))
  sen[3,7] <- as.numeric((0 <= ef300$AbsEffect.lower[2])|(0>=ef300$AbsEffect.upper[2]))
  
  p.sen[3,1] <- ef1$p[1]
  p.sen[3,2] <- ef10$p[1]
  p.sen[3,3] <- ef30$p[1]
  p.sen[3,4] <- ef50$p[1]
  p.sen[3,5] <- ef100$p[1]
  p.sen[3,6] <- ef200$p[1]
  p.sen[3,7] <- ef300$p[1]
  
  apee[3,1] <- apee_est(ef1$AbsEffect[1],ci1)
  apee[3,2] <- apee_est(ef10$AbsEffect[1],ci10)
  apee[3,3] <- apee_est(ef30$AbsEffect[1],ci30)
  apee[3,4] <- apee_est(ef50$AbsEffect[1],ci50)
  apee[3,5] <- apee_est(ef100$AbsEffect[1],ci100)
  apee[3,6] <- apee_est(ef200$AbsEffect[1],ci200)
  apee[3,7] <- apee_est(ef300$AbsEffect[1],ci300)
  
  
  return(list(pred.error = pred.error,beta.error=beta.error,inc.error=inc.error,p = p,p.sen=p.sen,apee=apee,
              effect=effect, spec=spec, sen=sen,cumeffect=cumeffect,
              specavg=specavg,cumeffectall=cumeffectall,avgeffectall=avgeffectall,senavg=senavg,
              effect.ci=effect.ci,cumeffect.ci=cumeffect.ci,
              avgeffectall.ci.lower=avgeffectall.ci.lower,avgeffectall.ci.upper=avgeffectall.ci.upper,
              cumeffectall.ci.lower=cumeffectall.ci.lower,cumeffectall.ci.upper=cumeffectall.ci.upper,
              realeffect=realeffect))
}




M=20
#realbeta <- c(0.5,-0.1,2,0.5,rep(0,2))
#realbeta <- c(0.1,0.5,0.5,0.5,rep(0,17))
#realbeta <- c(0.5,-0.1,2,0.5,0.3,1)
#realbeta <- c(0.1,-0.5,1,0.1,0.5,1)
#realbeta <- c(0.1,-0.5,1,0.1,rep(0,2))
#realbeta <- rep(1/20,21)#c(rep(1/3,4),0,0)#c(-0.0015,0.1,0.8,-0.06)#c(0.1,-0.5,1,0.5,rep(0,2))
meff <- 3



pred.error <- beta.error  <- inc.error <-
  effect <- cumeffect <- p <- spec <- specavg <- rep(0,3)
effect.ci <- cumeffect.ci <- matrix(0,nrow=3,ncol=2)
sen <- p.sen <- apee <- cumeffectall<-avgeffectall<- senavg <- avgeffectall.ci.lower <- avgeffectall.ci.upper<-cumeffectall.ci.lower <- cumeffectall.ci.upper <- matrix(0,nrow=3,ncol=7)
realeffect <- matrix(0,nrow=2,ncol=7)

numCores <- detectCores()
miss=1
results <- mclapply(rep(1,50), sim,mc.cores = numCores)




