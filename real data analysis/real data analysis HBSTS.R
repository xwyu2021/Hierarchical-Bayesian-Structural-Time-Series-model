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


policydate <- as.Date("2014-05-14")
enddate <- as.Date("2017-11-27") # before the second intervention

white_groups <- 1:4
mixed_groups <- 7:8
treat_groups <- c(14)
ethn_treat <- 14
donor_pool <- c(1:17,97)[-which(c(1:17,97) %in% c(treat_groups))]
nonwhite <- donor_pool[-which(donor_pool %in% white_groups)]
ethn.names <- c("Intercept","British","Irish","Gypsy/Other White",
                "Mixed Asian","Other Mixed","Indian","Pakistani",
                "Bangaldesh","Chinese","Other Asian","Black African",
                "Other black","Arab","Other")

#######################################################################
####################### load the processed data #######################
#######################################################################

## original data from https://www.understandingsociety.ac.uk/
## missing data imputation following Jeffery et al. 2024
## data1 is the processed data

data1 <- readRDS('finaldata_weight.RData')
survey.sub <-  data1 %>%
  dplyr::rename(c('id' = 'pidp', 'edu' = 'hiqual', 'rela' = 'mastat', 'jobStatus' = 'jbstat', 'ghq' = 'scghq1')) %>%
  dplyr::filter(interviewDate <= enddate)


length(unique(survey.sub$id)) # 59321
survey.sub$policyDate <- policydate
eth.df <- survey.sub %>% 
  dplyr::select(c('interviewDate','eth1','ghq','id','policyDate','wave','strata','psu'))
eth.df <- eth.df %>%
  dplyr::mutate(month = (lubridate::interval(start = policyDate, end = interviewDate) %/% months(1)))
eth.df <- eth.df %>%
  dplyr::mutate(season = (lubridate::interval(start = policyDate, end = interviewDate) %/% months(3)))

eth.df$strata = as.character(eth.df$strata)
eth.df$strata = as.factor(eth.df$strata)
eth.df$strata = as.numeric(eth.df$strata)
eth.df$psu = as.character(eth.df$psu)
eth.df$psu = as.factor(eth.df$psu)
eth.df$psu = as.numeric(eth.df$psu)
length(unique(eth.df$strata))
range(unique(eth.df$strata))

# -- extract treated Time Series --
treated.whole.df <- eth.df %>%
  dplyr::filter(eth1 %in% ethn_treat) %>%
  dplyr::mutate(season = (lubridate::interval(start = policyDate, end = interviewDate) %/% months(3)))
paste('there are',length(unique(treated.whole.df$id)),'individuals in the treated')
paste('there are',length(unique(treated.whole.df$interviewDate)),'dates in the treated')
paste('there are',length(unique(treated.whole.df$month)),'months in the treated')
paste('there are',length(unique(treated.whole.df$season)),'seasons in the treated')
treated.df.qua <- aggregate(treated.whole.df$ghq,by=list(treated.whole.df$season),mean)
colnames(treated.df.qua) <- c('season','ghq')

# -- extract control Time Series --
control.whole.df <- eth.df %>%
  dplyr::filter(eth1 %in% donor_pool)%>%
  dplyr::mutate(season = (lubridate::interval(start = policyDate, end = interviewDate) %/% months(3)))

control.df.qua <- dcast(control.whole.df,season ~ eth1, value.var='ghq',function(x) mean(x[!is.na(x)]))
paste(colSums(is.na(control.df.qua))*100/dim(control.df.qua)[1],'% missing values for control group', colnames(control.df.qua))
# merge ethnicity 3 to ethnicity 4
control.whole.df$new_ethn <- control.whole.df$eth1
control.whole.df$new_ethn[control.whole.df$new_ethn == 3] <- 4
control.df.qua <- dcast(control.whole.df,season ~ new_ethn, value.var='ghq',function(x) mean(x[!is.na(x)]))
paste(colSums(is.na(control.df.qua))*100/dim(control.df.qua)[1],'% missing values for control group', colnames(control.df.qua))


#######################################################################
################ Applying BSTS (Brodersen et al. 2015) ################
#######################################################################

bsts.df.qua <- merge(treated.df.qua, control.df.qua)
treated.zoo.qua <- read.zoo(treated.df.qua, index_column=1)
bsts.zoo.qua <-read.zoo(bsts.df.qua, index_column=1)
pre.period.qua <- c(min(index(bsts.zoo.qua)),0)
post.period.qua <- c(1,max(index(bsts.zoo.qua)))
post.y.qua <- post.period.response.qua <- as.vector(bsts.zoo.qua$ghq[index(bsts.zoo.qua) >0])
pre.time.qua <- length(c(-21:0))
post.time.qua <- length(c(1:14))

# --  Causalimpact package --
hold.qua <- bsts.zoo.qua
hold.qua[index(bsts.zoo.qua) >0,1] <- NA
ci <- CausalImpact(bsts.zoo.qua,pre.period = pre.period.qua,post.period = post.period.qua,
                   model.args = list(nseasons=1, prior.level.sd = 0.01))
plot(ci) 

# -- use our adapted code --
ce<- run_flex_bsts(data=bsts.zoo.qua,pre.period=pre.period.qua,
                   post.period=post.period.qua,
                   alpha=0.05,model.args=list(nseasons=1,
                                              prior.level.sd = 0.01,
                                              standardize.data=TRUE),
                   trend = 'AR1',modelsize=3,cycles = c(1),niter=10000)
plot(ce)

###################################################################################################
################ Applying HS-BSTS (standard Bsts with regularised horseshoe prior) ################
###################################################################################################

ce.rhs<- run_flex_bsts(data=bsts.white.qua,pre.period=pre.period.qua,
                   post.period=post.period.qua,
                   alpha=0.05,model.args=list(nseasons=4,
                                              prior.level.sd = 0.01,
                                              standardize.data=TRUE),
                   trend = 'AR1',modelsize=3,cycles = c(1),niter=10000)
plot(ce.rhs)


################################################
################ Applying HBSTS ################
################################################



pre.time.qua <- 22
post.time.qua <- 14
treated.multi <- treated.whole.df %>% dplyr::select(season,ghq,strata,psu,wave) %>% dplyr::arrange(season)
all.multi <- merge(treated.multi, control.df.qua, by='season')
colnames(all.multi) <- c('season','carib','British','Irish','Other White',
                         'mixed asian',
                         'other mixed','Indian',
                         'pakistani','Bangaldeshi','Chinese',
                         'other asian','African','other black','arab',
                         'Other')

pre.ys.qua <- as.vector(treated.multi$ghq[treated.multi$season <=0])
post.ys.qua <- as.vector(treated.multi$ghq[treated.multi$season >0])
pre.obs.qua <- sum(treated.multi$season <= 0)
post.obs.qua <- sum(treated.multi$season > 0)
pre.strata.qua <- as.vector(treated.multi$strata[treated.multi$season <=0])
post.strata.qua <- as.vector(treated.multi$strata[treated.multi$season >0])
obs.count <- as.vector(table(treated.multi$season))
index_pre <- cumsum(c(1,obs.count[1:pre.time.qua]))
index_post <- cumsum(c(1,obs.count[(pre.time.qua+1):(pre.time.qua+post.time.qua)]))

# ---- control group: White only -----
control.whole.df <- eth.df %>%
  dplyr::filter(eth1 %in% white_groups)%>%
  dplyr::mutate(season = (lubridate::interval(start = policyDate, end = interviewDate) %/% months(3)))
x.obs <- control.whole.df %>% dplyr::select(ghq,eth1,season)
x.obs <- x.obs %>% dplyr::arrange(season,eth1)
x.obs <- x.obs %>% dplyr::mutate(eth = dplyr::case_when(eth1 == 1 ~ 1,
                                                        eth1 == 2 ~ 2,
                                                        eth1 == 3 ~ 3,
                                                        eth1 == 4 ~ 3))
x.obs <- x.obs[,c(1,3,4)]
x.obs <- as.matrix(x.obs)
check <- as.data.frame(x.obs[,2:3])
x.count <- table(check)
x.cum.count <- t(matrix(cumsum(t(x.count)),nrow=3,ncol=pre.time.qua + post.time.qua))
x.obs <- as.vector(x.obs[,1])
K <- length(x.obs)
emp.mean.y <- as.vector(bsts.zoo.qua[,1])

pre.strata.qua <- as.vector(treated.multi$strata[treated.multi$season <=0])
post.strata.qua <- as.vector(treated.multi$strata[treated.multi$season >0])
x.strata <- control.whole.df %>% dplyr::select(strata,eth1,season)
x.strata <- x.strata %>% dplyr::arrange(season,eth1)
x.strata <- as.vector(x.strata[,c(1)])
NS <- length(unique(eth.df$strata))

pre.psu.qua <- as.vector(treated.multi$psu[treated.multi$season <=0])
post.psu.qua <- as.vector(treated.multi$psu[treated.multi$season >0])
x.psu <- control.whole.df %>% dplyr::select(psu,eth1,season)
x.psu <- x.psu%>% dplyr::arrange(season,eth1)
x.psu <- as.vector(x.psu[,c(1)])

allpsu <- c(pre.psu.qua,x.psu)
allpsu <- as.factor(allpsu)
allpsu <- as.numeric(allpsu)
pre.psu.qua <- allpsu[1:length(pre.psu.qua)]
x.psu <- allpsu[(length(pre.psu.qua)+1): (length(pre.psu.qua) + length(x.psu))]
NC <- length(unique(allpsu))

hmodmw <- run_hbsts_re(pre.time=pre.time.qua,post.time=post.time.qua,
                        pre.obs=pre.obs.qua,
                        pre.ys=pre.ys.qua,dimx=3,
                        post.obs =post.obs.qua,index_pre=index_pre,
                        K=K,x.obs=x.obs, x.cum.count=x.cum.count,index_post=index_post,
                        niter=10000,nchains=1,time_index=-21:14,
                        alpha=0.05,cycles=c(1,2,4,8),emp.mean.y=emp.mean.y,
                        trend = 'AR1',burn=1000,modelsize=2,adapt=0.9,
                        x_strata=x.strata,NS=NS,pre_strata = pre.strata.qua,
                        x_psu=x.psu,NC=NC,pre_psu = pre.psu.qua,ind=FALSE,hdi=FALSE)
hmodmw$summary
plot(hmodmw)

# ---- control group: all other 14 groups -----

control.whole.df <- eth.df %>%
  dplyr::filter(eth1 %in% donor_pool)%>%
  dplyr::mutate(season = (lubridate::interval(start = policyDate, end = interviewDate) %/% months(3)))
x.obs <- control.whole.df %>% dplyr::select(ghq,eth1,season)
x.obs <- x.obs %>% dplyr::arrange(season,eth1)
x.obs <- x.obs %>% dplyr::mutate(eth = dplyr::case_when(eth1 == 1 ~ 1,
                                                        eth1 == 2 ~ 2,
                                                        eth1 == 3 ~ 3,
                                                        eth1 == 4 ~ 3,
                                                        eth1 == 7 ~ 4, ## mix asian
                                                        eth1 == 8 ~ 5,
                                                        eth1 == 9 ~ 6,
                                                        eth1 == 10 ~ 7,## paki
                                                        eth1 == 11 ~ 8,
                                                        eth1 == 12 ~ 9,
                                                        eth1 == 13 ~ 10,
                                                        eth1 == 15 ~ 11, ## black african
                                                        eth1 == 16 ~ 12,
                                                        eth1 == 17 ~ 13,
                                                        eth1 == 97 ~ 14))

x.obs <- x.obs[,c(1,3,4)]
x.obs <- as.matrix(x.obs)
check <- as.data.frame(x.obs[,2:3])
x.count <- table(check)
x.cum.count <- t(matrix(cumsum(t(x.count)),nrow=14,ncol=pre.time.qua + post.time.qua))
x.obs <- as.vector(x.obs[,1])
K <- length(x.obs)

pre.strata.qua <- as.vector(treated.multi$strata[treated.multi$season <=0])
post.strata.qua <- as.vector(treated.multi$strata[treated.multi$season >0])
x.strata <- control.whole.df %>% dplyr::select(strata,eth1,season)
x.strata <- x.strata %>% dplyr::arrange(season,eth1)
x.strata <- as.vector(x.strata[,c(1)])
NS <- length(unique(eth.df$strata))

pre.psu.qua <- as.vector(treated.multi$psu[treated.multi$season <=0])
post.psu.qua <- as.vector(treated.multi$psu[treated.multi$season >0])
x.psu <- control.whole.df %>% dplyr::select(psu,eth1,season)
x.psu <- x.psu%>% dplyr::arrange(season,eth1)
x.psu <- as.vector(x.psu[,c(1)])

allpsu <- c(pre.psu.qua,x.psu)
allpsu <- as.factor(allpsu)
allpsu <- as.numeric(allpsu)
pre.psu.qua <- allpsu[1:length(pre.psu.qua)]
x.psu <- allpsu[(length(pre.psu.qua)+1): (length(pre.psu.qua) + length(x.psu))]
NC <- length(unique(allpsu))

emp.mean.y = as.vector(bsts.df.qua$ghq)
hmod <- run_hbsts_re(pre.time=pre.time.qua,post.time=post.time.qua,
                     pre.obs=pre.obs.qua,
                     pre.ys=pre.ys.qua,dimx=14,
                     post.obs =post.obs.qua,index_pre=index_pre,
                     K=K,x.obs=x.obs, x.cum.count=x.cum.count,index_post=index_post,
                     niter=8000,nchains=1,time_index=-21:14,
                     alpha=0.05,cycles=c(1,2,4,8),emp.mean.y=emp.mean.y,
                     trend = 'AR1',burn=3000,modelsize=3,adapt=0.9,
                     x_strata=x.strata,NS=NS,pre_strata = pre.strata.qua,
                     x_psu=x.psu,NC=NC,pre_psu = pre.psu.qua,ind=FALSE,hdi=FALSE)


hmod$summary
plot(hmod)

###########################################
################### plots #################
###########################################

## -- pointwise effects by quarter, half year, and year -- ##
# - quarter
quas <- seq(as.Date('2014-05-01'),as.Date('2017-11-01'), by = "3 months")
quas <- quas[-length(quas)]
postqua <- 14

pmat <- cbind(tail(hmodmw$series$point.effect.lower,postqua),tail(hmodmw$series$point.effect.upper,postqua)) 
ps <- apply(pmat,1,function(x) (x[1]<0 & x[2] <0) | (x[1] >0 & x[2] > 0))
efdat <- data.frame(effect = tail(hmodmw$series$point.effect,postqua),
                    upper = tail(hmodmw$series$point.effect.upper,postqua),
                    lower = tail(hmodmw$series$point.effect.lower,postqua),
                    p = ps,
                    time = quas)
efdat$significant <- as.factor(efdat$p)
ggplot(aes(x=time,y=effect),data=efdat) + 
  geom_point() + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(0.05))+
  theme_bw() + theme(axis.title = element_text(size = 18),
                     axis.text= element_text(size=18),
                     legend.text=element_text(size=13),
                     legend.title = element_text(size=13),
                     legend.position = 'none',
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     panel.grid.minor = element_blank())+
  geom_point(aes(shape=significant,colour = significant),size=3.5) +
  scale_color_manual(values=c('TRUE'='orange','FALSE'='grey2')) +labs(x='time (years)')

# - half year
half <- seq(as.Date('2014-05-01'),as.Date('2017-11-01'), by = "6 months")
half <- half[-length(half)]
pmat <- cbind(hmodmw$temporaleffect$temp_effect[[2]]$tempeffect.lower,hmodmw$temporaleffect$temp_effect[[2]]$tempeffect.upper)
ps <- apply(pmat,1,function(x) (x[1]<0 & x[2] <0) | (x[1] >0 & x[2] > 0))
efdat <- data.frame(effect = hmodmw$temporaleffect$temp_effect[[2]]$tempeffect,
                    upper = hmodmw$temporaleffect$temp_effect[[2]]$tempeffect.upper,
                    lower = hmodmw$temporaleffect$temp_effect[[2]]$tempeffect.lower,
                    p = ps,
                    time = half)
efdat$significant <- as.factor(efdat$p)

ggplot(aes(x=time,y=effect),data=efdat) + 
  geom_point() + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(0.05))+
  theme_bw() + theme(axis.title = element_text(size = 18),
                     axis.text= element_text(size=18),
                     legend.text=element_text(size=13),
                     legend.title = element_text(size=13),
                     legend.position = 'none',
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     panel.grid.minor = element_blank())+
  geom_point(aes(shape=significant,colour = significant),size=3.5) +
  scale_color_manual(values=c('TRUE'='orange','FALSE'='grey2')) + labs(x='time (years)')

# - year

year <- seq(as.Date('2014-05-14'),as.Date('2018-05-14'), by = "12 months")
year <- year[-length(year)]
pmat <- cbind(hmodmw$temporaleffect$temp_effect[[3]]$tempeffect.lower,hmodmw$temporaleffect$temp_effect[[3]]$tempeffect.upper)
ps <- apply(pmat,1,function(x) (x[1]<0 & x[2] <0) | (x[1] >0 & x[2] > 0))
efdat <- data.frame(effect = hmodmw$temporaleffect$temp_effect[[3]]$tempeffect,
                    upper = hmodmw$temporaleffect$temp_effect[[3]]$tempeffect.upper,
                    lower = hmodmw$temporaleffect$temp_effect[[3]]$tempeffect.lower,
                    p = ps,
                    time = year)
efdat$significant <- as.factor(efdat$p)
ggplot(aes(x=time,y=effect),data=efdat) + 
  geom_point() + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(0.05))+
  theme_bw() + theme(axis.title = element_text(size = 18),
                     axis.text= element_text(size=18),
                     legend.text=element_text(size=13),
                     legend.title = element_text(size=13),
                     legend.position = 'none',
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     panel.grid.minor = element_blank())+
  geom_point(aes(shape=significant,colour = significant),size=3.5) +
  scale_color_manual(values=c('TRUE'='orange','FALSE'='grey2')) + labs(x='time (years)')



## -- dynamic contribution of different covariates (groups) -- ##
fit <- rstan::extract(hmodmw$model$bsts.model)
PlotDynamicDistribution(curves=cbind(fit$mu,fit$mu_forecast),timestamps = -21:14,
                        xlab='time (quarters)',cex.lab=1.2,cex.axis=1.5)
PlotDynamicDistribution(curves=cbind(fit$f,fit$f_forecast),timestamps = -21:14,
                        xlab='time (quarters)',cex.lab=1.2,cex.axis=1.5)

PlotDynamicDistribution(curves=matrix(rep(fit$beta[,1],36),nrow=9000) *fit$x_scale[,,1],
                        timestamps = -21:14,xlab='time (quaters)',cex.lab=1.2,cex.axis=1.5)
lines(-21:14,colMeans(matrix(rep(fit$beta[,1],36),nrow=9000) *fit$x_scale[,,1]),
      col=2)


PlotDynamicDistribution(curves=matrix(rep(fit$beta[,2],36),nrow=9000) *fit$x_scale[,,2],
                        timestamps = -21:14,xlab='time (quaters)',cex.lab=1.2,cex.axis=1.5)
lines(-21:14,colMeans(matrix(rep(fit$beta[,2],36),nrow=9000) *fit$x_scale[,,2]),
      col=2)

PlotDynamicDistribution(curves=matrix(rep(fit$beta[,3],36),nrow=9000) *fit$x_scale[,,3],
                        timestamps = -21:14,xlab='time (quaters)',cex.lab=1.2,cex.axis=1.5)
lines(-21:14,colMeans(matrix(rep(fit$beta[,3],36),nrow=9000) *fit$x_scale[,,3]),
      col=2)










