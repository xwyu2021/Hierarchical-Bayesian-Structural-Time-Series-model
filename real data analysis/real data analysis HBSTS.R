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
library(data.table)


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
  dplyr::mutate(season = (lubridate::interval(start = policyDate, end = interviewDate) %/% months(1)))

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
  dplyr::mutate(season = (lubridate::interval(start = policyDate, end = interviewDate) %/% months(1)))
paste('there are',length(unique(treated.whole.df$id)),'individuals in the treated')
paste('there are',length(unique(treated.whole.df$interviewDate)),'dates in the treated')
paste('there are',length(unique(treated.whole.df$month)),'months in the treated')
paste('there are',length(unique(treated.whole.df$season)),'seasons in the treated')
treated.df.qua <- aggregate(treated.whole.df$ghq,by=list(treated.whole.df$season),mean)
colnames(treated.df.qua) <- c('season','ghq')

# -- extract control Time Series --
control.whole.df <- eth.df %>%
  dplyr::filter(eth1 %in% donor_pool)%>%
  dplyr::mutate(season = (lubridate::interval(start = policyDate, end = interviewDate) %/% months(1)))

control.df.qua <- dcast(as.data.table(control.whole.df),
                        season ~ eth1, value.var='ghq',function(x) mean(x[!is.na(x)]))
paste(colSums(is.na(control.df.qua))*100/dim(control.df.qua)[1],'% missing values for control group', colnames(control.df.qua))
# merge ethnicity 3 to ethnicity 4
control.whole.df$new_ethn <- control.whole.df$eth1
control.whole.df$new_ethn[control.whole.df$new_ethn == 3] <- 4
control.df.qua <- dcast(as.data.table(control.whole.df),season ~ new_ethn, value.var='ghq',function(x) mean(x[!is.na(x)]))
paste(colSums(is.na(control.df.qua))*100/dim(control.df.qua)[1],'% missing values for control group', colnames(control.df.qua))


missing_name <- colnames(control.df.qua)[which(as.vector(colSums(is.na(control.df.qua)) >0))]
missing_name <- as.vector(missing_name)
control.df.qua <- as.data.frame(control.df.qua)
for(i in missing_name){
  indices <- which(is.nan(control.df.qua[,i]))
  for(j in indices){
    if(j == 0){
      control.df.qua[j,i] <- ifelse(is.na(control.df.qua[1,i]),control.df.qua[2,i],control.df.qua[1,i])
    }else{
      if(j == length(control.df.qua$season)){
        control.df.qua[j,i] <- ifelse(is.na(control.df.qua[j-1,i]),control.df.qua[j-2,i],control.df.qua[j-1,i])
      }else{
        control.df.qua[j,i] <- mean(c(control.df.qua[j-1,i],control.df.qua[j+1,i]))
      }}
  }
}
paste(colSums(is.na(control.df.qua))*100/dim(control.df.qua)[1],'% missing values for control group', colnames(control.df.qua))

#######################################################################
################ Applying BSTS (Brodersen et al. 2015) ################
#######################################################################

bsts.df.qua <- merge(treated.df.qua, control.df.qua)
treated.zoo.qua <- read.zoo(treated.df.qua, index_column=1)
bsts.zoo.qua <-read.zoo(bsts.df.qua, index_column=1)
bsts.white.qua <- bsts.zoo.qua[,colnames(bsts.zoo.qua) %in% c('ghq','1','2','4')]
pre.period.qua <- c(min(index(bsts.zoo.qua)),-1)
post.period.qua <- c(0,max(index(bsts.zoo.qua)))
post.y.qua <- post.period.response.qua <- as.vector(bsts.zoo.qua$ghq[index(bsts.zoo.qua) >=0])
pre.time.qua <- pre.period.qua[2] - pre.period.qua[1] + 1
post.time.qua <- post.period.qua[2] - post.period.qua[1] + 1
#bsts.zoo.norm <- log(bsts.zoo.qua + 1)
# --  Causalimpact package --
hold.qua <- bsts.zoo.qua
hold.qua[index(bsts.zoo.qua) >=0,1] <- NA
ci <- CausalImpact(bsts.zoo.qua,pre.period = pre.period.qua,post.period = post.period.qua,
                   model.args = list(nseasons=1, prior.level.sd = 0.01))
plot(ci) 

# -- use our adapted code --
ce<- run_flex_bsts(data=bsts.zoo.qua,pre.period=pre.period.qua,
                   post.period=post.period.qua,
                   alpha=0.05,model.args=list(nseasons=1,
                                              prior.level.sd = 0.01,
                                              standardize.data=TRUE),
                   trend = 'AR1',modelsize=3,cycles = c(1,2,4,8),niter=10000)
plot(ce)

RMSE(ce$series$point.pred[1:65],bsts.zoo.qua$ghq[1:65]) #0.956
RMSE(ce$series$point.pred[2:65],bsts.zoo.qua$ghq[2:65]) #0.958
RMSE(ce$series$point.pred[3:65],bsts.zoo.qua$ghq[3:65]) #0.966
mean(abs(ce$series$point.pred[1:65] - bsts.zoo.qua$ghq[1:65])) # 0.709
mean(abs(ce$series$point.pred[1:65]-bsts.zoo.qua$ghq[1:65])/(abs(ce$series$point.pred[1:65]) + abs(bsts.zoo.qua$ghq[1:65]))) # 0.03109709

#diff <- quantile(rowMeans(matrix(rep(emp.mean.y[23:36],each= 9931),ncol=14) -(exp(ce$model$posterior.samples[,23:36])-1) ),c(0.025,0.975))
#mean(emp.mean.y[23:36] - colMeans(exp(ce$model$posterior.samples[,23:36])-1))
## not make much difference normalise or not
###################################################################################################
################ Applying HS-BSTS (standard Bsts with regularised horseshoe prior) ################
###################################################################################################

ce.rhs<- run_flex_bsts(data=bsts.white.qua,pre.period=pre.period.qua,
                   post.period=post.period.qua,
                   alpha=0.05,model.args=list(nseasons=1,
                                              prior.level.sd = 0.01,
                                              standardize.data=TRUE),
                   trend = 'AR1',modelsize=3,cycles = c(1),niter=10000)
plot(ce.rhs)


################################################
################ Applying HBSTS ################
################################################




treated.multi <- treated.whole.df %>% dplyr::select(season,ghq,strata,psu,wave) %>% dplyr::arrange(season)
all.multi <- merge(treated.multi, control.df.qua, by='season')
colnames(all.multi) <- c('season','carib','British','Irish','Other White',
                         'mixed asian',
                         'other mixed','Indian',
                         'pakistani','Bangaldeshi','Chinese',
                         'other asian','African','other black','arab',
                         'Other')

pre.ys.qua <- as.vector(treated.multi$ghq[treated.multi$season <0])
post.ys.qua <- as.vector(treated.multi$ghq[treated.multi$season >=0])
#pre.ys.qua <- log(as.vector(treated.multi$ghq[treated.multi$season <=0]) + 1)
#post.ys.qua <- log(as.vector(treated.multi$ghq[treated.multi$season >0]) + 1)
pre.obs.qua <- sum(treated.multi$season < 0)
post.obs.qua <- sum(treated.multi$season >= 0)
pre.strata.qua <- as.vector(treated.multi$strata[treated.multi$season <0])
post.strata.qua <- as.vector(treated.multi$strata[treated.multi$season >=0])
obs.count <- as.vector(table(treated.multi$season))
index_pre <- cumsum(c(1,obs.count[1:pre.time.qua]))
index_post <- cumsum(c(1,obs.count[(pre.time.qua+1):(pre.time.qua+post.time.qua)]))

# ---- control group: White only -----
control.whole.df <- eth.df %>%
  dplyr::filter(eth1 %in% white_groups)%>%
  dplyr::mutate(season = (lubridate::interval(start = policyDate, end = interviewDate) %/% months(1)))
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
#x.obs <- log(x.obs+1)
K <- length(x.obs)
emp.mean.y <- as.vector(bsts.zoo.qua[,1])

p.mat <- as.matrix(x.count)
p.mat[p.mat!=0] <- 1


hmodmw <- run_hbsts(pre.time=pre.time.qua,post.time=post.time.qua,
                       pre.obs=pre.obs.qua,
                       pre.ys=pre.ys.qua,dimx=3,
                       post.obs =post.obs.qua,index_pre=index_pre,
                       K=K,x.obs=x.obs, x.cum.count=x.cum.count,
                       presence_matrix = p.mat,
                       index_post=index_post,
                       niter=10000,nchains=1,time_index=pre.period.qua[1]:post.period.qua[2],
                       alpha=0.05,cycles=c(1,2,4,8),emp.mean.y=emp.mean.y,
                       trend = 'AR1',burn=3000,modelsize=3,adapt=0.9,ind=FALSE,hdi=FALSE)

RMSE(hmodmw$series$point.pred[1:65],emp.mean.y[1:65]) #  1.004975
RMSE(hmodmw$series$point.pred[2:65],emp.mean.y[2:65]) # 0.7318236
RMSE(hmodmw$series$point.pred[3:65],emp.mean.y[3:65]) #  0.6809745
mean(abs(hmodmw$series$point.pred[1:65] - emp.mean.y[1:65])) # 0.6092669
mean(abs(hmodmw$series$point.pred[1:65]-as.vector(emp.mean.y[1:65]))/(abs(hmodmw$series$point.pred[1:65]) + abs(as.vector(emp.mean.y[1:65])))) # 0.02889844

pre.strata.qua <- as.vector(treated.multi$strata[treated.multi$season <0])
post.strata.qua <- as.vector(treated.multi$strata[treated.multi$season >=0])
x.strata <- control.whole.df %>% dplyr::select(strata,eth1,season)
x.strata <- x.strata %>% dplyr::arrange(season,eth1)
x.strata <- as.vector(x.strata[,c(1)])
NS <- length(unique(eth.df$strata))

pre.psu.qua <- as.vector(treated.multi$psu[treated.multi$season <0])
post.psu.qua <- as.vector(treated.multi$psu[treated.multi$season >=0])
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
                        niter=10000,nchains=1,time_index=pre.period.qua[1]:post.period.qua[2],
                        alpha=0.05,cycles=c(1,2,4,8),emp.mean.y=emp.mean.y,
                        trend = 'AR1',burn=1000,modelsize=2,adapt=0.9,presence_matrix = p.mat,
                        x_strata=x.strata,NS=NS,pre_strata = pre.strata.qua,
                        x_psu=x.psu,NC=NC,pre_psu = pre.psu.qua,ind=FALSE,hdi=FALSE)
hmodmw$summary
plot(hmodmw)

# ---- control group: all other 14 groups -----

control.whole.df <- eth.df %>%
  dplyr::filter(eth1 %in% donor_pool)%>%
  dplyr::mutate(season = (lubridate::interval(start = policyDate, end = interviewDate) %/% months(1)))
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
p.mat <- as.matrix(x.count)
p.mat[p.mat!=0] <- 1

pre.strata.qua <- as.vector(treated.multi$strata[treated.multi$season <0])
post.strata.qua <- as.vector(treated.multi$strata[treated.multi$season >=0])
x.strata <- control.whole.df %>% dplyr::select(strata,eth1,season)
x.strata <- x.strata %>% dplyr::arrange(season,eth1)
x.strata <- as.vector(x.strata[,c(1)])
NS <- length(unique(eth.df$strata))

pre.psu.qua <- as.vector(treated.multi$psu[treated.multi$season <0])
post.psu.qua <- as.vector(treated.multi$psu[treated.multi$season >=0])
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
                     niter=8000,nchains=1,time_index=pre.period.qua[1]:post.period.qua[2],
                     alpha=0.05,cycles=c(1),emp.mean.y=emp.mean.y,
                     trend = 'AR1',burn=3000,modelsize=3,adapt=0.9,
                     x_strata=x.strata,NS=NS,pre_strata = pre.strata.qua,presence_matrix = p.mat,
                     x_psu=x.psu,NC=NC,pre_psu = pre.psu.qua,ind=FALSE,hdi=FALSE)


hmod$summary
plot(hmod)

###########################################
################### plots #################
###########################################

bstsmod <- ce
hmod <- hmodmw

## --- (1) comparative counterfactual plot (bsts vs hbsts vs observed), shown in main text
allplot <- data.frame(pred.mean = c(colMeans(((bstsmod$model$posterior.samples))),colMeans(((hmod$model$posterior.samples))),bsts.zoo.mon$ghq),
                      low=c(apply(((bstsmod$model$posterior.samples)),2, quantile,0.025),apply(((hmod$model$posterior.samples)),2, quantile,0.025),bsts.zoo.mon$ghq),
                      upper=c(apply(((bstsmod$model$posterior.samples)),2,quantile,0.975),apply(((hmod$model$posterior.samples)),2, quantile,0.975),bsts.zoo.mon$ghq),
                      #time=c(-64:42,-64:42),
                      time=rep(seq.Date(from=as.Date('2009-01-10'),by='1 month',length.out=107),3),
                      type=rep(c('Prediction (BSTS)', 'Prediction (HBSTS)','Observation'),each=107))


ggplot(aes(x=time,y=pred.mean,color=type,group=type,fill=type),data=allplot) + 
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(aes(ymin = low, ymax = upper), alpha = 0.25, colour = NA) +
  geom_vline(xintercept = as.Date('2014-05-10'),linetype = "dashed")+
  # ggplot2::scale_y_continuous(limits = c(10, 14), breaks = seq(10, 14, by = 1)) +
  ggplot2::labs(x = 'Time (year-month)', 
                y = 'GHQ-12 score') +
  scale_x_date(date_labels = "%Y-%m", 
               breaks = seq.Date(from=as.Date('2009-01-10'),by='8 month',length.out=107))+
  ggplot2::scale_colour_manual(values = c('black', 'orange2','blue'))+#,"#009966")) +
  ggplot2::scale_fill_manual(values = c('black', 'orange2','blue'))+
  ggplot2::scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5))+
  my.theme(legend.title = element_blank(),
           legend.position = 'bottom',
           text = element_text(size = 16),
           axis.text.x = element_text(angle = 60, hjust = 1,size=14)); 


## --- (2) pointwise plot for bsts
bsts.quarter <- data.frame(effect = bstsmod$temporaleffect$temp_effect[[2]]$tempeffect,
                           lower = bstsmod$temporaleffect$temp_effect[[2]]$tempeffect.lower,
                           upper = bstsmod$temporaleffect$temp_effect[[2]]$tempeffect.upper,
                           p = bstsmod$temporaleffect$p_vals[[2]],
                           time = 1:22)
bsts.quarter <- bsts.quarter %>%
  mutate(sig = ifelse(p < 0.05, "Significant", "Not Significant"))
ggplot(bsts.quarter, aes(x = time, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y = 0
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +       # Error bars
  geom_point(aes(color = sig), size = 3) +                             # Points colored by significance
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(
    x = "Time",
    y = "Effect",
    color = "P < 0.05"
  )+my.theme(legend.title = element_blank(),
             legend.position = 'bottom',
             text = element_text(size = 14));


bsts.halfyear <- data.frame(effect = bstsmod$temporaleffect$temp_effect[[3]]$tempeffect,
                            lower = bstsmod$temporaleffect$temp_effect[[3]]$tempeffect.lower,
                            upper = bstsmod$temporaleffect$temp_effect[[3]]$tempeffect.upper,
                            p = bstsmod$temporaleffect$p_vals[[3]],
                            time = 1:11)
bsts.halfyear <- bsts.halfyear %>%
  mutate(sig = ifelse(p < 0.05, "Significant", "Not Significant"))
ggplot(bsts.halfyear, aes(x = time, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y = 0
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +       # Error bars
  geom_point(aes(color = sig), size = 3) +                             # Points colored by significance
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(
    x = "Time",
    y = "Effect",
    color = "P < 0.05"
  )+my.theme(legend.title = element_blank(),
             legend.position = 'bottom',
             text = element_text(size = 14));




bsts.year <- data.frame(effect = bstsmod$temporaleffect$temp_effect[[4]]$tempeffect,
                        lower = bstsmod$temporaleffect$temp_effect[[4]]$tempeffect.lower,
                        upper = bstsmod$temporaleffect$temp_effect[[4]]$tempeffect.upper,
                        p = bstsmod$temporaleffect$p_vals[[4]],
                        time = 1:6)
bsts.year <- bsts.year %>%
  mutate(sig = ifelse(p < 0.05, "Significant", "Not Significant"))
ggplot(bsts.year, aes(x = time, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y = 0
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05) +       # Error bars
  geom_point(aes(color = sig), size = 3) +                             # Points colored by significance
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(
    x = "Time",
    y = "Effect",
    color = "P < 0.05"
  )+my.theme(legend.title = element_blank(),
             legend.position = 'bottom',
             text = element_text(size = 14));



## --- (3) pointwise plot for hbsts, appear in the supplement 



hbstsplot <- data.frame(pred.mean = c(colMeans(((hmod$model$posterior.samples))),bsts.zoo.mon$ghq),
                        low=c(apply(((hmod$model$posterior.samples)),2, quantile,0.025),bsts.zoo.mon$ghq),
                        upper=c(apply(((hmod$model$posterior.samples)),2,quantile,0.975),bsts.zoo.mon$ghq),
                        # time=c(-64:42,-64:42)
                        time=rep(seq.Date(from=as.Date('2009-01-10'),by='1 month',length.out=107),2),
                        type=rep(c('Prediction','Observation'),each=107))

ggplot(aes(x=time,y=pred.mean,color=type,group=type,fill=type),data=hbstsplot) + 
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(aes(ymin = low, ymax = upper), alpha = 0.25, colour = NA) +
  geom_vline(xintercept = as.Date('2014-05-10'),linetype = "dashed")+
  # ggplot2::scale_y_continuous(limits = c(10, 14), breaks = seq(10, 14, by = 1)) +
  ggplot2::labs(x = 'Time (year-month)', 
                y = 'GHQ-12 score') +
  scale_x_date(date_labels = "%Y-%m", 
               breaks = seq.Date(from=as.Date('2009-01-10'),by='8 month',length.out=107))+
  ggplot2::scale_colour_manual(values = c('black', 'orange2'))+#,"#009966")) +
  ggplot2::scale_fill_manual(values = c('black', 'orange2'))+
  ggplot2::scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5))+
  my.theme(legend.title = element_blank(),
           legend.position = 'bottom',
           text = element_text(size = 14),
           axis.text.x = element_text(angle = 60, hjust = 1,size=12)); 


hbsts.mon <- data.frame(effect = hmod$series$point.effect[66:107],
                        lower = hmod$series$point.effect.lower[66:107],
                        upper = hmod$series$point.effect.upper[66:107],
                        sig = as.numeric(hmod$series$point.effect.lower[66:107] >0),
                        time = seq.Date(from=as.Date('2014-07-13'),by='1 month',length.out = 42))
hbsts.mon <- hbsts.mon %>%
  mutate(sig = ifelse(sig == 1, "Significant", "Not Significant"))
hbsts.mon$sig <- as.factor(hbsts.mon$sig)
ggplot(hbsts.mon, aes(x = time, y = effect)) +
  geom_line()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y = 0
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.25,fill='gray') +       # Error bars
  geom_point(aes(color = sig), size = 3) +    
  scale_x_date(date_labels = "%Y-%m", 
               breaks = seq.Date(from=as.Date('2014-09-13'),by='3 month',length.out=14))+# Points colored by significance
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(
    x = "Time (year-month)",
    y = "Effect",
    color = "P < 0.05"
  )+my.theme(legend.title = element_blank(),
             legend.position = 'bottom',
             text = element_text(size = 14),
             axis.text.x = element_text(angle = 60, hjust = 1,size=12));




hbsts.mon.cum <- data.frame(effect = hmod$series$cum.effect[66:107],
                            lower = hmod$series$cum.effect.lower[66:107],
                            upper = hmod$series$cum.effect.upper[66:107],
                            sig = as.numeric(hmod$series$cum.effect.lower[66:107] >0),
                            time = seq.Date(from=as.Date('2014-07-13'),by='1 month',length.out = 42))
hbsts.mon.cum <- hbsts.mon.cum %>%
  mutate(sig = ifelse(sig == 1, "Significant", "Not Significant"))
hbsts.mon.cum$sig <- as.factor(hbsts.mon.cum$sig)
ggplot(hbsts.mon.cum, aes(x = time, y = effect)) +
  geom_line()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y = 0
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.25,fill='gray') +       # Error bars
  geom_point(aes(color = sig), size = 3) +    
  scale_x_date(date_labels = "%Y-%m", 
               breaks = seq.Date(from=as.Date('2014-09-13'),by='3 month',length.out=14))+# Points colored by significance
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(
    x = "Time (year-month)",
    y = "Effect",
    color = "P < 0.05"
  )+my.theme(legend.title = element_blank(),
             legend.position = 'bottom',
             text = element_text(size = 14),
             axis.text.x = element_text(angle = 60, hjust = 1,size=12));


samples <- hmod$model$posterior.samples[,66:107]
obsy <- matrix(rep(emp.mean.y[66:107],dim(samples)[1]),byrow = T,nrow=dim(samples)[1])
samp.eff <- obsy - samples
qua.effect <- sapply(seq(1, 40, by = 3), function(i){
  mean(samp.eff[,i:(i+2)])
})
qua.low <- sapply(seq(1, 40, by = 3), function(i){
  quantile(rowMeans(samp.eff[,i:(i+2)]),0.025)
})
qua.up <- sapply(seq(1, 40, by = 3), function(i){
  quantile(rowMeans(samp.eff[,i:(i+2)]),0.975)
})


dates <- seq.Date(from = as.Date("2014-09-12"), by = "3 months", length.out = 14)
labels <- format(dates, "%Y-%m")
hbsts.quarter <- data.frame(effect = qua.effect,lower=qua.low,upper=qua.up,
                            sig = as.numeric(qua.low >0), 
                            time=seq.Date(from = as.Date("2014-09-12"), by = "3 months", length.out = 14))
hbsts.quarter <- hbsts.quarter %>%
  mutate(sig = ifelse(sig == 1, "Significant", "Not Significant"))
hbsts.quarter$sig <- as.factor(hbsts.quarter$sig)



ggplot(hbsts.quarter, aes(x = time, y = effect)) +
  geom_line()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y = 0
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.25,fill='gray') +       # Error bars
  geom_point(aes(color = sig), size = 3) +    
  scale_x_date(date_labels = "%Y-%m", breaks = hbsts.quarter$time)+# Points colored by significance
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(
    x = "Time (year-month)",
    y = "Effect",
    color = "P < 0.05"
  )+my.theme(legend.title = element_blank(),
             legend.position = 'bottom',
             text = element_text(size = 14),
             axis.text.x = element_text(angle = 60, hjust = 1,size=12));


half.effect <- sapply(seq(1, 37, by = 6), function(i){
  mean(samp.eff[,i:(i+5)])
})
half.low <- sapply(seq(1, 37, by = 6), function(i){
  quantile(rowMeans(samp.eff[,i:(i+5)]),0.025)
})
half.up <- sapply(seq(1, 37, by = 6), function(i){
  quantile(rowMeans(samp.eff[,i:(i+5)]),0.975)
})

dates <- seq.Date(from = as.Date("2014-12-13"), by = "6 months", length.out = 7)
labels <- format(dates, "%Y-%m")

hbsts.half <- data.frame(effect = half.effect,lower=half.low,upper=half.up,
                         sig = as.numeric(half.low >0), time=seq.Date(from = as.Date("2014-12-13"), by = "6 months", length.out = 7))
hbsts.half<- hbsts.half %>%
  mutate(sig = ifelse(sig == 1, "Significant", "Not Significant"))
hbsts.half$sig <- as.factor(hbsts.half$sig)


ggplot(hbsts.half, aes(x = time, y = effect)) +
  geom_line()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y = 0
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, fill = "gray") +       # Error bars
  geom_point(aes(color = sig), size = 3) +   
  scale_x_date(date_labels = "%Y-%m", breaks = hbsts.half$time)+# Points colored by significance
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(
    x = "Time (year-month)",
    y = "Effect",
    color = "P < 0.05"
  )+my.theme(legend.title = element_blank(),
             legend.position = 'bottom',
             text = element_text(size = 14),
             axis.text.x = element_text(angle = 60, hjust = 1,size=12));


year.effect <- sapply(seq(1, 37, by = 12), function(i){
  mean(samp.eff[,i:min((i+11),42)])
})

year.low <- sapply(seq(1, 37, by = 12), function(i){
  quantile(rowMeans(samp.eff[,i:min((i+11),42)]),0.025)
})
year.up <- sapply(seq(1, 37, by = 12), function(i){
  quantile(rowMeans(samp.eff[,i:min((i+11),42)]),0.975)
})

x <- as.Date(c("2015-06-13", "2016-6-13", "2017-6-13","2017-11-27"))

formatted <- format(x, "%Y-%m")
hbsts.year <- data.frame(effect = year.effect,lower=year.low,upper=year.up,
                         sig = as.numeric(year.low >0), time=as.Date(c("2015-06-13", "2016-6-13", "2017-6-13","2017-11-27")))
hbsts.year<- hbsts.year %>%
  mutate(sig = ifelse(sig == 1, "Significant", "Not Significant"))
hbsts.year$sig <- as.factor(hbsts.year$sig)


ggplot(hbsts.year, aes(x = time, y = effect)) +
  geom_line()+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, fill = "gray") +   
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y = 0
  # Error bars
  geom_point(aes(color = sig), size = 3) + 
  scale_x_date(date_labels = "%Y-%m", breaks = hbsts.year$time)+
  # Points colored by significance
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(
    x = "Time (year-month)",
    y = "Effect",
    color = "P < 0.05"
  )+my.theme(legend.title = element_blank(),
             legend.position = 'bottom',
             text = element_text(size = 14),axis.text.x = element_text(angle = 60, hjust = 1,size=12));



## --- (4) inclusion probabilities (not in paper)

fit <- rstan::extract(hmod$model$bsts.model)
mean(fit$sigma_mu)
mean(fit$sigma_obs)
colMeans(fit$sigma_obsx)
mean(fit$sigma_global)
comp <- data.frame(
  group = c("White British", "White Irish", "Other White"),
  probability = colMeans(fit$rho!=0)
)
comp$group <- factor(comp$group, levels = comp$group[order(comp$probability)])

ggplot(comp, aes(x = group, y = probability)) +
  geom_col(fill = "steelblue") +
  coord_flip() +  # flips x and y: horizontal bars
  labs(x = NULL, y = "Posterior inclusion probability") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.line = element_line()
  )


fit <- rstan::extract(hmod14$model$bsts.model)

comp <- data.frame(
  group = c('White British','White Irish','Other White',
            'Mixed Asian',
            'Other Mixed','Indian',
            'Pakistani','Bangaldeshi','Chinese',
            'other Asian','African','Other Black','Arab',
            'Other'),
  probability = colMeans(fit$rho!=0)
)
comp$group <- factor(comp$group, levels = comp$group[order(comp$probability)])

ggplot(comp, aes(x = group, y = probability)) +
  geom_col(fill = "steelblue") +
  coord_flip() +  # flips x and y: horizontal bars
  labs(x = NULL, y = "Posterior inclusion probability") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.line = element_line()
  )



