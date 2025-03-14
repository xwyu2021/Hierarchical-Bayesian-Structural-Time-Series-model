#########################################################
## simulated data are saved, here we simply load them ###
#########################################################
library(coda)
library(bayesplot)
library(RColorBrewer)
library(viridis)

###############################
### convergence diagnostics ###
###############################


## Geweke Diagnostic Test
geweke.diag <- function(chain, first = 0.1, last = 0.5, width = 0.1) {
  # Validate inputs
  if (!is.numeric(chain)) stop("Chain must be numeric")
  if (first <= 0 || first >= 1) stop("'first' must be between 0 and 1")
  if (last <= 0 || last >= 1) stop("'last' must be between 0 and 1")
  if (width <= 0 || width >= 1) stop("'width' must be between 0 and 1")
  
  n <- length(chain)
  
  # Calculate indices for first and last portions
  n1 <- floor(first * n)
  n2 <- floor((1 - last) * n)
  
  # Extract first and last portions
  first_portion <- chain[1:n1]
  last_portion <- chain[n2:n]
  
  # Calculate spectral densities using spectrum0.ar()
  spec1 <- spectrum0.ar(first_portion)
  spec2 <- spectrum0.ar(last_portion)
  
  # Calculate means
  mean1 <- mean(first_portion)
  mean2 <- mean(last_portion)
  
  # Calculate standard error of the difference in means
  se <- sqrt(spec1$spec/length(first_portion) + spec2$spec/length(last_portion))
  
  # Calculate z-score
  z <- (mean1 - mean2) / se
  
  # Calculate two-sided p-value explicitly
  # Two-sided p-value = 2 * P(Z > |z|) = 2 * (1 - P(Z ≤ |z|))
  two_sided_p <- 2 * (1 - pnorm(abs(z)))
  
  # Create result object
  result <- list(
    z = z,
    p.value = two_sided_p,
    first.mean = mean1,
    last.mean = mean2,
    first.proportion = first,
    last.proportion = last,
    spec1 = spec1$spec,
    spec2 = spec2$spec
  )
  
  class(result) <- "geweke.diag"
  return(result)
}

# Print method for geweke.diag objects
print.geweke.diag <- function(x, digits = 3, ...) {
  cat("\nGeweke Convergence Diagnostic\n\n")
  cat("Z-score:", round(x$z, digits), "\n")
  cat("Two-sided P-value:", round(x$p.value, digits), "\n")
  cat("First window mean:", round(x$first.mean, digits), "\n")
  cat("Last window mean:", round(x$last.mean, digits), "\n")
}

# load the simulated dataset (first experiment)
compare.data <- readRDS('simulated data/simulation_compare_hbststrend.RData')
hsim.ar <- compare.data$hsim.ar
hsim.ll <- compare.data$hsim.ll
hsim.llt <- compare.data$hsim.llt
hsim.gllt <- compare.data$hsim.gllt

modar <- rstan::extract(hsim.ar$model$bsts.model)
geweke.diag(as.matrix(modar$sigma_mu))
geweke.diag(as.matrix(modar$sigma_obs))
geweke.diag(as.matrix(modar$sigma_y))
geweke.diag(as.matrix(modar$sigma_xx))
geweke.diag(modar$kappa[,1])
geweke.diag(as.matrix(modar$beta0))

modll <- rstan::extract(hsim.ll$model$bsts.model)
geweke.diag(as.matrix(modll$sigma_mu))
geweke.diag(as.matrix(modll$sigma_obs))
geweke.diag(as.matrix(modll$sigma_y))
geweke.diag(as.matrix(modll$sigma_xx))
geweke.diag(modll$kappa[,1])
geweke.diag(as.matrix(modll$beta0))



modllt <- rstan::extract(hsim.llt$model$bsts.model)
geweke.diag(as.matrix(modllt$sigma_mu))
geweke.diag(as.matrix(modllt$sigma_obs))
geweke.diag(as.matrix(modllt$sigma_y))
geweke.diag(as.matrix(modllt$sigma_xx))
geweke.diag(modllt$kappa[,1])
geweke.diag(as.matrix(modllt$beta0))



modgllt <- rstan::extract(hsim.gllt$model$bsts.model)
geweke.diag(as.matrix(modgllt$sigma_mu))
geweke.diag(as.matrix(modgllt$sigma_obs))
geweke.diag(as.matrix(modgllt$sigma_y))
geweke.diag(as.matrix(modgllt$sigma_xx))
geweke.diag(modgllt$kappa[,1])
geweke.diag(as.matrix(modgllt$beta0))



y_pred <- as.matrix(modar$alpha)
ppc_dens_overlay(simemp.mean.y[1:35],y_pred) # simemp.mean.y is the real observed outcome after intervention (see simulation.R)

stan_dens(hsim.ar$model$bsts.model,c('sigma_mu','sigma_y','sigma_obs','sigma_xx'),nrow=1)
stan_trace(hsim.ar$model$bsts.model,c('sigma_mu','sigma_y','sigma_obs','sigma_xx'),nrow=1) + ylab('')
check_hmc_diagnostics(hsim.ar$model$bsts.model)
check_divergences(hsim.ar$model$bsts.model)
check_treedepth(hsim.ar$model$bsts.model)
check_energy(hsim.ar$model$bsts.model)
get_bfmi(hsim.ar$model$bsts.model)
get_low_bfmi_chains(hsim.ar$model$bsts.model)

Rhat(as.matrix(modar$sigma_mu,ncol=2))
ess_bulk(as.matrix(modar$sigma_mu,ncol=2))
ess_tail(as.matrix(modar$sigma_mu,ncol=2))




###############################################################
### robustness and effectiveness of detecting an effect #######
###############################################################

## e.g. first experiement
nomiss <- readRDS('simulated data/sim100_50t_5positive/sim100obs_positive.RData')
mnar <-readRDS('simulated data/sim100_50t_5positive/sim100MCARobs_positive.RData')
mcar <- readRDS('simulated data/sim100_50t_5positive/sim100MCARobs_positive.RData')

### --- plot of probability of detecting an effect --- ###

k <- cbind(colSums(1-nomiss$spec),apply(nomiss$sen,3,colSums))
n_k <- 200-k
point.hbsts <- point.hsbsts <- point.ssbsts <- matrix(0,nrow = 2,ncol=8)
for(i in 1:8){
  point.hbsts[,i] <- quantile(rbeta(10000,1+k[1,i],1+n_k[1,i]),c(0.025,0.975))
  point.hsbsts[,i] <- quantile(rbeta(10000,1+k[2,i],1+n_k[2,i]),c(0.025,0.975))
  point.ssbsts[,i] <- quantile(rbeta(10000,1+k[3,i],1+n_k[3,i]),c(0.025,0.975))
}


df.nomiss <- data.frame(mean = as.vector(t(k))/200,
                        lower = c(point.hbsts[1,],point.hsbsts[1,],point.ssbsts[1,]),
                        upper = c(point.hbsts[2,],point.hsbsts[2,],point.ssbsts[2,]),
                        model = rep(c('HBSTS','HS-BSTS','SS-BSTS'),each=8),
                        increment = rep(c(0,0.01,0.1,0.3,0.5,1,2,3),3))
df.nomiss <- df.nomiss[df.nomiss$increment != 0.01,]

ggplot(df.nomiss, aes(x=increment, y=mean,group=model,color=model)) + 
  geom_line(size=0.6) +
  geom_point(size=3.5)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.3,
                position=position_dodge(0.05),size=0.6) + 
  theme_bw() + theme(axis.title = element_text(size = 20),
                     axis.text= element_text(size=20),
                     legend.position='none',
                     axis.text.x = element_text(angle = 60, hjust = 1,size=18),
                     panel.grid.minor = element_blank())+
  scale_color_manual(values=c('HBSTS'='orange','HS-BSTS'='skyblue2','SS-BSTS'='grey2')) +
  labs(x='effect size',y='proportion of intervals excluding 0') +
  scale_x_continuous(breaks = c(0,0.1,0.3,0.5,1,2,3),
                     labels = c(0,0.1,0.3,0.5,1,2,3)) + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),
                     labels = c(0,0.2,0.4,0.6,0.8,1),
                     limits = c(0,1))


k <- cbind(colSums(1-mcar$spec),apply(mcar$sen,3,colSums))
n_k <- 200-k
point.hbsts <- point.hsbsts <- point.ssbsts <- matrix(0,nrow = 2,ncol=8)
for(i in 1:8){
  point.hbsts[,i] <- quantile(rbeta(10000,1+k[1,i],1+n_k[1,i]),c(0.025,0.975))
  point.hsbsts[,i] <- quantile(rbeta(10000,1+k[2,i],1+n_k[2,i]),c(0.025,0.975))
  point.ssbsts[,i] <- quantile(rbeta(10000,1+k[3,i],1+n_k[3,i]),c(0.025,0.975))
}
df.mcar <- data.frame(mean = as.vector(t(k))/200,
                      lower = c(point.hbsts[1,],point.hsbsts[1,],point.ssbsts[1,]),
                      upper = c(point.hbsts[2,],point.hsbsts[2,],point.ssbsts[2,]),
                      model = rep(c('HBSTS','HS-BSTS','SS-BSTS'),each=8),
                      increment = rep(c(0,0.01,0.1,0.3,0.5,1,2,3),3))
df.mcar <- df.mcar[df.mcar$increment != 0.01,]
ggplot(df.mcar, aes(x=increment, y=mean,group=model,color=model)) + 
  geom_line(size=0.6) +
  geom_point(size=3.5)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(0.05),size=0.6) + 
  theme_bw() + theme(axis.title = element_text(size = 18),
                     axis.text= element_text(size=18),
                     legend.position='none',
                     axis.text.x = element_text(angle = 60, hjust = 1,size=18),
                     panel.grid.minor = element_blank())+
  scale_color_manual(values=c('HBSTS'='orange','HS-BSTS'='skyblue2','SS-BSTS'='grey2')) +
  labs(x='effect size',y='proportion of intervals excluding 0') +
  scale_x_continuous(breaks = c(0,0.1,0.3,0.5,1,2,3),
                     labels = c(0,0.1,0.3,0.5,1,2,3))+ 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),
                     labels = c(0,0.2,0.4,0.6,0.8,1),
                     limits = c(0,1))


k <- cbind(colSums(1-mnar$spec),apply(mnar$sen,3,colSums))
n_k <- 200-k
point.hbsts <- point.hsbsts <- point.ssbsts <- matrix(0,nrow = 2,ncol=8)
for(i in 1:8){
  point.hbsts[,i] <- quantile(rbeta(10000,1+k[1,i],1+n_k[1,i]),c(0.025,0.975))
  point.hsbsts[,i] <- quantile(rbeta(10000,1+k[2,i],1+n_k[2,i]),c(0.025,0.975))
  point.ssbsts[,i] <- quantile(rbeta(10000,1+k[3,i],1+n_k[3,i]),c(0.025,0.975))
}

df.mnar <- data.frame(mean = as.vector(t(k))/200,
                      lower = c(point.hbsts[1,],point.hsbsts[1,],point.ssbsts[1,]),
                      upper = c(point.hbsts[2,],point.hsbsts[2,],point.ssbsts[2,]),
                      model = rep(c('HBSTS','HS-BSTS','SS-BSTS'),each=8),
                      increment = rep(c(0,0.01,0.1,0.3,0.5,1,2,3),3))

df.mnar <- df.mnar[df.mnar$increment !=0.01,]
ggplot(df.mnar, aes(x=increment, y=mean,group=model,color=model)) + 
  geom_line(size=0.6) +
  geom_point(size=3.5)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(0.05),size=0.6) + 
  theme_bw() + theme(axis.title = element_text(size = 18),
                     axis.text= element_text(size=18),
                     legend.position='none',
                     axis.text.x = element_text(angle = 60, hjust = 1,size=18),
                     panel.grid.minor = element_blank())+
  scale_color_manual(values=c('HBSTS'='orange','HS-BSTS'='skyblue2','SS-BSTS'='grey2')) +
  labs(x='effect size',y='proportion of intervals excluding 0') +
  scale_x_continuous(breaks = c(0,0.1,0.3,0.5,1,2,3),
                     labels = c(0,0.1,0.3,0.5,1,2,3)) + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),
                     labels = c(0,0.2,0.4,0.6,0.8,1),
                     limits = c(0,1))



### --- plot of APEE --- ###

apee1 <- abs(mcar$avgeffectall[,1,1] -  mcar$realeffect[,1,1])/mcar$realeffect[,1,1]
apee2 <- abs(mcar$avgeffectall[,2,1] -  mcar$realeffect[,1,1])/mcar$realeffect[,1,1]
apee3 <- abs(mcar$avgeffectall[,3,1] -  mcar$realeffect[,1,1])/mcar$realeffect[,1,1]
apee.whole <- data.frame(mean = c(mean(apee1),as.vector(apply(mcar$apee,2,colMeans))[2:7],
                                  mean(apee2),as.vector(apply(mcar$apee,2,colMeans))[9:14],
                                  mean(apee3),as.vector(apply(mcar$apee,2,colMeans))[16:21]),
                         lower=c(quantile(apee1,0.025),as.vector(apply(mcar$apee,2,apply,2,quantile,0.025))[2:7],
                                 quantile(apee2,0.025),as.vector(apply(mcar$apee,2,apply,2,quantile,0.025))[9:14],
                                 quantile(apee3,0.025),as.vector(apply(mcar$apee,2,apply,2,quantile,0.025))[16:21]),
                         upper = c(quantile(apee1,0.975),as.vector(apply(mcar$apee,2,apply,2,quantile,0.975))[2:7],
                                   quantile(apee2,0.975),as.vector(apply(mcar$apee,2,apply,2,quantile,0.975))[9:14],
                                   quantile(apee3,0.975),as.vector(apply(mcar$apee,2,apply,2,quantile,0.975))[16:21]),
                         model = rep(c('HBSTS','HS-BSTS','SS-BSTS'), each =7),
                         increment = rep(c(0,0.1,0.3,0.5,1,2,3),3))

ggplot(apee.whole,aes(x=increment,y=mean,group=model,color=model)) + 
  geom_point(size=2,aes(color=model)) + 
  geom_errorbar(aes(ymin=lower,ymax=upper,color=model),width=.1,
                position=position_dodge(0.05),size=0.5)+
  theme_bw() + theme(axis.title = element_text(size = 20),
                     axis.text= element_text(size=16),
                     legend.position=c(0.8,0.8),
                     legend.text = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     axis.text.x = element_text(angle = 60, hjust = 1),
                     panel.grid.minor = element_blank(),
                     strip.text = element_text(size = 12))+
  #scale_color_manual(values=c('not missing'='orange','MCAR'='skyblue2','MNAR'='grey2')) +
  labs(x='effect size',y='absolute percentage estimation error',) +
  scale_x_continuous(breaks = c(0,0.1,0.3,0.5,1,2,3),
                     labels = c(0,0.1,0.3,0.5,1,2,3))+
  scale_color_manual(values=c('HBSTS'='orange','HS-BSTS'='skyblue2','SS-BSTS'='grey2'))




###########################################
## when having larger random errors: ######
###########################################

color_palette <- brewer.pal(10,'Set3')
vid_color <- viridis(10,option='turbo')
mcar1 <- readRDS('simulated data/sim200_moreerror/sim200MCARobs_new.RData')
mcar2 <- readRDS('simulated data/sim200_moreerror/sim200MCARobs_err02_new.RData')
mcar3 <- readRDS('simulated data/sim200_moreerror/sim200MCARobs_err03_new.RData')
mcar4 <- readRDS('simulated data/sim200_moreerror/sim200MCARobs_err04_new.RData')
mcar5 <- readRDS('simulated data/sim200_moreerror/sim200MCARobs_err05_new.RData')
mcar6 <- readRDS('simulated data/sim200_moreerror/sim200MCARobs_err06_new.RData')
mcar7 <- readRDS('simulated data/sim200_moreerror/sim200MCARobs_err07_new.RData')
mcar8 <- readRDS('simulated data/sim200_moreerror/sim200MCARobs_err08_new.RData')
mcar9 <- readRDS('simulated data/sim200_moreerror/sim200MCARobs_err09_new.RData')
mcar10 <- readRDS('simulated data/sim200_moreerror/sim200MCARobs_err1_new.RData')


## -- APEE -- e.g. when effect size = 3##
apee1 <-abs(mcar1$avgeffectall[,1,7] -  mcar1$realeffect[,1,7])/mcar1$realeffect[,1,7]
apee2 <-abs(mcar2$avgeffectall[,1,7] -  mcar2$realeffect[,1,7])/mcar2$realeffect[,1,7]
apee3 <- abs(mcar3$avgeffectall[,1,7] -  mcar3$realeffect[,1,7])/mcar3$realeffect[,1,7]
apee4 <- abs(mcar4$avgeffectall[,1,7] -  mcar4$realeffect[,1,7])/mcar4$realeffect[,1,7]
apee5 <- abs(mcar5$avgeffectall[,1,7] -  mcar5$realeffect[,1,7])/mcar5$realeffect[,1,7]
apee6 <- abs(mcar6$avgeffectall[,1,7] -  mcar6$realeffect[,1,7])/mcar6$realeffect[,1,7]
apee7 <- abs(mcar7$avgeffectall[,1,7] -  mcar7$realeffect[,1,7])/mcar7$realeffect[,1,7]
apee8 <- abs(mcar8$avgeffectall[,1,7] -  mcar8$realeffect[,1,7])/mcar8$realeffect[,1,7]
apee9 <- abs(mcar9$avgeffectall[,1,7] -  mcar9$realeffect[,1,7])/mcar9$realeffect[,1,7]
apee10 <- abs(mcar10$avgeffectall[,1,7] -  mcar10$realeffect[,1,7])/mcar10$realeffect[,1,7]

apee1hs <-abs(mcar1$avgeffectall[,2,7] -  mcar1$realeffect[,1,7])/mcar1$realeffect[,1,7]
apee2hs <-abs(mcar2$avgeffectall[,2,7] -  mcar2$realeffect[,1,7])/mcar2$realeffect[,1,7]
apee3hs <- abs(mcar3$avgeffectall[,2,7] -  mcar3$realeffect[,1,7])/mcar3$realeffect[,1,7]
apee4hs <- abs(mcar4$avgeffectall[,2,7] -  mcar4$realeffect[,1,7])/mcar4$realeffect[,1,7]
apee5hs <- abs(mcar5$avgeffectall[,2,7] -  mcar5$realeffect[,1,7])/mcar5$realeffect[,1,7]
apee6hs <- abs(mcar6$avgeffectall[,2,7] -  mcar6$realeffect[,1,7])/mcar6$realeffect[,1,7]
apee7hs <- abs(mcar7$avgeffectall[,2,7] -  mcar7$realeffect[,1,7])/mcar7$realeffect[,1,7]
apee8hs <- abs(mcar8$avgeffectall[,2,7] -  mcar8$realeffect[,1,7])/mcar8$realeffect[,1,7]
apee9hs <- abs(mcar9$avgeffectall[,2,7] -  mcar9$realeffect[,1,7])/mcar9$realeffect[,1,7]
apee10hs <- abs(mcar10$avgeffectall[,2,7] -  mcar10$realeffect[,1,7])/mcar10$realeffect[,1,7]

apee1ss <-abs(mcar1$avgeffectall[,3,7] -  mcar1$realeffect[,1,7])/mcar1$realeffect[,1,7]
apee2ss <-abs(mcar2$avgeffectall[,3,7] -  mcar2$realeffect[,1,7])/mcar2$realeffect[,1,7]
apee3ss <- abs(mcar3$avgeffectall[,3,7] -  mcar3$realeffect[,1,7])/mcar3$realeffect[,1,7]
apee4ss <- abs(mcar4$avgeffectall[,3,7] -  mcar4$realeffect[,1,7])/mcar4$realeffect[,1,7]
apee5ss <- abs(mcar5$avgeffectall[,3,7] -  mcar5$realeffect[,1,7])/mcar5$realeffect[,1,7]
apee6ss <- abs(mcar6$avgeffectall[,3,7] -  mcar6$realeffect[,1,7])/mcar6$realeffect[,1,7]
apee7ss <- abs(mcar7$avgeffectall[,3,7] -  mcar7$realeffect[,1,7])/mcar7$realeffect[,1,7]
apee8ss <- abs(mcar8$avgeffectall[,3,7] -  mcar8$realeffect[,1,7])/mcar8$realeffect[,1,7]
apee9ss <- abs(mcar9$avgeffectall[,3,7] -  mcar9$realeffect[,1,7])/mcar9$realeffect[,1,7]
apee10ss <- abs(mcar10$avgeffectall[,3,7] -  mcar10$realeffect[,1,7])/mcar10$realeffect[,1,7]


apeedf <- data.frame(apee.mean = c(mean(apee1),mean(apee2),mean(apee3),mean(apee4),mean(apee5),
                                   mean(apee6),mean(apee7),mean(apee8),mean(apee9),mean(apee10),
                                   mean(apee1hs),mean(apee2hs),mean(apee3hs),mean(apee4hs),mean(apee5hs),
                                   mean(apee6hs),mean(apee7hs),mean(apee8hs),mean(apee9hs),mean(apee10hs),
                                   mean(apee1ss),mean(apee2ss),mean(apee3ss),mean(apee4ss),mean(apee5ss),
                                   mean(apee6ss),mean(apee7ss),mean(apee8ss),mean(apee9ss),mean(apee10ss)),
                     apee.upper = c(quantile(apee1,0.975),quantile(apee2,0.975),quantile(apee3,0.975),
                                    quantile(apee4,0.975),quantile(apee5,0.975),quantile(apee6,0.975),
                                    quantile(apee7,0.975),quantile(apee8,0.975),quantile(apee9,0.975),
                                    quantile(apee10,0.975),
                                    quantile(apee1hs,0.975),quantile(apee2hs,0.975),quantile(apee3hs,0.975),
                                    quantile(apee4hs,0.975),quantile(apee5hs,0.975),quantile(apee6hs,0.975),
                                    quantile(apee7hs,0.975),quantile(apee8hs,0.975),quantile(apee9hs,0.975),
                                    quantile(apee10hs,0.975),
                                    quantile(apee1ss,0.975),quantile(apee2ss,0.975),quantile(apee3ss,0.975),
                                    quantile(apee4ss,0.975),quantile(apee5ss,0.975),quantile(apee6ss,0.975),
                                    quantile(apee7ss,0.975),quantile(apee8ss,0.975),quantile(apee9ss,0.975),
                                    quantile(apee10ss,0.975)),
                     apee.lower =c(quantile(apee1,0.025),quantile(apee2,0.025),quantile(apee3,0.025),
                                   quantile(apee4,0.025),quantile(apee5,0.025),quantile(apee6,0.025),
                                   quantile(apee7,0.025),quantile(apee8,0.025),quantile(apee9,0.025),
                                   quantile(apee10,0.025),
                                   quantile(apee1hs,0.025),quantile(apee2hs,0.025),quantile(apee3hs,0.025),
                                   quantile(apee4hs,0.025),quantile(apee5hs,0.025),quantile(apee6hs,0.025),
                                   quantile(apee7hs,0.025),quantile(apee8hs,0.025),quantile(apee9hs,0.025),
                                   quantile(apee10hs,0.025),
                                   quantile(apee1ss,0.025),quantile(apee2ss,0.025),quantile(apee3ss,0.025),
                                   quantile(apee4ss,0.025),quantile(apee5ss,0.025),quantile(apee6ss,0.025),
                                   quantile(apee7ss,0.025),quantile(apee8ss,0.025),quantile(apee9ss,0.025),
                                   quantile(apee10ss,0.025)),
                     error = rep(seq(0.1,1,by=0.1),3),
                     model = rep(c('HBSTS','HS-BSTS','SS-BSTS'),each=10))

ggplot(apeedf,aes(x=error,y=apee.mean,group=model,color=model)) + 
  geom_point(size=2,aes(color=model)) + 
  geom_errorbar(aes(ymin=apee.lower,ymax=apee.upper),width=.02,
                position=position_dodge(0.02),size=0.5)+
  theme_bw() + theme(axis.title = element_text(size = 20),
                     axis.text= element_text(size=16),
                     legend.position=c(0.1,0.85),
                     legend.text= element_text(size=12),
                     legend.title = element_text(size=14),
                     axis.text.x = element_text(angle = 55, hjust = 1),
                     panel.grid.minor = element_blank(),
                     strip.text = element_text(size = 12))+
  #scale_color_manual(values=c('not missing'='orange','MCAR'='skyblue2','MNAR'='grey2')) +
  labs(x=expression(sigma[mu]),y='absolute percentage estimation error') +
  scale_x_continuous(breaks = seq(0.1,1,by=0.1),
                     labels = seq(0.1,1,by=0.1))+
  scale_y_continuous(breaks = seq(0.1,1,by=0.1),
                     labels = seq(0.1,1,by=0.1),
                     limits = c(0,1))+
  scale_color_manual(values=c('HBSTS'='orange','HS-BSTS'='skyblue2','SS-BSTS'='grey2'))


## -- C.I. -- ##
mcarci <- data.frame(error = rep(seq(0.1,1,0.1),each=400),
                     width = c(mcar1$avgeffectall.ci.upper[,1,2]-mcar1$avgeffectall.ci.lower[,1,2],
                               mcar2$avgeffectall.ci.upper[,1,2]-mcar2$avgeffectall.ci.lower[,1,2],
                               mcar3$avgeffectall.ci.upper[,1,2]-mcar3$avgeffectall.ci.lower[,1,2],
                               mcar4$avgeffectall.ci.upper[,1,2]-mcar4$avgeffectall.ci.lower[,1,2],
                               mcar5$avgeffectall.ci.upper[,1,2]-mcar5$avgeffectall.ci.lower[,1,2],
                               mcar6$avgeffectall.ci.upper[,1,2]-mcar6$avgeffectall.ci.lower[,1,2],
                               mcar7$avgeffectall.ci.upper[,1,2]-mcar7$avgeffectall.ci.lower[,1,2],
                               mcar8$avgeffectall.ci.upper[,1,2]-mcar8$avgeffectall.ci.lower[,1,2],
                               mcar9$avgeffectall.ci.upper[,1,2]-mcar9$avgeffectall.ci.lower[,1,2],
                               mcar10$avgeffectall.ci.upper[,1,2]-mcar10$avgeffectall.ci.lower[,1,2]))

mcarci$error <- as.factor(mcarci$error)

ggplot(mcarci,aes(x=width,fill = error,color=error)) + 
  geom_density(alpha=0.2,size=1) + 
  theme_minimal() + 
  labs(title = '',x = 'C.I. width',y= 'density') + 
  scale_fill_manual(values = color_palette,expression(paste(sigma[mu]))) +
  scale_color_manual(values = color_palette,expression(paste(sigma[mu]))) +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        legend.position = c(0.85,0.7)) + 
  scale_x_continuous(breaks = seq(0,18,3),
                     labels = seq(0,18,3))+
  scale_y_continuous(breaks = seq(0,6,1),
                     labels = seq(0,6,1),limits = c(0,6))



## -- SD --##
plot(x=seq(0.1,1,0.1),y=c(sd(mcar1$avgeffectall[,1,1]),
                          sd(mcar2$avgeffectall[,1,2]),
                          sd(mcar3$avgeffectall[,1,2]),
                          sd(mcar4$avgeffectall[,1,2]),
                          sd(mcar5$avgeffectall[,1,2]),
                          sd(mcar6$avgeffectall[,1,2]),
                          sd(mcar7$avgeffectall[,1,2]),
                          sd(mcar8$avgeffectall[,1,2]),
                          sd(mcar9$avgeffectall[,1,2]),
                          sd(mcar10$avgeffectall[,1,2])),type='l',
     xlab=expression(sigma[mu]),ylab='sd of estimated average treatment effect',
     ylim=c(0.2,10),cex.lab=1.5,cex.axis = 1.5,lwd=2)
lines(x=seq(0.1,1,0.1),y=c(sd(mcar1$avgeffectall[,2,2]),
                           sd(mcar2$avgeffectall[,2,2]),
                           sd(mcar3$avgeffectall[,2,2]),
                           sd(mcar4$avgeffectall[,2,2]),
                           sd(mcar5$avgeffectall[,2,2]),
                           sd(mcar6$avgeffectall[,2,2]),
                           sd(mcar7$avgeffectall[,2,2]),
                           sd(mcar8$avgeffectall[,2,2]),
                           sd(mcar9$avgeffectall[,2,2]),
                           sd(mcar10$avgeffectall[,2,2])),col=2,lwd=2)
lines(x=seq(0.1,1,0.1),y=c(sd(mcar1$avgeffectall[,3,2]),
                           sd(mcar2$avgeffectall[,3,2]),
                           sd(mcar3$avgeffectall[,3,2]),
                           sd(mcar4$avgeffectall[,3,2]),
                           sd(mcar5$avgeffectall[,3,2]),
                           sd(mcar6$avgeffectall[,3,2]),
                           sd(mcar7$avgeffectall[,3,2]),
                           sd(mcar8$avgeffectall[,3,2]),
                           sd(mcar9$avgeffectall[,3,2]),
                           sd(mcar10$avgeffectall[,3,2])),col=4,lwd=2)

lines(x=seq(0.1,1,0.1),y=c(sd(mcar1$avgeffectall[,1,7]),
                           sd(mcar2$avgeffectall[,1,7]),
                           sd(mcar3$avgeffectall[,1,7]),
                           sd(mcar4$avgeffectall[,1,7]),
                           sd(mcar5$avgeffectall[,1,7]),
                           sd(mcar6$avgeffectall[,1,7]),
                           sd(mcar7$avgeffectall[,1,7]),
                           sd(mcar8$avgeffectall[,1,7]),
                           sd(mcar9$avgeffectall[,1,7]),
                           sd(mcar10$avgeffectall[,1,7])),col=1,lty=2,lwd=2)
lines(x=seq(0.1,1,0.1),y=c(sd(mcar1$avgeffectall[,2,7]),
                           sd(mcar2$avgeffectall[,2,7]),
                           sd(mcar3$avgeffectall[,2,7]),
                           sd(mcar4$avgeffectall[,2,7]),
                           sd(mcar5$avgeffectall[,2,7]),
                           sd(mcar6$avgeffectall[,2,7]),
                           sd(mcar7$avgeffectall[,2,7]),
                           sd(mcar8$avgeffectall[,2,7]),
                           sd(mcar9$avgeffectall[,2,7]),
                           sd(mcar10$avgeffectall[,2,7])),col=2,lty=2,lwd=2)
lines(x=seq(0.1,1,0.1),y=c(sd(mcar1$avgeffectall[,3,7]),
                           sd(mcar2$avgeffectall[,3,7]),
                           sd(mcar3$avgeffectall[,3,7]),
                           sd(mcar4$avgeffectall[,3,7]),
                           sd(mcar5$avgeffectall[,3,7]),
                           sd(mcar6$avgeffectall[,3,7]),
                           sd(mcar7$avgeffectall[,3,7]),
                           sd(mcar8$avgeffectall[,3,7]),
                           sd(mcar9$avgeffectall[,3,7]),
                           sd(mcar10$avgeffectall[,3,7])),col=4,lty=2,lwd=2)
text(0.7,1,'effect size=0.1',cex=1.4, pos=3)
text(0.7,6,'effect size=3',cex=1.4, pos=3)














##




     