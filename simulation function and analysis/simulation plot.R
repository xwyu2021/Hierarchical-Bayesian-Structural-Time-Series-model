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



y_pred <- as.matrix(modar$alpha_forecast_unscale)
ppc_dens_overlay(simemp.mean.y,y_pred) # simemp.mean.y is the real observed outcome after intervention (see simulation.R)

p<- stan_dens(
  hsim$model$bsts.model,
  pars = c("sigma_mu", "sigma_y", "sigma_obs", "sigma_obsx[1]"),nrow=2
)

p + facet_wrap(
  ~parameter,nrow=2,
  labeller = as_labeller(c(
    "sigma_mu"      = "sigma[mu]",       # σ_μ
    "sigma_y"       = "sigma",        # σ_y
    "sigma_obs"     = "sigma[y]",      # σ_obs
    "sigma_obsx[1]" = "sigma[x[1]]"   # σ_obsx1
  ), label_parsed),scales='free'
) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),  # facet label font size
    axis.title = element_text(size = 14),                 # axis title size
    axis.text  = element_text(size = 12)                  # axis tick size
  )

p <- stan_trace(hsim$model$bsts.model,c("sigma_mu", "sigma_y", "sigma_obs", "sigma_obsx[1]"),nrow=2) + ylab('')

p + facet_wrap(
  ~parameter,nrow=2,
  labeller = as_labeller(c(
    "sigma_mu"      = "sigma[mu]",       # σ_μ
    "sigma_y"       = "sigma",        # σ_y
    "sigma_obs"     = "sigma[y]",      # σ_obs
    "sigma_obsx[1]" = "sigma[x[1]]"   # σ_obsx1
  ), label_parsed),scales='free'
) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),  # facet label font size
    axis.title = element_text(size = 14),                 # axis title size
    axis.text  = element_text(size = 12)                  # axis tick size
  )


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
########### sensitivity to expected model size ################
###############################################################
e1 <- readRDS('simulated_data/simulation1/smallsim100_5N1.RData')
e2 <- readRDS('simulated_data/simulation1/smallsim100_5N2.RData')
e3 <- readRDS('simulated_data/simulation1/smallsim100_5N3.RData')
e4 <- readRDS('simulated_data/simulation1/smallsim100_5N4.RData')
e5 <- readRDS('simulated_data/simulation1/smallsim100_5N5.RData')

d1 <- readRDS('simulated_data/simulation1/smallsim100_20N1.RData')
d2 <- readRDS('simulated_data/simulation1/smallsim100_20N2.RData')
d3 <- readRDS('simulated_data/simulation1/smallsim100_20N3.RData')
d4 <- readRDS('simulated_data/simulation1/smallsim100_20N4.RData')
d5 <- readRDS('simulated_data/simulation1/smallsim100_20N5.RData')
d6 <- readRDS('simulated_data/simulation1/smallsim100_20N6.RData')
d7 <- readRDS('simulated_data/simulation1/smallsim100_20N7.RData')
d8 <- readRDS('simulated_data/simulation1/smallsim100_20N8.RData')
d9 <- readRDS('simulated_data/simulation1/smallsim100_20N9.RData')
d10 <- readRDS('simulated_data/simulation1/smallsim100_20N10.RData')
d11 <- readRDS('simulated_data/simulation1/smallsim100_20N11.RData')
d12 <- readRDS('simulated_data/simulation1/smallsim100_20N12.RData')
d13 <- readRDS('simulated_data/simulation1/smallsim100_20N13.RData')
d14 <- readRDS('simulated_data/simulation1/smallsim100_20N14.RData')
d15 <- readRDS('simulated_data/simulation1/smallsim100_20N15.RData')
d16 <- readRDS('simulated_data/simulation1/smallsim100_20N16.RData')
d17 <- readRDS('simulated_data/simulation1/smallsim100_20N17.RData')
d18 <- readRDS('simulated_data/simulation1/smallsim100_20N18.RData')
d19 <- readRDS('simulated_data/simulation1/smallsim100_20N19.RData')
d20 <- readRDS('simulated_data/simulation1/smallsim100_20N20.RData')
# RMSE(y): $pred.error
# MAE: $absolute.error
# SMAPE: $SMAPE
# RMSE(beta): $beta.error
# RMSE(inclusion prob): $inc.error
df.pred <- data.frame(error = c(colMeans(e1$pred.error),
                                colMeans(e2$pred.error),
                                colMeans(e3$pred.error),
                                colMeans(e4$pred.error),
                                colMeans(e5$pred.error),
                                colMeans(d1$pred.error),
                                colMeans(d2$pred.error),
                                colMeans(d3$pred.error),
                                colMeans(d4$pred.error),
                                colMeans(d5$pred.error),
                                colMeans(d6$pred.error),
                                colMeans(d7$pred.error),
                                colMeans(d8$pred.error),
                                colMeans(d9$pred.error),
                                colMeans(d10$pred.error),
                                colMeans(d11$pred.error),
                                colMeans(d12$pred.error),
                                colMeans(d13$pred.error),
                                colMeans(d14$pred.error),
                                colMeans(d15$pred.error),
                                colMeans(d16$pred.error),
                                colMeans(d17$pred.error),
                                colMeans(d18$pred.error),
                                colMeans(d19$pred.error),
                                colMeans(d20$pred.error)),
                      Model = rep(c('HBSTS','HS-BSTS',
                                   # 'SS-BSTS'),25),
                                   'BSTS'),25),
                      size = rep(c(1:5,1:20),each=3),
                      GroupSize = c(rep(5,15),rep(20,60)))
df.pred$error  <- df.pred$error * 200
ggplot(df.pred[df.pred$Model != 'HS-BSTS',], aes(x = size, y = error,
                    linetype = factor(GroupSize),
                    group = interaction(Model, GroupSize),
                    color = Model)) + 
  geom_line(size = 0.6) +
  geom_point(size = 3.5) +
  scale_linetype_manual(values = c("5" = "solid", "20" = "dotdash")) +
  scale_color_manual(values = c('HBSTS'='orange', 
                                #'HS-BSTS'='skyblue2', 
                                #'SS-BSTS'='grey2')) +
                                'BSTS'='grey2')) +
  labs(x = 'Expected model size', 
       y = 'RMSE', 
       color = 'Model', 
       linetype = 'Number of control groups') +
  guides(color = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2, nrow = 1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 18),
        legend.position = 'top',
        legend.box = "vertical",          # stack legends vertically
        legend.box.just = "center",       # center the legends
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 18, face = "bold"))



###############################################################
### robustness and effectiveness of detecting an effect #######
###############################################################

## e.g. first experiment
nomiss <- readRDS('simulated_data/simulation2 H20/smallsim100_20.RData')

### --- plot of probability of detecting an effect --- ###
R <- dim(nomiss$pred.error)[1]
k <- cbind(colSums(1-nomiss$spec),apply(nomiss$sen,3,colSums))
n_k <- R-k
point.hbsts <- point.hsbsts <- point.ssbsts <- matrix(0,nrow = 2,ncol=8)
for(i in 1:8){
  point.hbsts[,i] <- quantile(rbeta(10000,1+k[1,i],1+n_k[1,i]),c(0.025,0.975))
  point.hsbsts[,i] <- quantile(rbeta(10000,1+k[2,i],1+n_k[2,i]),c(0.025,0.975))
  point.ssbsts[,i] <- quantile(rbeta(10000,1+k[3,i],1+n_k[3,i]),c(0.025,0.975))
}


df.nomiss <- data.frame(mean = as.vector(t(k))/R,
                        lower = c(point.hbsts[1,],point.hsbsts[1,],point.ssbsts[1,]),
                        upper = c(point.hbsts[2,],point.hsbsts[2,],point.ssbsts[2,]),
                        Model = rep(c('HBSTS','HS-BSTS','BSTS'),each=8),
                        increment = rep(c(0,0.01,0.03,0.05,0.07,0.1,0.3,0.5),3))
#df.nomiss <- df.nomiss[df.nomiss$increment != 0.01,]

ggplot(df.nomiss[df.nomiss$Model!= 'HS-BSTS',], aes(x=increment, y=mean,group=Model,color=Model)) + 
  geom_line(size=0.6) +
  geom_point(size=3.5)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.02,
                position=position_dodge(0.01),size=0.6) + 
  theme_bw() + theme(axis.title = element_text(size = 18),
                     axis.text= element_text(size=18),
                     #legend.position='none',
                     legend.position=c(0.8,0.2),
                     legend.text = element_text(size=20),
                     legend.title = element_text(size=20),
                     axis.text.x = element_text(angle = 60, hjust = 1,size=18),
                     panel.grid.minor = element_blank())+
  scale_color_manual(values=c('HBSTS'='orange',
                              #'HS-BSTS'='skyblue2',
                              'BSTS'='grey2')) +
  labs(x='Effect size',y='Proportion of intervals excluding 0') +
  scale_x_continuous(breaks = c(0,0.01,0.03,0.05,0.07,0.1,0.3,0.5),#c(0,0.1,0.3,0.5,1,2,3),
                     labels = c(0,0.01,0.03,0.05,0.07,0.1,0.3,0.5))+#c(0,0.1,0.3,0.5,1,2,3)) + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),
                     labels = c(0,0.2,0.4,0.6,0.8,1),
                     limits = c(0,1))




### --- plot of APEE --- ###
mcar <- readRDS('simulated_data/simulation2 H20/smallsim100mcar_20.RData')
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

sigma001 <- readRDS('simulated_data/simulation2 H20/smallsim100_20.RData')
sigma05 <-  readRDS('simulated_data/simulation2 H20/smallsim100mu05_20.RData')
sigma1 <-  readRDS('simulated_data/simulation2 H20/smallsim100mu1_20.RData')
sigma15 <-  readRDS('simulated_data/simulation2 H20/smallsim100mu15_20.RData')
sigma2 <-  readRDS('simulated_data/simulation2 H20/smallsim100mu2_20.RData')
apee.sigma <- data.frame(mean = c(apply(sigma001$apee,2,colMeans)[3,],
                                  apply(sigma05$apee,2,colMeans)[3,],
                                  apply(sigma1$apee,2,colMeans)[3,],
                                  apply(sigma15$apee,2,colMeans)[3,],
                                  apply(sigma2$apee,2,colMeans)[3,]),
                         lower = c(apply(sigma001$apee,2,apply,2,quantile,0.025)[3,],
                                   apply(sigma05$apee,2,apply,2,quantile,0.025)[3,],
                                   apply(sigma1$apee,2,apply,2,quantile,0.025)[3,],
                                   apply(sigma15$apee,2,apply,2,quantile,0.025)[3,],
                                   apply(sigma2$apee,2,apply,2,quantile,0.025)[3,]),
                         upper = c(apply(sigma001$apee,2,apply,2,quantile,0.975)[3,],
                                   apply(sigma05$apee,2,apply,2,quantile,0.975)[3,],
                                   apply(sigma1$apee,2,apply,2,quantile,0.975)[3,],
                                   apply(sigma15$apee,2,apply,2,quantile,0.975)[3,],
                                   apply(sigma2$apee,2,apply,2,quantile,0.975)[3,]),
                         Model =  rep(c('HBSTS','HS-BSTS','BSTS'), 5),
                         Sigma = rep(c(0.01,0.5,1,1.5,2),each=3))

ggplot(apee.sigma[apee.sigma$Model!= 'HS-BSTS',],aes(x=Sigma,y=mean,group=Model,color=Model)) + 
  geom_point(size=2,aes(color=Model)) + 
  geom_errorbar(aes(ymin=lower,ymax=upper,color=Model),width=.1,
                position=position_dodge(0.03),size=0.5)+
  theme_bw() + theme(axis.title = element_text(size = 20),
                     axis.text= element_text(size=16),
                     legend.position=c(0.15,0.8),
                     legend.text = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     axis.text.x = element_text(angle = 0, hjust = 1),
                     panel.grid.minor = element_blank(),
                     strip.text = element_text(size = 12))+
  #scale_color_manual(values=c('not missing'='orange','sigma001'='skyblue2','MNAR'='grey2')) +
  labs(x=expression(sigma[mu]),y='Absolute percentage estimation error',) +
  scale_x_continuous(breaks = c(0.01,0.5,1,1.5,2),
                     labels = c(0.01,0.5,1,1.5,2))+
  scale_color_manual(values=c('HBSTS'='orange',#'HS-BSTS'='skyblue2',
                              'BSTS'='grey2'))

t50 <-  readRDS('simulated_data/simulation2 H20/smallsim50_20.RData')
t100 <-  readRDS('simulated_data/simulation2 H20/smallsim100_20.RData')
t200 <-   readRDS('simulated_data/simulation2 H20/smallsim200_20.RData')

apee.t <- data.frame(mean = c(apply(t50$apee,2,colMeans)[1,],
                              apply(t100$apee,2,colMeans)[1,],
                              apply(t200$apee,2,colMeans)[1,]),
                     lower = c(apply(t50$apee,2,apply,2,quantile,0.025)[1,],
                               apply(t100$apee,2,apply,2,quantile,0.025)[1,],
                               apply(t200$apee,2,apply,2,quantile,0.025)[1,]),
                     upper = c(apply(t50$apee,2,apply,2,quantile,0.975)[1,],
                               apply(t100$apee,2,apply,2,quantile,0.975)[1,],
                               apply(t200$apee,2,apply,2,quantile,0.975)[1,]),
                     Model =  rep(c('HBSTS','HS-BSTS','BSTS'), 3),
                     Time = rep(c(50,100,200),each=3))


ggplot(apee.t[apee.t$Model!='HS-BSTS',],aes(x=Time,y=mean,group=Model,color=Model)) + 
  geom_point(size=2,aes(color=Model)) + 
  geom_errorbar(aes(ymin=lower,ymax=upper,color=Model),width=10,
                position=position_dodge(3),size=0.5)+
  theme_bw() + theme(axis.title = element_text(size = 20),
                     axis.text= element_text(size=16),
                     legend.position=c(0.6,0.9),
                     legend.direction = "horizontal",
                     legend.text = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     axis.text.x = element_text(angle = 0, hjust = 1),
                     panel.grid.minor = element_blank(),
                     strip.text = element_text(size = 12))+
  #scale_color_manual(values=c('not missing'='orange','sigma001'='skyblue2','MNAR'='grey2')) +
  labs(x='Length of study',y='Absolute percentage estimation error',) +
  scale_x_continuous(breaks = c(50,100,200),
                     labels = c(50,100,200))+
  scale_color_manual(values=c('HBSTS'='orange',#'HS-BSTS'='skyblue2',
                              'BSTS'='grey2'))+
  guides(color = guide_legend(nrow = 1)) 

###########################################
## when having larger random errors: ######
###########################################

color_palette <- brewer.pal(10,'Set3')
vid_color <- viridis(10,option='turbo')
#[,1,3]: hbsts;[,2,3]: hsbsts;[,3,3]: bsts
len <- c(length(sigma001$avgeffectall.ci.upper[,2,3]),
         length(sigma05$avgeffectall.ci.upper[,2,3]),
         length(sigma1$avgeffectall.ci.upper[,2,3]),
         length(sigma15$avgeffectall.ci.upper[,2,3]),
         length(sigma2$avgeffectall.ci.upper[,2,3]))
         
         
ciplot <- data.frame(error = rep(c(0.01,0.5,1,1.5,2),len),
                     width = c(sigma001$avgeffectall.ci.upper[,3,3]-sigma001$avgeffectall.ci.lower[,3,3],
                               sigma05$avgeffectall.ci.upper[,3,3]-sigma05$avgeffectall.ci.lower[,3,3],
                               sigma1$avgeffectall.ci.upper[,3,3]-sigma1$avgeffectall.ci.lower[,3,3],
                               sigma15$avgeffectall.ci.upper[,3,3]-sigma15$avgeffectall.ci.lower[,3,3],
                               sigma2$avgeffectall.ci.upper[,3,3]-sigma2$avgeffectall.ci.lower[,3,3]))

ciplot$error <- as.factor(ciplot$error)

ggplot(ciplot,aes(x=width,fill = error,color=error)) + 
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
        legend.position = c(0.85,0.7)) 














     