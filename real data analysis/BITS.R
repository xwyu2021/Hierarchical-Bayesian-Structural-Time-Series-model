###########################
### load processed data ###
###########################

## raw data are available from https://www.understandingsociety.ac.uk/
## missing data imputation following Jeffery et al. 2024
## this dataset contains more covariates than data1 (finaldata_weight.RData)

testdf <- readRDS("~/Desktop/LSE/data/finaldata_imputed.RData")
aMat <- readRDS("~/Desktop/LSE/data/aMat_england.rds") # this dataset require special liscience from understanding society
testdf$LAD21NM <- as.factor(datax$LAD21NM)
aMatNumber <- data.frame(LAD21NM = rownames(aMat), LADNM21_id = 1:309)
testdf$aMatNumber <- aMatNumber$LADNM21_id[match(datax$LAD21NM,aMatNumber$LAD21NM)]

data1 <- readRDS('finaldata_weight.RData')
testdf$psu <- as.factor(data1$psu)

#testdf <- readRDS('final_imputed_noc.RData')
#testdf$psu <- as.factor(data1$psu)


testdf$policyDate <-  as.Date("2014-05-14")
testdf$mediaDate <- as.Date("2017-11-28") 
testdf <- testdf[testdf$interviewDate < testdf$mediaDate,]

testdf <- testdf %>%  dplyr::mutate(quater = (lubridate::interval(start = policyDate, end = interviewDate) %/% months(1)))
testdf <- testdf %>% dplyr::mutate(quaterSince = if_else(policyDate <= interviewDate,quater + 1,0))
testdf$quater <- as.numeric(testdf$quater) #  -64-42 just time counts
testdf$quaterSince <- as.factor(testdf$quaterSince) 
testdf$quaterSince <- as.numeric(testdf$quaterSince) #1:44
testdf$time <- as.factor(testdf$quater)
testdf$time <- as.numeric(testdf$time) # 1-107
testdf <- testdf %>% mutate(experiod = case_when(interviewDate < policyDate ~ "0", interviewDate >= policyDate ~ "1"))
testdf$experiod <- as.factor(testdf$experiod)
testdf$experiod <- as.numeric(testdf$experiod)
testdf$experiod <- as.factor(testdf$experiod)

datasub <- testdf[testdf$eth1 %in% c(1,14),] ## select the two groups to be compared
datasub$group <- 'treated'
datasub$group[datasub$eth1 == 1] <- 'control' ## white group is the reference

datasub <- datasub %>% mutate(experiod = case_when(interviewDate < policyDate ~ "0", interviewDate >= policyDate ~ "1"))
datasub$experiod <- as.factor(datasub$experiod)
datasub$experiod <- as.numeric(datasub$experiod) ## 1,2

datasub$experiod.carib <- datasub$experiod
datasub$experiod.carib[datasub$group == 'control'] <- 1
datasub$experiod.carib <- as.factor(datasub$experiod.carib)
datasub$experiod <- as.factor(datasub$experiod)

datasub<- datasub %>% mutate(quater.carib = case_when(group == "treated" ~ quater, TRUE ~ 1))
datasub$quater.carib <- as.factor(datasub$quater.carib)
datasub$quater.carib <- as.numeric(datasub$quater.carib)


datasub$quaterSince.carib <- datasub$quaterSince
datasub$quaterSince.carib[datasub$group == 'control'] <- 1
datasub$quaterSince.carib <- as.factor(datasub$quaterSince.carib)
datasub$quaterSince.carib <- as.numeric(datasub$quaterSince.carib)


pc.u <- 1; pc.alpha <- 0.01; pc.u.phi <- 0.5; pc.alpha.phi <- 2/3
hyper.pc <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
hyper.pc.space <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)), 
                       phi = list(prior = 'pc', param = c(pc.u.phi, pc.alpha.phi)))
hyper_pc_time <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
hyper_pc_space_struc <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)), phi = list(prior = 'pc', param = c(pc.u.phi, pc.alpha.phi)))


itsFormula.mon <- as.formula(
  scghq1 ~ group + 
    quater + experiod + quaterSince + 
    quater.carib +experiod.carib + quaterSince.carib + 
    age + sex  + decileIMD  + urban + nchild+  ukborn_imp +
    hiqual_imp + tenure_imp+  jbstat_imp +mastat_imp  + incomequin_imp + health_imp +
    f(strata, model = 'iid', hyper= hyper.pc) + 
    f(psu, model = 'iid', hyper= hyper.pc) + 
    f(time,
      model = 'rw1',
      hyper= hyper.pc) + 
    f(aMatNumber, graph = aMat, model = 'bym2', hyper= hyper_pc_space_struc, scale.model = TRUE, adjust.for.con.comp = TRUE))


my.summary <- function(x, CI = 0.95) {
  lowerCI <- (1 - CI)/2
  upperCI <- 1 - lowerCI
  
  qntl <- quantile(x, probs = c(lowerCI, 0.5, upperCI))
  data.frame(mean = mean(x), variance = var(x), lower = qntl[1], median = qntl[2], upper = qntl[3])
  
}
my.theme<-function(...){
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black'),
        # legend.title=element_blank(),
        legend.text.align=0,
        legend.key=element_rect(fill=NA),
        ...)
}




fit1 <- inla(itsFormula.mon,family='gaussian',
             data = testdf,control.compute=list(config = TRUE),control.predictor = list(compute = TRUE))

fit1$summary.fixed




data_pred <- datasub[datasub$group == 'treated' & datasub$experiod !='1',]

data_pred$quaterSince.carib <- 1
data_pred$experiod.carib <- '1'
data_pred$scghq1 <- NA 
data_extend <- rbind(datasub,data_pred)
fit2 <- inla(itsFormula.mon,family='gaussian',
             data = data_extend,control.compute=list(config = TRUE),
             control.predictor = list(compute = TRUE,link=1))
n_samples <- 1000  
samples <- inla.posterior.sample(n_samples, fit2)
start.idx <- dim(datasub)[1] + 1
end.idx <- dim(data_extend)[1]
theta.predictor.a <-
  lapply(X = samples,
         FUN = function(x) {x$latent[startsWith(rownames(x$latent), 'Predictor:')]}) %>%
  unlist() %>%
  matrix(., ncol = n_samples)
colnames(theta.predictor.a) <- paste0('theta:', 1:n_samples)
pred.samples <- theta.predictor.a[start.idx:end.idx,]
pred.original <- theta.predictor.a[1:(start.idx-1),]
pred.new <- pred.original
pred.new[datasub$group == 'treated' & datasub$experiod !='1',] <- pred.samples
new.predictor <- 
  cbind(datasub %>% 
          dplyr::select(scghq1,group,time,experiod,
                        age, sex, urban,decileIMD,
                        nchild,ukborn_imp, mastat_imp, hiqual_imp, jbstat_imp,
                        tenure_imp, incomequin_imp, health_imp,
                        strata,quater),#,,aMatNumber,month),#month), 
        pred.new)

origin.predictor <- cbind(datasub %>% 
                            dplyr::select(scghq1,group,time,experiod,
                                          age, sex, urban,decileIMD,
                                          nchild,ukborn_imp, mastat_imp, hiqual_imp, jbstat_imp,
                                          tenure_imp, incomequin_imp, health_imp,
                                          strata,quater),#,,aMatNumber,month),#month), 
                          pred.original)

all.predictor <- rbind(origin.predictor,new.predictor[datasub$group == 'treated' & datasub$experiod !='1',])
all.predictor$newgroup <- all.predictor$group
num.pred <- length(which(datasub$group == 'treated' & datasub$experiod !='1'))
all.predictor$newgroup[(dim(all.predictor)[1]-num.pred+1):dim(all.predictor)[1]] <- 'treated-counterfactual'

all.plotdata <-
  all.predictor %>%
  dplyr::summarise(dplyr::across(dplyr::starts_with('theta:'), mean),
                   .by = c('newgroup', 'time')) %>%
  dplyr::mutate(dplyr::select(., starts_with('theta:')) %>% 
                  apply(., 1, my.summary) %>% 
                  lapply(., data.frame) %>%
                  do.call(rbind, .)) %>%
  dplyr::select(-starts_with('theta:'))
all.plot<-
  ggplot2::ggplot(data = all.plotdata, aes(x = time, y = mean, group = newgroup, colour = newgroup, fill = newgroup)) +
  ggplot2::geom_vline(xintercept = 43, colour = 'black', linetype = 'dashed') +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
  ggplot2::scale_colour_manual(values = c('blue3', 'orange2',"darkgreen"))+#,"#009966")) +
  ggplot2::scale_fill_manual(values = c('blue3', 'orange2',"darkgreen"))+#,"#009966")) +
  ggplot2::labs(x = 'time (months)', 
                y = 'scaled outcome') +
  my.theme(legend.title = element_blank(),
           legend.position = 'bottom',
           text = element_text(size = 12)); all.plot

scaledall.predictor <- as.matrix(all.predictor)
scalecounter <- scaledall.predictor[(dim(all.predictor)[1]-num.pred+1):dim(all.predictor)[1],19:1018]
scalecounter <- apply(scalecounter, 2, as.numeric)
scalefit <-matrix(rep(datasub$scghq1[which(datasub$group == 'treated' & datasub$experiod !='1')],1000),
                  ncol=1000,byrow = FALSE)
scalediff <- scalefit - scalecounter

scalediff <- cbind(scaledall.predictor[(dim(all.predictor)[1]-num.pred+1):dim(all.predictor)[1],4],scalediff)
colnames(scalediff)[1] <- 'experiod'
ckdf <- scalediff[as.numeric(which(scalediff[,1] == '2')),2:dim(scalediff)[2]]
ckdf <- apply(scalediff[,2:1001], 2, as.numeric)
ckdf <- as.data.frame(ckdf)
ckdf$experiod <- scalediff[,1]

avg <- ckdf %>%
  dplyr::summarise(dplyr::across(dplyr::starts_with('theta:'), mean),
                   .by=c('experiod')) %>%
  dplyr::mutate(dplyr::select(., starts_with('theta:')) %>% 
                  apply(., 1, my.summary) %>% 
                  lapply(., data.frame) %>%
                  do.call(rbind, .)) %>%
  dplyr::select(-starts_with('theta:'))

avg
