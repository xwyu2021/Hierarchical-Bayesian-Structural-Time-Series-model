###########################
### load processed data ###
###########################

## raw data are available from https://www.understandingsociety.ac.uk/
## missing data imputation following Jeffery et al. 2024
## this dataset contains more covariates than data1 (finaldata_weight.RData)

datax <- readRDS("~/Desktop/LSE/data/finaldata_imputed.RData")
aMat <- readRDS("aMat_england.rds") # this dataset require special liscience from understanding society
datax$LAD21NM <- as.factor(datax$LAD21NM)
aMatNumber <- data.frame(LAD21NM = rownames(aMat), LADNM21_id = 1:309)
datax$aMatNumber <- aMatNumber$LADNM21_id[match(datax$LAD21NM,aMatNumber$LAD21NM)]

datax$policyDate <-  as.Date("2014-05-14")
datax$mediaDate <- as.Date("2017-11-28") 
datax <- datax %>%  dplyr::mutate(quater = (lubridate::interval(start = policyDate, end = interviewDate) %/% months(3)))
datax <- datax %>% dplyr::mutate(quaterSince = if_else(policyDate <= interviewDate,quater + 1,0))
datax <- datax %>%  dplyr::mutate(quater2 = (lubridate::interval(start = mediaDate, end = interviewDate) %/% months(3)))
datax <- datax %>% dplyr::mutate(quaterSince2 = ifelse(mediaDate <= interviewDate,quater2+1,0))
datax$quater <- datax$quater - min(datax$quater) + 1
datax$quater <- as.factor(datax$quater)
datax$quaterSince <- as.factor(datax$quaterSince)
datax$quaterSince2 <- as.factor(datax$quaterSince2)
datax$quater <- as.numeric(datax$quater)
datax$quaterSince <- as.numeric(datax$quaterSince)
datax$quaterSince2 <- as.numeric(datax$quaterSince2)

datax <- datax %>% mutate(quater.carib = case_when(eth3 == "carib" ~ quater, TRUE ~ 1))
datax <- datax %>% mutate(quater.afric = case_when(eth3 == "afric" ~ quater, TRUE ~ 1))
datax <- datax %>% mutate(quater.india = case_when(eth3 == "india" ~ quater, TRUE ~ 1))
datax <- datax %>% mutate(quater.pakis = case_when(eth3 == "pakis" ~ quater, TRUE ~ 1))
datax <- datax %>% mutate(quater.bangl = case_when(eth3 == "bangl" ~ quater, TRUE ~ 1))
datax <- datax %>% mutate(quater.mixed = case_when(eth3 == "mixed" ~ quater, TRUE ~ 1))
datax <- datax %>% mutate(quater.asian = case_when(eth3 == "asian" ~ quater, TRUE ~ 1))
datax <- datax %>% mutate(quater.other = case_when(eth3 == "other" ~ quater, TRUE ~ 1))

datax$quater.carib <- as.factor(datax$quater.carib)
datax$quater.afric <- as.factor(datax$quater.afric)
datax$quater.india <- as.factor(datax$quater.india)
datax$quater.pakis <- as.factor(datax$quater.pakis)
datax$quater.bangl <- as.factor(datax$quater.bangl)
datax$quater.mixed <- as.factor(datax$quater.mixed)
datax$quater.asian <- as.factor(datax$quater.asian)
datax$quater.other <- as.factor(datax$quater.other)

datax$quater.carib <- as.numeric(datax$quater.carib)
datax$quater.afric <- as.numeric(datax$quater.afric)
datax$quater.india <- as.numeric(datax$quater.india)
datax$quater.pakis <- as.numeric(datax$quater.pakis)
datax$quater.bangl <- as.numeric(datax$quater.bangl)
datax$quater.mixed <- as.numeric(datax$quater.mixed)
datax$quater.asian <- as.numeric(datax$quater.asian)
datax$quater.other <- as.numeric(datax$quater.other)


datax <- datax %>% mutate(quaterSince.carib = case_when(eth3 == "carib" ~ quaterSince, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince.afric = case_when(eth3 == "afric" ~ quaterSince, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince.india = case_when(eth3 == "india" ~ quaterSince, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince.pakis = case_when(eth3 == "pakis" ~ quaterSince, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince.bangl = case_when(eth3 == "bangl" ~ quaterSince, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince.mixed = case_when(eth3 == "mixed" ~ quaterSince, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince.asian = case_when(eth3 == "asian" ~ quaterSince, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince.other = case_when(eth3 == "other" ~ quaterSince, TRUE ~ 1))

datax$quaterSince.carib <- as.factor(datax$quaterSince.carib)
datax$quaterSince.afric <- as.factor(datax$quaterSince.afric)
datax$quaterSince.india <- as.factor(datax$quaterSince.india)
datax$quaterSince.pakis <- as.factor(datax$quaterSince.pakis)
datax$quaterSince.bangl <- as.factor(datax$quaterSince.bangl)
datax$quaterSince.mixed <- as.factor(datax$quaterSince.mixed)
datax$quaterSince.asian <- as.factor(datax$quaterSince.asian)
datax$quaterSince.other <- as.factor(datax$quaterSince.other)

datax$quaterSince.carib <- as.numeric(datax$quaterSince.carib)
datax$quaterSince.afric <- as.numeric(datax$quaterSince.afric)
datax$quaterSince.india <- as.numeric(datax$quaterSince.india)
datax$quaterSince.pakis <- as.numeric(datax$quaterSince.pakis)
datax$quaterSince.bangl <- as.numeric(datax$quaterSince.bangl)
datax$quaterSince.mixed <- as.numeric(datax$quaterSince.mixed)
datax$quaterSince.asian <- as.numeric(datax$quaterSince.asian)
datax$quaterSince.other <- as.numeric(datax$quaterSince.other)


### collect datax from imputed and start from here!
datax <- datax %>% mutate(quaterSince2.carib = case_when(eth3 == "carib" ~ quaterSince2, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince2.afric = case_when(eth3 == "afric" ~ quaterSince2, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince2.india = case_when(eth3 == "india" ~ quaterSince2, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince2.pakis = case_when(eth3 == "pakis" ~ quaterSince2, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince2.bangl = case_when(eth3 == "bangl" ~ quaterSince2, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince2.asian = case_when(eth3 == "asian" ~ quaterSince2, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince2.other = case_when(eth3 == "other" ~ quaterSince2, TRUE ~ 1))
datax <- datax %>% mutate(quaterSince2.mixed = case_when(eth3 == "mixed" ~ quaterSince2, TRUE ~ 1))

datax$quaterSince2.carib <- as.factor(datax$quaterSince2.carib)
datax$quaterSince2.afric <- as.factor(datax$quaterSince2.afric)
datax$quaterSince2.india <- as.factor(datax$quaterSince2.india)
datax$quaterSince2.pakis <- as.factor(datax$quaterSince2.pakis)
datax$quaterSince2.bangl <- as.factor(datax$quaterSince2.bangl)
datax$quaterSince2.asian <- as.factor(datax$quaterSince2.asian)
datax$quaterSince2.other <- as.factor(datax$quaterSince2.other)
datax$quaterSince2.mixed <- as.factor(datax$quaterSince2.mixed)

datax$quaterSince2.carib <- as.numeric(datax$quaterSince2.carib)
datax$quaterSince2.afric <- as.numeric(datax$quaterSince2.afric)
datax$quaterSince2.india <- as.numeric(datax$quaterSince2.india)
datax$quaterSince2.pakis <- as.numeric(datax$quaterSince2.pakis)
datax$quaterSince2.bangl <- as.numeric(datax$quaterSince2.bangl)
datax$quaterSince2.asian <- as.numeric(datax$quaterSince2.asian)
datax$quaterSince2.other <- as.numeric(datax$quaterSince2.other)
datax$quaterSince2.mixed <- as.numeric(datax$quaterSince2.mixed)

pc.u <- 1; pc.alpha <- 0.01; pc.u.phi <- 0.5; pc.alpha.phi <- 2/3
hyper.pc <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
hyper.pc.space <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)), 
                       phi = list(prior = 'pc', param = c(pc.u.phi, pc.alpha.phi)))
hyper_pc_time <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
hyper_pc_space_struc <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)), phi = list(prior = 'pc', param = c(pc.u.phi, pc.alpha.phi)))

## compare a group vs white group: e.g. caribbean vs white
subdatax <- datax[datax$eth3 %in% c('aawhit','carib'),]
subdatax <- subdatax%>% dplyr::mutate(eth4 = dplyr::case_when(eth3 %in% c('aawhit') ~ 'White',
                                                               eth3 %in% c('carib') ~ 'carib',
                                                               TRUE ~ 'mixed'))
subdatax$eth4 <- as.factor(subdatax$eth4)
subdatax <- subdatax %>% dplyr::mutate(eth4 = eth4 %>% relevel(.,ref='White'))
subdatax<- subdatax[subdatax$experiod %in% c(1,2),]
subdatax$experiod <- as.numeric(subdatax$experiod)
subdatax$experiod <- as.factor(subdatax$experiod)

modelq <- inla(
  scghq1 ~ eth3 + quater + experiod + quaterSince + quaterSince2 +
    quater.carib + #quater.afric + quater.india + quater.pakis + quater.bangl + quater.asian + quater.mixed + quater.other + 
    experiod.carib + #experiod.afric + experiod.india + experiod.pakis + experiod.bangl + experiod.asian + experiod.mixed + experiod.other +
    quaterSince.carib + #quaterSince.afric + quaterSince.india + quaterSince.pakis + quaterSince.bangl + quaterSince.asian + quaterSince.mixed + quaterSince.other +
    quaterSince2.carib + #quaterSince2.afric + quaterSince2.india + quaterSince2.pakis + quaterSince2.bangl + quaterSince2.asian + quaterSince2.mixed + quaterSince2.other +
    age + sex + urban + decileIMD + nchild + 
    ukborn_imp + mastat_imp + hiqual_imp + jbstat_imp + tenure_imp + incomequin_imp + health_imp+
    f(strata, model = 'iid', hyper= hyper.pc) + 
    f(year,
      model = 'rw1',
      hyper= hyper.pc) + 
    f(aMatNumber, graph = aMat, model = 'bym2', hyper= hyper_pc_space_struc, scale.model = TRUE, adjust.for.con.comp = TRUE),
  family='gaussian',
  data = subdatax,control.compute=list(config = TRUE))
modelq$summary.fixed

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
n.sims <- 1000
samp.all <- inla.posterior.sample(n = n.sims, result = modelq, intern = TRUE)
theta.predictor.a <-
  lapply(X = samp.all,
         FUN = function(x) {x$latent[startsWith(rownames(x$latent), 'Predictor:')]}) %>%
  unlist() %>%
  matrix(., ncol = n.sims)
colnames(theta.predictor.a) <- paste0('theta:', 1:n.sims)
theta.predictor <- 
  cbind(subdatax %>% 
          dplyr::select(scghq1,eth4,quater,experiod,quaterSince,quater),#month), 
        theta.predictor.a)
before.after.temporal.data <-
  theta.predictor %>%
  dplyr::summarise(dplyr::across(dplyr::starts_with('theta:'), mean),
                   .by = c('eth4', 'quater')) %>%
  dplyr::mutate(dplyr::select(., starts_with('theta:')) %>% 
                  apply(., 1, my.summary) %>% 
                  lapply(., data.frame) %>%
                  do.call(rbind, .)) %>%
  dplyr::select(-starts_with('theta:'))
levels(before.after.temporal.data$eth4)[levels(before.after.temporal.data$eth4) == "carib"] = "Black Caribbean"
before.after.temporal.data$quater = before.after.temporal.data$quater - 22
before.after.temporal.plot <-
  ggplot2::ggplot(data = before.after.temporal.data[before.after.temporal.data$eth4 %in% c('White','Black Caribbean'),],
                  aes(x = quater, y = mean, group = eth4, colour = eth4, fill = eth4)) +
  ggplot2::geom_vline(xintercept = 0, colour = 'black', linetype = 'dashed') +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
  ggplot2::scale_colour_manual(values = c('orange2','blue3'))+
  ggplot2::scale_fill_manual(values = c('orange2','blue3'))+
  ggplot2::scale_y_continuous(limits = c(10, 14), breaks = seq(10, 14, by = 1)) +
  ggplot2::labs(x = 'time (quarters)', 
                y = 'Self reported mental ill health (GHQ-12)') +
  my.theme(legend.title = element_blank(),
           legend.position = 'bottom',
           text = element_text(size = 12)); before.after.temporal.plot








