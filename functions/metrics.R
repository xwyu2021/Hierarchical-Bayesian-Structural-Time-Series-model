###################
# metrics #########
###################

# -- rmse
RMSE <- function(x,y){
  return(sqrt(mean((x-y)^2)))
}

# -- Dickey-Savage density ratio tests
DSratio <- function(priorsd,posterior.samples){
  dnorm(0,mean=0,sd=priorsd)/mean(which(round(posterior.samples,digit=3) == 0))
}

# -- symmetric mean absolute percentage error (sMAPE
sMAPE <- function(T,predicty,obsy){
  a <- sum(abs(predicty - obsy)/((abs(obsy)-abs(predicty))/2))
  return((a*100)/T) 
}
