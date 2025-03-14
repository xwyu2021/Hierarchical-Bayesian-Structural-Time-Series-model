//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
  int<lower=1> T; // Number of timesteps
  int<lower=0> T_forecast;
  int<lower=0> N; //total number of observations
  int<lower=0> N_forecast;
  int<lower=0> M; // number of control TSs
  real y[N]; // for real data analysis
  int index_obs_t[T+1]; // number of observations at each time t, include the last index = N+1
  int<lower=0> K; // total number of observations for x
  real x_obs[K]; // covariates obs
  int index_x_t[T+T_forecast,M];
  int index_obs_forecast[T_forecast+1];
  real<lower=0> sigma_itcpt; // sigma of intercept beta0
  real<lower=0> scale_global;
  real<lower=1> nu_global; # d.f. for half-t prior for tau
  real<lower=1> nu_local; # d.f. for half-t prior for lambdas
  real<lower=0> slab_scale; # slab scale for the regularised horseshoe
  real<lower=0> slab_df; # slab d.f. for regularised horseshoe
  real<lower=0> sdy;
  real<lower=0> c_df;
  real<lower=0> tr_df;
}

parameters {
  vector[T] mu_err; // error of trend
  vector[T] delta_err; // error of trend
  real<lower=0> sigma2_mu; // sigma^2 of error of trend
  real<lower=0> sigma2_delta;
  real<lower=-1,upper=1> phi; // parameter of ar1
  real D;
  real<lower=0> sigma2_y; // sigma^2 of y
  real<lower=0> sigma2_obs; // sigma^2 of yit
  matrix<lower=0,upper=36>[T+T_forecast,M] x; // est mean of obs control series
  vector[M] z; // normal(0,1) r.v. for beta
  real<lower=0> aux1_global;
  real<lower=0> aux2_global;
  vector<lower=0>[M] aux1_local;
  vector<lower=0>[M] aux2_local;
  real<lower=0> caux; # unscaled c
  real beta_aux;
  real<lower=0,upper=36> alpha[T];//+T_forecast]; // mean
  real<lower=-1,upper=1> kappa[M];
  real<lower=0> sigma2_xx;//[M];
  //real<lower=0> sigma2_alpha;
}

transformed parameters {
  real<lower=0> sigma_xx;//[M];
  sigma_xx=sqrt(sigma2_xx);
  real<lower=0> sigma_y; # noise std
  sigma_y = sqrt(sigma2_y)*slab_scale;
  real<lower=0,upper=1> sigma_mu;
  sigma_mu = sqrt(sigma2_mu)*sdy;
  real<lower=0,upper=1> sigma_delta;
  sigma_delta = sqrt(sigma2_delta)*sdy;
  real<lower=0> sigma_obs;
  sigma_obs = sqrt(sigma2_obs);
  vector[T] mu; 
  mu[1] = mu_err[1]*sigma_mu;
  vector[T] delta; 
  delta[1] = delta_err[1]*sigma_delta;
  for (t in 1:(T-1)) {
    mu[t+1] = mu[t] + delta[t] + (mu_err[t]* sigma_mu);
    delta[t+1] = D + phi *(delta[t] - D) + delta_err[t] * sigma_delta;
  }
  real<lower=0> tau; // global shrinkage
  vector<lower=0>[M] lambda; // local shrinkage
  vector<lower=0>[M] lambda_tilde;
  real<lower=0> c; // slab scale
  lambda = aux1_local .* sqrt(aux2_local);
  tau = aux1_global * sqrt(aux2_global) * scale_global * sigma_y;
  c = sqrt(caux)*slab_scale;
  lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2*square(lambda)));//+ (sigma_y^2*varx)));
  vector[M] beta; // regression coeffs
  beta = z .* lambda_tilde*tau;
  real beta0;
  beta0 = beta_aux * sigma_itcpt;
 
  real x_mean[M];
  real<lower=0> x_sd[M];
  matrix[T+T_forecast,M] x_scale;
  for (k in 1:M) {
    x_mean[k] = mean(x[1:T, k]);
    x_sd[k] = sd(x[1:T,k]);
    x_scale[,k] = (x[,k]-x_mean[k])/x_sd[k];
  }
  vector[T] f; // regression component
  f = x_scale[1:T,]*beta+beta0;
  vector[M] rho;
  for(i in 1:M){
    if((1/(beta[i]^2)) >= (abs(beta[i]) * sum(x_scale[,i]^2))){
      rho[i] = 0;
    }else{
      rho[i] = abs(beta[i]) - 1/(beta[i]^2 * sum(x_scale[,i]^2));
      if(beta[i] < 0){
         rho[i] = -rho[i];
      }
    }
  }
  real alpha_mean;
  alpha_mean = mean(alpha);
  real<lower=0> alpha_sd;
  alpha_sd = sd(alpha);
  real alpha_scale[T];
  for(t in 1:T){
     alpha_scale[t]=(alpha[t] - alpha_mean)/alpha_sd; 
  }
}

model {
  // priors
  sigma2_y ~ inv_gamma(0.5*slab_df,0.5*slab_df);
  sigma2_mu ~ inv_gamma(tr_df/2,tr_df/2);
  sigma2_delta ~ inv_gamma(tr_df/2,tr_df/2);
  sigma2_obs ~ inv_gamma(0.001,0.001);
  sigma2_xx ~ inv_gamma(tr_df/2,tr_df/2);
 

  mu_err ~ std_normal();
  delta_err ~ std_normal();
  phi ~ std_normal();
  D ~ std_normal();
  
  kappa ~ std_normal();
  z ~ std_normal();
  aux1_local ~ std_normal();
  aux2_local ~ inv_gamma(0.5*nu_local,0.5*nu_local);
  aux1_global ~ std_normal();
  aux2_global ~ inv_gamma(0.5*nu_global,0.5*nu_global);
  caux ~ inv_gamma(0.5*slab_df,0.5*slab_df);
  beta_aux ~ std_normal();

  
  for(i in 1:(T+T_forecast)){
    for(m in 1:M){
      if(m == 1){
        if(i==1){
          x_obs[1:index_x_t[1,1]] ~ normal(x[1,1],sigma_xx);
        }else{
          x_obs[(index_x_t[i-1,M]+1):index_x_t[i,1]] ~ normal(x[i,1],sigma_xx);
        }
      }else{
        x_obs[(index_x_t[i,m-1]+1):index_x_t[i,m]] ~ normal(x[i,m],sigma_xx);
      }
    }
  }
  
  alpha ~ cauchy(0,1);
  for(i in 1:T){
    y[index_obs_t[i]:(index_obs_t[i+1]-1)] ~ normal(alpha[i],sigma_obs);
  }
  for(i in 1:T){
    alpha_scale[i] ~ normal(mu[i] + f[i], sigma_y);
  }
}

generated quantities {

  real mu_forecast[T_forecast];
  mu_forecast[1] = normal_rng(mu[T]+ delta[T],sigma_mu);
  real delta_forecast[T_forecast];
  delta_forecast[1] = normal_rng(D+phi*(delta[T]-D),sigma_delta);
  for(i in 2:T_forecast){
    delta_forecast[i] = normal_rng(D+phi*(delta_forecast[i-1]-D),sigma_delta);
    mu_forecast[i] = normal_rng(mu_forecast[i-1] + delta_forecast[i-1],sigma_mu);
  }
  real f_forecast[T_forecast];
  real alpha_forecast[T_forecast];
  real alpha_forecast_unscale[T_forecast];
  for(i in 1:T_forecast){
    f_forecast[i] = beta0 + (x_scale[i+T,]*rho);
    alpha_forecast[i] = mu_forecast[i] + normal_rng(f_forecast[i],sigma_y);
    alpha_forecast_unscale[i] = (alpha_forecast[i] * alpha_sd) + alpha_mean;
  }
  real log_lik_y[N];
  for(i in 1:T){
    for(j in index_obs_t[i]:(index_obs_t[i+1]-1)){
      log_lik_y[j]= normal_lpdf(y[j]|alpha[i],sigma_obs);
    }
  }
}






