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
  //int<lower=0> Npost;
  int<lower=0> M; // number of control TSs
  real<lower=0,upper=36> y[N]; // for real data analysis
  //real<lower=0,upper=36> postys[Npost];
  int index_obs_t[T+1]; // number of observations at each time t, include the last index = N+1
  //matrix[T+T_forecast,M] x; // covariate mean
  int<lower=0> K; // total number of observations for x
  real<lower=0,upper=36> x_obs[K]; // covariates obs
  int index_x_t[T+T_forecast,M];
  int index_obs_forecast[T_forecast+1];
  real<lower=0> sigma_itcpt; // sigma of intercept beta0
  real<lower=0> scale_global;
  real<lower=1> nu_global; # d.f. for half-t prior for tau
  real<lower=1> nu_local; # d.f. for half-t prior for lambdas
  real<lower=0> slab_scale; # slab scale for the regularised horseshoe
  real<lower=0> slab_df; # slab d.f. for regularised horseshoe
  real<lower=0> sdy;
  //vector<lower=0>[M] sdx;
  real<lower=0> c_df;
  real<lower=0> tr_df;
  int x_strata[K]; // covariates obs
  int<lower=1> NS; // number of different stratas
  int pre_strata[N];
  int x_psu[K]; // covariates obs
  int<lower=1> NC; // number of different stratas
  int pre_psu[N];
}

parameters {
  vector[T] mu_err; // error of trend
  real<lower=0> sigma2_mu; // sigma^2 of error of trend
  real<lower=0> sigma2_y; // sigma^2 of y
  real<lower=0> sigma2_obs; // sigma^2 of yit
  real<lower=-1,upper=1> phi; // parameter of ar1
  //matrix<lower=0>[T+T_forecast,M] sigma2_x;
  //real<lower=0> sigma2_x;//[M];
  matrix[T+T_forecast,M] x; // est mean of obs control series
  vector[M] z; // normal(0,1) r.v. for beta
  real<lower=0> aux1_global;
  real<lower=0> aux2_global;
  vector<lower=0>[M] aux1_local;
  vector<lower=0>[M] aux2_local;
  real<lower=0> caux; # unscaled c
  real beta_aux;
  real alpha[T];//+T_forecast]; // mean
  real<lower=-1,upper=1> kappa[M];
  real<lower=0> sigma2_xx;//[M];
  //real<lower=0> sigma2_alpha;
  //real alpha_pred[T_forecast];
  //matrix[T+T_forecast,M] zeta_err;//[M];
  //real<lower=0> sigma2_zeta;
  //real mu_alpha;
  vector[NS] eps;
  real<lower=0> sigma2_gamma; // random effects from strata sampling
  vector[NC] eps_c;
  real<lower=0> sigma2_cluster; // random effects from cluster sampling
}

transformed parameters {
  //real<lower=0> sigma_alpha;
  //sigma_alpha = sqrt(sigma2_alpha);
  real<lower=0> sigma_gamma;
  sigma_gamma= sqrt(sigma2_gamma)*0.1;
  real<lower=0> sigma_cluster;
  sigma_cluster= sqrt(sigma2_cluster)*0.1;
  real<lower=0> sigma_xx;//[M];
  //for(i in 1:M){
    //sigma_xx[i] = sqrt(sigma2_xx[i]);
  //}
  sigma_xx=sqrt(sigma2_xx);
  real<lower=0,upper=1.2> sigma_y; # noise std
  sigma_y = sqrt(sigma2_y)*slab_scale;
  real<lower=0,upper=1> sigma_mu;
  sigma_mu = sqrt(sigma2_mu)*sdy;
  real<lower=0> sigma_obs;
  sigma_obs = sqrt(sigma2_obs);
  //real<lower=0> sigma_x;//[M];
  //sigma_x = sqrt(sigma2_x);
  vector[T] mu; 
  mu[1] = mu_err[1] * sigma_mu;
  for (t in 2:T) {
    //mu[t] = mu[t-1] + (mu_err[t] * sigma_mu);
    mu[t] = phi*mu[t-1] + (mu_err[t] * sigma_mu);
  }
  
  real<lower=0> tau; // global shrinkage
  vector<lower=0>[M] lambda; // local shrinkage
  vector<lower=0>[M] lambda_tilde;
  real<lower=0> c; // slab scale
  lambda = aux1_local .* sqrt(aux2_local);
  //lambda = aux1_local .* sqrt(aux2_local);
  tau = aux1_global * sqrt(aux2_global) * scale_global * sigma_y;
  c = slab_scale * sqrt(caux);
  //c = 2 * sqrt(caux);
  lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2*square(lambda)));//+ (sigma_y^2*varx)));
  vector[M] beta; // regression coeffs
  beta = z .* lambda_tilde*tau;
  real beta0;
  beta0 = beta_aux * sigma_itcpt;
  //matrix[T+T_forecast,M] zeta;
  //real<lower=0> sigma_zeta;
  //sigma_zeta = sqrt(sigma2_zeta);
  //for(i in 1:(T+T_forecast)){
    //for(m in 1:M){
      //if(i == 1){
        //zeta[i,m] = zeta_err[i,m];
      //}else{
        //zeta[i,m] = kappa[m] * zeta[i-1,m] + zeta_err[i,m] * sigma_zeta;
      //}
    //}
  //}
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
  vector[NC] eps_c2;
  eps_c2 = sigma_cluster * eps_c;
  vector[NS] eps2;
  eps2 = sigma_gamma * eps;
}

model {
  // priors
  sigma2_y ~ inv_gamma(0.5*slab_df,0.5*slab_df);
  sigma2_mu ~ inv_gamma(tr_df/2,tr_df/2);
  sigma2_obs ~ inv_gamma(0.001,0.001);
  //sigma2_x ~ inv_gamma(0.001,0.001);
  //sigma2_xx ~ inv_gamma(0.001,0.001);
  sigma2_xx ~ inv_gamma(tr_df/2,tr_df/2);
  //sigma2_alpha ~ inv_gamma(0.001,0.001);
  sigma2_gamma ~ inv_gamma(0.5*slab_df,0.5*slab_df);
  sigma2_cluster ~ inv_gamma(0.5*slab_df,0.5*slab_df);
  //sigma2_gamma ~ inv_gamma(0.001,0.001);
  //sigma2_cluster ~ inv_gamma(0.001,0.001);
  //sigma_gamma ~ lognormal(0,1);
  //sigma_cluster ~ lognormal(0,1);
  eps ~ std_normal();
  eps_c ~ std_normal();
  //sigma2_zeta ~ inv_gamma(0.001,0.001);
  mu_err ~ std_normal();
  phi ~ std_normal();
  kappa ~ std_normal();
  z ~ std_normal();
  aux1_local ~ std_normal();
  aux2_local ~ inv_gamma(0.5*nu_local,0.5*nu_local);
  aux1_global ~ std_normal();
  aux2_global ~ inv_gamma(0.5*nu_global,0.5*nu_global);
  caux ~ inv_gamma(0.5*c_df,0.5*c_df);
  //caux ~ inv_gamma(2,2);
  beta_aux ~ std_normal();
  for(i in 1:(T+T_forecast)){
    if(i == 1){
      x[i,] ~ normal(0,sigma_xx);
    }else{
      for(m in 1:M){
        x[i,m] ~ normal(kappa[m] * x[i-1,m],sigma_xx);//[m]);
      }
    }
  }
  
  
   for(i in 1:(T+T_forecast)){
    for(m in 1:M){
      if(m == 1){
        if(i==1){
          for(k in 1:index_x_t[1,1]){
            x_obs[k] ~ normal(x[1,1]+eps2[x_strata[k]]+eps_c2[x_psu[k]],sigma_obs);
          }
        }else{
          for(k in (index_x_t[i-1,M]+1):index_x_t[i,1]){
            x_obs[k] ~ normal(x[i,1]+eps2[x_strata[k]]+eps_c2[x_psu[k]],sigma_obs);
          }
        }
      }else{
        for(k in (index_x_t[i,m-1]+1):index_x_t[i,m]){
          x_obs[k] ~ normal(x[i,m]+eps2[x_strata[k]]+eps_c2[x_psu[k]],sigma_obs);
        }
      }
    }
  }
  alpha ~ cauchy(0,1);

  for(i in 1:T){
    for(k in index_obs_t[i]:(index_obs_t[i+1]-1)){
        y[k] ~ normal(alpha[i]+eps2[pre_strata[k]]+ eps_c2[pre_psu[k]],sigma_obs);
    }
  }
  for(i in 1:T){
    alpha_scale[i] ~ normal(mu[i] + f[i], sigma_y);
  }

}

generated quantities {
  real alpha_forecast[T_forecast];
  real mu_forecast[T_forecast];
  mu_forecast[1] = normal_rng(phi*mu[T],sigma_mu);
  for(i in 2:T_forecast){
    mu_forecast[i] = normal_rng(phi*mu_forecast[i-1],sigma_mu);
  }
  real f_forecast[T_forecast];
  real alpha_forecast_unscale[T_forecast];
  for(i in 1:T_forecast){
    f_forecast[i] = beta0 + (x_scale[i+T,]*rho);
    alpha_forecast[i] = normal_rng(mu_forecast[i] + f_forecast[i],sigma_y);
    alpha_forecast_unscale[i] = (alpha_forecast[i] * alpha_sd) + alpha_mean;
  }
  real ys[index_obs_forecast[T_forecast+1]-1];
  for(i in 1:T_forecast){
    for(j in index_obs_forecast[i]:(index_obs_forecast[i+1]-1)){
      ys[j] = normal_rng(alpha_forecast_unscale[i],sigma_obs);
      
    }
  }
}






