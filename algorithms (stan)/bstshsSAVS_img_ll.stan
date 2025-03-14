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
  int<lower=0> M; // number of control TSs
  //real<lower=0,upper=36> y[T]; // observed value
  real y[T]; // observed value
  matrix[T+T_forecast,M] x; // covariate excluding intercept
  //real<lower=0> df; //degree of freedom
  //real<lower=0> ss; //sum of squares
  //int<lower=0,upper=M> ems; //expected model size
  //real kappa; // weight
  real<lower=0> sigma_itcpt; // sigma of intercept beta0
  real<lower=0> scale_global; # prior scale for the halft-t prior for tau
  real<lower=1> nu_global; # d.f. for half-t prior for tau
  real<lower=1> nu_local; # d.f. for half-t prior for lambdas
  real<lower=0> slab_scale; # slab scale for the regularised horseshoe
  real<lower=0> slab_df; # slab d.f. for regularised horseshoe
  real<lower=0> c_df;
  //cov_matrix[M] scale_beta;
  int Meff; // indicate whether specifying eff
  real<lower=0> sdy;
  //vector<lower=0>[M] sdx;
  real<lower=0> tr_df;
}

//transformed data {
  //matrix[M,M] L;
  //L = cholesky_decompose(scale_beta);
//}

parameters {
  vector[T] mu_err; // error of trend
  //real<lower=0> sigma2_mu; // sigma^2 of error of trend
  real<lower=0> mu_aux;
  vector[M] z; // normal(0,1) r.v. for beta
  real<lower=0> aux1_global;
  real<lower=0> aux2_global;
  vector<lower=0>[M] aux1_local;
  vector<lower=0>[M] aux2_local;
  real<lower=0> caux; # unscaled c
  real beta_aux;
  real<lower=0> y_aux;
  //real<lower=0> sigma_y; # noise std
}

transformed parameters {
  real<lower=0,upper=1.2> sigma_y; # noise std
  sigma_y = sqrt(y_aux)*slab_scale;
  vector[T] mu; 
  mu[1] = mu_err[1];
  real<lower=0,upper=1> sigma_mu;
  sigma_mu = sqrt(mu_aux)*sdy;
  for (t in 2:T) {
     mu[t] = mu[t-1] + (mu_err[t]* sigma_mu);
  }
  real<lower=0> tau; // global shrinkage
  vector<lower=0>[M] lambda; // local shrinkage
  lambda = aux1_local .* sqrt(aux2_local);//./ sdx;
  if(Meff == 1){
    tau = aux1_global * sqrt(aux2_global) * scale_global * sigma_y;
  }else{
    tau = aux1_global * sqrt(aux2_global); // tau ~cauchy(0,1)
  }
  real<lower=0> c; // slab scale
  c = slab_scale * sqrt(caux);
  //c = sdy * sqrt(caux);
  //c = 0.1*slab_scale * sqrt(caux)/sqrt(0.2);
  vector<lower=0>[M] lambda_tilde;
  lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2*square(lambda)));//+ (sigma_y^2*varx)));
  vector[M] beta; // regression coeffs
  //beta = (L * z) .* lambda_tilde*tau;
  beta = z .* lambda_tilde*tau;
  //beta = z .* lambda * tau;
  real beta0; // intercept in regression
  beta0 = sigma_itcpt * beta_aux;
  vector[T] f; // regression component
  f = beta0 + x[1:T,]*beta;
  //f = x[1:T,]*rho;
  vector[M] rho;
  for(i in 1:M){
    if((1/(beta[i]^2)) >= (abs(beta[i]) * sum(x[,i]^2))){
      rho[i] = 0;
    }else{
      rho[i] = abs(beta[i]) - 1/(beta[i]^2 * sum(x[,i]^2));
      if(beta[i] < 0){
         rho[i] = -rho[i];
      }
    }
  }
  
}

model {
  // priors
  //sigma2_mu ~ inv_gamma(0.01,0.01);
  //y_aux ~ inv_gamma(0.5,0.5*5);//
  y_aux ~ inv_gamma(0.5*slab_df,0.5*slab_df);
  mu_aux ~ inv_gamma(tr_df/2,tr_df/2);// no difference from using 0.5,0.5 or 0.01,0.01
  mu_err ~ std_normal();
  z ~ std_normal();
  //z ~ multi_normal(rep_vector(0,M),scale_beta);
  // auxilary variables for shrinkage parameters and slab variance
  aux1_local ~ std_normal();
  aux2_local ~ inv_gamma(0.5*nu_local,0.5*nu_local);
  aux1_global ~ std_normal();
  aux2_global ~ inv_gamma(0.5*nu_global,0.5*nu_global);
  caux ~ inv_gamma(0.5*c_df,0.5*c_df);
  beta_aux ~ std_normal();
  
  
  for(t in 1:T){
    y[t]-mu[t] ~ normal(f[t],sigma_y);
  }
  
}

generated quantities {
  real y_fitted[T];
  for(i in 1:T){
    y_fitted[i] = normal_rng(mu[i] + f[i],sigma_y);
  }
  real y_forecast[T_forecast];
  real mu_forecast[T_forecast];
  mu_forecast[1] = normal_rng(mu[T],sigma_mu);
  real f_forecast[T_forecast];
  for(i in 2:T_forecast){
    mu_forecast[i] = normal_rng(mu_forecast[i-1],sigma_mu);
  }
  for(i in 1:T_forecast){
    f_forecast[i] = beta0 + x[i+T,]*rho;
    //f_forecast[i] = x[i+T,]*rho;
    y_forecast[i] = mu_forecast[i] + normal_rng(f_forecast[i],sigma_y);
  }
  
}






