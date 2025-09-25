data {
  int<lower=1> T; // Number of timesteps
  int<lower=0> T_forecast;
  int<lower=0> M; // number of control TSs
  real y[T]; // observed value
  matrix[T+T_forecast,M] x; // covariate excluding intercept
  real<lower=0> sigma_itcpt; // sigma of intercept beta0
  real<lower=0> scale_global; // prior scale for the halft-t prior for tau
  real<lower=1> nu_global; // d.f. for half-t prior for tau
  real<lower=1> nu_local; // d.f. for half-t prior for lambdas
  real<lower=0> slab_scale; // slab scale for the regularised horseshoe
  real<lower=0> slab_df;// slab d.f. for regularised horseshoe
  real<lower=0> c_df;
  //cov_matrix[M] scale_beta;
  int Meff; // indicate whether specifying eff
  real<lower=0> sdy;
  //vector<lower=0>[M] sdx;
  real<lower=0> tr_df;
}

parameters {
  vector[T] mu_err; // error of trend
  //real<lower=0> sigma2_mu; // sigma^2 of error of trend
  real<lower=0> mu_aux;
  vector[M] z; // normal(0,1) r.v. for beta
  real<lower=0> aux1_global;
  real<lower=0> aux2_global;
  vector<lower=0>[M] aux1_local;
  vector<lower=0>[M] aux2_local;
  real<lower=0> caux; // unscaled c
  real beta_aux;
  real<lower=0> y_aux;
  real<lower=-1,upper=1> phi;
}

transformed parameters {
  real<lower=0,upper=1.2> sigma_y; // noise std
  vector[T] mu; 
  real<lower=0,upper=1> sigma_mu;
  real<lower=0> tau; // global shrinkage
  vector<lower=0>[M] lambda; // local shrinkage
  real<lower=0> c; // slab scale
  vector<lower=0>[M] lambda_tilde;
  vector[M] beta; // regression coeffs
  vector[T] f; // regression component
  real beta0; // intercept in regression
  vector[M] rho;
  
  sigma_y = sqrt(y_aux)*slab_scale;
  mu[1] = mu_err[1];
  
  sigma_mu = sqrt(mu_aux)*sdy;
  for (t in 2:T) {
     mu[t] = phi*mu[t-1] + (mu_err[t]* sigma_mu);
  }
  
  
  lambda = aux1_local .* sqrt(aux2_local);//./ sdx;
  if(Meff == 1){
    tau = aux1_global * sqrt(aux2_global) * scale_global * sigma_y;
  }else{
    tau = aux1_global * sqrt(aux2_global); // tau ~cauchy(0,1)
  }
  
  c = slab_scale * sqrt(caux);
  //c = sdy * sqrt(caux);
  //c = 0.1*slab_scale * sqrt(caux)/sqrt(0.2);
  
  lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2*square(lambda)));//+ (sigma_y^2*varx)));
  
  //beta = (L * z) .* lambda_tilde*tau;
  beta = z .* lambda_tilde*tau;
  //beta = z .* lambda * tau;
  
  beta0 = sigma_itcpt * beta_aux;
  
  f = beta0 + x[1:T,]*beta;
  //f = x[1:T,]*rho;
  
  for(i in 1:M){
    if((1/(beta[i]^2)) >= (fabs(beta[i]) * sum(x[,i] .* x[,i]))){
      rho[i] = 0;
    }else{
      rho[i] = fabs(beta[i]) - 1/(beta[i]^2 * sum(x[,i] .* x[,i]));
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
  phi ~ std_normal();
  z ~ std_normal();
  //z ~ multi_normal(rep_vector(0,M),scale_beta);
  // auxilary variables for shrinkage parameters and slab variance
  aux1_local ~ std_normal();
  aux2_local ~ inv_gamma(0.5*nu_local,0.5*nu_local);
  aux1_global ~ std_normal();
  aux2_global ~ inv_gamma(0.5*nu_global,0.5*nu_global);
  caux ~ inv_gamma(0.5*c_df,0.5*c_df);
  //caux ~ inv_gamma(2,2);
  beta_aux ~ std_normal();
  
  
  for(t in 1:T){
    y[t]-mu[t] ~ normal(f[t],sigma_y);
  }
  
}

generated quantities {
  real y_fitted[T];
  real y_forecast[T_forecast];
  real mu_forecast[T_forecast];
  real f_forecast[T_forecast];
  for(i in 1:T){
    y_fitted[i] = normal_rng(mu[i] + f[i],sigma_y);
  }
  
  mu_forecast[1] = normal_rng(mu[T]*phi,sigma_mu);
  
  for(i in 2:T_forecast){
    mu_forecast[i] = normal_rng(mu_forecast[i-1]*phi,sigma_mu);
  }
  for(i in 1:T_forecast){
    f_forecast[i] = beta0 + x[i+T,]*rho;
    f_forecast[i] = x[i+T,]*rho;
    y_forecast[i] = mu_forecast[i] + normal_rng(f_forecast[i],sigma_y);
  }
  
}






