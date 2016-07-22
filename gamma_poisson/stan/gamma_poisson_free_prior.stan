data {
  int<lower=0> N;  // total number of observations
  int<lower=0> y[N];       // Poisson draws
  real prior_gamma_mean;
  real prior_beta_mean;
  real<lower=0> prior_gamma_var;
  real<lower=0> prior_beta_var;
}
parameters {
  real<lower=0> lambda[N];
  real<lower=0> prior_gamma;
  real<lower=0> prior_beta;
}
transformed parameters {
  real<lower=0> gamma_alpha;
  real<lower=0> gamma_beta;
  real<lower=0> beta_alpha;
  real<lower=0> beta_beta;

  gamma_alpha <- (prior_gamma_mean ^ 2) / prior_gamma_var;
  gamma_beta <- prior_gamma_mean / prior_gamma_var;
  
  beta_alpha <- (prior_beta_mean ^ 2) / prior_beta_var;
  beta_beta <-prior_gamma_mean / prior_beta_var;
}
model {
  prior_gamma ~ gamma(gamma_alpha, gamma_beta);
  prior_beta ~ gamma(beta_alpha, beta_beta);
  for (n in 1:N) {
    lambda[n] ~ gamma(prior_gamma, prior_beta);
    y[n] ~ poisson(lambda[n]);
  }
}
