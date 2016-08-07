data {
  int<lower=0> N;  // total number of observations
  int<lower=0> y[N];       // Poisson draws
  real<lower=0> gamma_prior_alpha;
  real<lower=0> gamma_prior_beta;
  real<lower=0> beta_prior_alpha;
  real<lower=0> beta_prior_beta;
}
parameters {
}
transformed parameters {
}
model {
  prior_gamma ~ gamma(gamma_prior_alpha, gamma_prior_beta);
  prior_beta ~ gamma(beta_prior_alpha, beta_prior_beta);
  for (n in 1:N) {
    lambda[n] ~ gamma(prior_gamma, prior_beta);
    y[n] ~ poisson(lambda[n]);
  }
}
