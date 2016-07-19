data {
  int<lower=0> N;  // total number of observations
  int<lower=0> y[N];       // Poisson draws

  // Prior parameters
  real<lower=0> prior_gamma;
  real<lower=0> prior_beta;
}
parameters {
  real<lower=0> lambda[N];
}
transformed parameters {
}
model {
  for (n in 1:N) {
    lambda[n] ~ gamma(prior_gamma, prior_beta);
    y[n] ~ poisson(lambda[n]);
  }
}
