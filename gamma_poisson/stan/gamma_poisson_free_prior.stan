data {
  int<lower=0> N;  // total number of observations
  int<lower=0> y[N];       // Poisson draws
}
parameters {
  real<lower=0> prior_gamma;
  real<lower=0> prior_beta;
  real<lower=0> lambda[N];
}
transformed parameters {
}
model {
  prior_gamma ~ lognormal(0, 50);
  prior_beta ~ lognormal(0, 50);
  for (n in 1:N) {
    lambda[n] ~ gamma(prior_gamma, prior_beta);
    y[n] ~ poisson(lambda[n]);
  }
}
