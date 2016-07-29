# Debugging code assuming you have already run run_stan.R in the current environment.



#####################################
# Check posterior

gamma_est <- stan_dat$prior_gamma
beta_est <- stan_dat$prior_beta


lambda_fixed <- extract(fixed_stan_sim)

obs <- 4
lambda <- lambda_fixed$lambda[, obs]
true_sd <- sqrt(y[obs] + gamma_est) / (1 + beta_est)
true_sd_sd <- sqrt(SampleSDVariance(length(lambda), gamma_est, beta_est))
draws <- rgamma(length(lambda), shape=y[obs] + gamma_est, rate=1 + beta_est)
# plot(sort(lambda), sort(draws)); abline(0, 1)
mean(lambda)
mean(draws)
sd(lambda)
sd(draws)
sd(draws2)


cdf_draws_list <- list()
for (obs in 1:ncol(lambda_fixed$lambda)) {
  lambda <- lambda_fixed$lambda[, obs]
  cdf_draws_list[[obs]] <- pgamma(lambda, shape=y[obs] + gamma_est, rate=1 + beta_est)
}
cdf_draws <- do.call(c, cdf_draws_list)
hist(cdf_draws, 200)

sd_df_list <- list()
for (obs in 1:ncol(lambda_fixed$lambda)) {
  lambda <- lambda_fixed$lambda[, obs]
  alpha <- y[obs] + gamma_est
  beta <- 1 + beta_est
  true_sd <- sqrt(alpha) / beta
  true_sd_sd <- sqrt(SampleSDVariance(length(lambda), alpha, beta))
  sd_df_list[[obs]] <- data.frame(obs=obs, true_sd=true_sd, true_sd_sd=true_sd_sd, mcmc_sd=sd(lambda))
}
sd_df <- do.call(rbind, sd_df_list)





GammaMoment <- function(k, shape, rate) {
  exp(lgamma(k + shape) - lgamma(shape) - k * log(rate))
}

draws <- rgamma(10000, shape=gamma_est, rate=beta_est)
GammaMoment(1, gamma_est, beta_est)
mean(draws)

GammaMoment(4, gamma_est, beta_est)
mean(draws ^ 4)

SampleVarianceVariance <- function(n, shape, rate) {
  m1 <- GammaMoment(1, shape, rate)
  m2 <- GammaMoment(2, shape, rate)
  m3 <- GammaMoment(3, shape, rate)
  m4 <- GammaMoment(4, shape, rate)
  
  mu2 <- m2 - (m1^2)
  mu4 <- m4 - 4 * m3 * m1 + 6 * m2 * (m1^2) - 3 * (m1^4)

  var_s2 <- mu4 * ((n - 1)^2) / (n^3) - (n - 1) * (n - 3) * (mu2 ^ 2) / (n^3)  
  return(var_s2)
}

SampleSDVariance <- function(n, shape, rate) {
  s2_var <- SampleVarianceVariance(n, shape, rate)
  # Delta method:
  s2_mean <- GammaMoment(2, shape, rate) - (GammaMoment(1, shape, rate)^2)
  return(s2_var / (4 * s2_mean))
}


vars_list <- list()
n_sim <- 5000
for (sim in 1:10000) {
  draws <- rgamma(n_sim, shape=gamma_est, rate=beta_est)
  vars_list[[sim]] <- var(draws)  
}
vars <- unlist(vars_list)


sd(sqrt(vars))
sqrt(SampleSDVariance(n_sim, gamma_est, beta_est))


################################
# Check derivatives
epsilon <- 1
obs <- 1

library(numDeriv)

t_vec <- rep(0, n_obs)
NegBinLogLikArgs <- function(theta) {
  gamma <- theta[1]
  beta <- theta[2]
  t_vec <- theta[3:length(theta)]
  log_lik <-
    sum(gamma * log(beta - t_vec) - lgamma(gamma)) -
    sum((gamma + y) * log(beta + 1 - t_vec)) +
    sum(lgamma(gamma + y))
  return(log_lik)
}

NegBinLogLikManual <- function(theta) {
  par <- DecodeTheta(theta)
  return(NegBinLogLikArgs(c(par$gamma, par$beta, rep(0, n_obs))))
}

gamma_est <- prior_mle$gamma
beta_est <- prior_mle$beta
M_aa[1, 1] <- sum(trigamma(gamma_est + y) - trigamma(gamma_est))
M_aa[2, 1] <- M_aa[1, 2] <- n_obs * (1 / beta_est - 1 / (1 + beta_est))
M_aa[2, 2] <- sum(gamma_est / (beta_est ^ 2) - (gamma_est + y) / ((1 + beta_est) ^ 2))

auto_hess <- hessian(NegBinLogLikArgs, c(gamma_est, beta_est, rep(0, n_obs)))
auto_hess[1:2, 1:2]
M_aa

auto_hess[1:2, 3:(2 + n_obs)][1, ] - M_at[1,]
auto_hess[1:2, 3:(2 + n_obs)][2, ] - M_at[2,] 

# Look at manual derivatives.

vals <- list()
for (epsilon in seq(0, 1, length.out=100)) {
  t_vec <- rep(0, n_obs)
  t_vec[obs] <- epsilon
  manual_prior_mle_optim <- optim(c(0, 0), NegBinLogLikManual, control=list(fnscale=-1))
  manual_prior_mle <- DecodeTheta(manual_prior_mle_optim$par)
  vals[[length(vals) + 1]] <- c(epsilon, manual_prior_mle$gamma, manual_prior_mle$beta)
}
vals <- do.call(rbind, vals)

plot(vals[, 1], vals[, 2])
plot(vals[, 1], vals[, 3])

ind_diff <- 10
(vals[ind_diff, 2] - vals[1, 2]) / (vals[ind_diff, 1] - vals[1, 1])
(vals[ind_diff, 3] - vals[1, 3]) / (vals[ind_diff, 1] - vals[1, 1])

dalpha_dt[, obs] * epsilon
