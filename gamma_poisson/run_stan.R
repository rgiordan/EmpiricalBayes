library(rstan)
library(boot) # For inverse logit

project_directory <- file.path(Sys.getenv("GIT_REPO_LOC"), "EmpiricalBayes/gamma_poisson")
data_directory <- file.path(project_directory, "data/")
stan_directory <- file.path(project_directory, "stan/")

set.seed(42)

# Load and cache compiled STAN models
LoadStanModel <- function(stan_model_name) {
  model_file <- file.path(stan_directory, paste(stan_model_name, "stan", sep="."))
  model_file_rdata <- file.path(stan_directory, paste(stan_model_name, "Rdata", sep="."))
  if (file.exists(model_file_rdata)) {
    print("Loading pre-compiled Stan model.")
    stan_model_env <- environment()
    load(model_file_rdata, env=stan_model_env)
    return(stan_model_env$model)
  } else {
    print("Compiling Stan model.")
    model_file <- file.path(stan_directory, paste(stan_model_name, "stan", sep="."))
    model <- stan_model(model_file)
    save(model, file=model_file_rdata)
    return(model)
  }
}



#############################
# Simualate some data

n_obs <- 100

set.seed(42)
true_params <- list()
true_params$prior_gamma <- 30;
true_params$prior_beta <- 3;

true_params$lambda <- rep(NaN, n_obs)
true_params$lambda <- rgamma(n_obs, true_params$prior_gamma, true_params$prior_beta)
y <- rpois(n_obs, lambda=true_params$lambda)


##################################
# Get an empirical Bayes prior

DecodeTheta <- function(theta) {
  # Note: the negative binomial in R is opposite wikipedia notation
  # Here, p will follow R's convention so that p is the probability of
  # the non-terminating probability.

  r <- exp(theta[1])
  p <- inv.logit(theta[2])
  gamma <- r
  beta <- p / (1 - p)
  list(r=r, p=p, gamma=gamma, beta=beta)
}

NegBinLogLik <- function(theta) {
  par <- DecodeTheta(theta)
  log_lik <- sum(dnbinom(y, par$r, par$p, log=TRUE)) 
  cat(log_lik, "\n")
  return(log_lik)
}

prior_mle_optim <- optim(c(0, 0), NegBinLogLik, control=list(fnscale=-1))
prior_mle <- DecodeTheta(prior_mle_optim$par)


gamma_est <- prior_mle$gamma
beta_est <- prior_mle$beta


######################################
# STAN

free_model <- LoadStanModel("gamma_poisson_free_prior")
fixed_model <- LoadStanModel("gamma_poisson_fixed_prior")

# Stan data.
stan_dat <- list(N = length(y),
                 y = y,
                 prior_gamma_mean=gamma_est,
                 prior_gamma_var=20,
                 prior_beta_mean=beta_est,
                 prior_beta_var=20,
                 prior_gamma = gamma_est,
                 prior_beta = beta_est)

# Some knobs we can tweak.  Note that we need many iterations to accurately assess
# the prior sensitivity in the MCMC noise.
chains <- 1
iters <- 50000
seed <- 42

# Draw the draws and save.
stan_draws_file <- file.path(data_directory, "gamma_poisson_mcmc_draws.Rdata")
free_stan_sim <- sampling(free_model, data = stan_dat, seed = seed, chains = chains, iter = iters)
fixed_stan_sim <- sampling(fixed_model, data = stan_dat, seed = seed, chains = chains, iter = iters)


###############################
# Condition to get moments

lambda_fixed <- extract(fixed_stan_sim)
lambda_free <- extract(free_stan_sim)
gamma_draws <- lambda_free$prior_gamma
beta_draws <- lambda_free$prior_beta

moment_list <- list()
for (y_obs in unique(y)) {
  lambda_mean_draws <- (gamma_draws + y_obs) / (beta_draws + 1)
  lambda_var_draws <- (gamma_draws + y_obs) / ((beta_draws + 1)^2)
  lambda_var <- mean(lambda_var_draws) + var(lambda_mean_draws)
  moment_list[[length(moment_list) + 1]] <-
    data.frame(y=y_obs, lambda_mean=mean(lambda_mean_draws), lambda_var=lambda_var)
}
free_moment_df <- do.call(rbind, moment_list)


##################################
# Save

save(free_stan_sim, fixed_stan_sim, stan_dat, true_params, free_moment_df, file=stan_draws_file)


