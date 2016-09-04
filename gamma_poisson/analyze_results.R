library(ggplot2)
library(dplyr)
library(reshape2)
library(rstan)
library(Matrix)

project_directory <- file.path(Sys.getenv("GIT_REPO_LOC"), "EmpiricalBayes/gamma_poisson")
data_directory <- file.path(project_directory, "data/")
stan_directory <- file.path(project_directory, "stan/")

stan_draws_file <- file.path(data_directory, "gamma_poisson_mcmc_draws.Rdata")
stan_results <- environment()
load(stan_draws_file, env=stan_results)

lambda_fixed <- extract(stan_results$fixed_stan_sim)
lambda_free <- extract(stan_results$free_stan_sim)

true_lambda <- stan_results$true_params$lambda
lambda_true_df <- data.frame(g=1:length(true_lambda), "true_lambda"=true_lambda)
lambda_true_df$y <- stan_results$stan_dat$y
gamma_est <- stan_results$stan_dat$prior_gamma
beta_est <- stan_results$stan_dat$prior_beta
lambda_true_df$fixed_mean <- with(lambda_true_df, (y + gamma_est) / (1 + beta_est))
lambda_true_df$fixed_var <- with(lambda_true_df, (y + gamma_est) / ((1 + beta_est)^2))
lambda_true_df <-
  inner_join(lambda_true_df, rename(free_moment_df, free_mean=lambda_mean, free_var=lambda_var), by="y")

if (FALSE) {
  ggplot(lambda_true_df) +
    geom_point(aes(x=free_mean, y=fixed_mean)) +
    geom_abline(aes(slope=1, intercept=0))

  ggplot(lambda_true_df) +
    geom_point(aes(x=free_var, y=fixed_var)) +
    geom_abline(aes(slope=1, intercept=0))
}

lambda_free_df <-
  melt(lambda_free$lambda) %>% rename(g=Var2, draw=iterations, lambda=value) %>%
  mutate(method="free")
lambda_fixed_df <-
  melt(lambda_fixed$lambda) %>% rename(g=Var2, draw=iterations, lambda=value) %>%
  mutate(method="fixed")
lambda_df <- rbind(lambda_free_df, lambda_fixed_df)

lambda_summary <-
  group_by(lambda_df, method, g) %>%
  summarize(mean=mean(lambda), sd=sd(lambda)) %>%
  dcast(g ~ method, value.var="mean") %>%
  inner_join(lambda_true_df, by="g")

# ggplot(lambda_summary) +
#     geom_point(aes(x=fixed, y=true_lambda, color="true")) +
#     geom_point(aes(x=fixed, y=free, color="free")) +
#     geom_abline(aes(intercept=0, slope=1))


##############################
# Accuracy

lambda_grouped <-
  dcast(lambda_df, g + draw ~ method, value.var="lambda") %>%
  melt(id.vars=c("g", "draw")) %>%
  rename(method=variable, lambda=value) %>%
  inner_join(lambda_true_df, by="g")

rmse <-
  group_by(lambda_grouped, g, method) %>%
  summarize(abs_bias=abs(mean(lambda - true_lambda)),
            sd=sd(lambda - true_lambda),
            rmse=sqrt(mean((lambda - true_lambda) ^ 2))) %>%
  melt(id.vars=c("g", "method")) %>%
  rename(measure=variable)

#dcast(rmse, g ~ measure + method)
group_by(rmse, method, measure) %>%
  summarize(mean=mean(value)) %>%
  ungroup() %>% arrange(measure, method)



############################
# Hyperparameter posteriors

prior_df <- data.frame(gamma=stan_results$lambda_free$prior_gamma,
                       beta=stan_results$lambda_free$prior_beta)
beta_est <- stan_results$stan_dat$prior_beta
gamma_est <- stan_results$stan_dat$prior_gamma
stan_results$stan_dat$prior_gamma_mean

gamma_prior_sd <-
  sqrt(stan_results$stan_dat$gamma_prior_alpha) / stan_results$stan_dat$gamma_prior_beta
beta_prior_sd <-
  sqrt(stan_results$stan_dat$beta_prior_alpha) / stan_results$stan_dat$beta_prior_beta

if (FALSE) {
  ggplot(prior_df) +
    geom_histogram(aes(x=gamma), bins=100) +
    geom_vline(aes(xintercept=gamma_est), lwd=3) +
    geom_vline(aes(xintercept=gamma_est - 2 * gamma_prior_sd), lwd=1) +
    geom_vline(aes(xintercept=gamma_est + 2 * gamma_prior_sd), lwd=1)

  ggplot(prior_df) +
    geom_histogram(aes(x=beta), bins=100) +
    geom_vline(aes(xintercept=beta_est), lwd=3) +
    geom_vline(aes(xintercept=beta_est - 2 * beta_prior_sd), lwd=1) +
    geom_vline(aes(xintercept=beta_est + 2 * beta_prior_sd), lwd=1)

    ggplot(prior_df) +
    geom_histogram(aes(x=beta), bins=100) +
    geom_vline(aes(xintercept=beta_est), lwd=3)
}

########################
# Linear response

# Second derivatives
n_g <- max(lambda_fixed_df$g)
n_draws <- max(lambda_fixed_df$draw)

y <- stan_results$stan_dat$y
gamma_est <- stan_results$stan_dat$prior_gamma
beta_est <- stan_results$stan_dat$prior_beta

# Note that the leading N in M_aa cancels the leading N in ell_alpha
M_aa <- matrix(NaN, 2, 2)
M_aa[1, 1] <- mean(trigamma(gamma_est + y)) - trigamma(gamma_est)
M_aa[2, 1] <- M_aa[1, 2] <- 1 / beta_est - 1 / (1 + beta_est)
M_aa[2, 2] <- (gamma_est + mean(y)) / ((1 + beta_est) ^ 2) - gamma_est / (beta_est ^ 2)
M_aa <- M_aa * n_g

# Add the prior Hessian if you're comparing with the MAP and not the maximum marignal likelihood
M_prior_aa <- matrix(0, 2, 2)
gamma_prior_alpha <- stan_results$stan_dat$gamma_prior_alpha
beta_prior_alpha <- stan_results$stan_dat$beta_prior_alpha
M_prior_aa[1, 1] <- -1 * (gamma_prior_alpha - 1) / (gamma_est ^ 2)
M_prior_aa[2, 2] <- -1 * (beta_prior_alpha - 1) / (beta_est ^ 2)

M_at <- matrix(NaN, 2, n_g)
M_at[1, ] <- 1 / (1 + beta_est) - 1 / beta_est
M_at[2, ] <- gamma_est / (beta_est ^ 2) - (gamma_est + t(y)) / ((1 + beta_est) ^ 2)

# Since we're using MLE estimates rather than MAP use the corresponding objective function.
M_hess <- matrix(NaN, 2, 2)
if (stan_results$map_estimate) {
  M_hess <- M_aa + M_prior_aa
} else {
  M_hess <- M_aa
}


##########################
# Exact version

dmom_dalpha <- matrix(NaN, n_g, 2)

# d / d gamma
dmom_dalpha[, 1] <- 1 / (1 + beta_est)

# d / d beta
dmom_dalpha[, 2] <- -1 * (gamma_est + y) / ((1 + beta_est) ^ 2)

dmom_dt <- Diagonal(x=(gamma_est + y) / ((1 + beta_est) ^ 2))

lr_cov_corr_exact <- -1 * dmom_dalpha %*% solve(M_hess, t(dmom_dalpha))

stopifnot(max(abs(diag(dmom_dt) - lambda_true_df$fixed_var)) < 1e-12)


###########################
# Summaries

lambda_stats <-
  select(lambda_grouped, g, draw, method, lambda) %>%
  group_by(g, method) %>%
  summarize(sd=sd(lambda), mean=mean(lambda))

fixed_sd <-
  filter(lambda_stats, method == "fixed") %>%
  arrange(g) %>%
  `[[`("sd")

lr_var <- diag(lr_cov_corr_exact) + fixed_sd ^ 2
lr_diff_exact <- sqrt(lr_var) - fixed_sd



########################
# Plots

# Compare corrections to actual differences.
lambda_stats_correction <-
  data.frame(g=1:n_g, lr_diff_exact=lr_diff_exact) %>%
  inner_join(lambda_true_df, by="g") %>%
  mutate(mcmc_diff=sqrt(free_var) - sqrt(fixed_var))

ggplot(lambda_stats_correction) +
  geom_point(aes(x=mcmc_diff, y=lr_diff_exact)) +
  geom_abline(aes(slope=1, intercept=0)) +
  expand_limits(x=0, y=0)


if (FALSE) {

# Compare fixed vs free
ggplot(lambda_true_df) +
  geom_point(aes(x=free_mean, y=fixed_mean)) +
  geom_abline(aes(slope=1, intercept=0)) +
  expand_limits(x=0, y=0)

max_sd <- max(lambda_true_df$free_var)
ggplot(lambda_true_df) +
  geom_point(aes(x=sqrt(free_var), y=sqrt(fixed_var), color="MCMC")) +
  geom_abline(aes(slope=1, intercept=0)) +
  expand_limits(x=0, y=0) +
  expand_limits(x=max_sd, y=max_sd)

# Free has more variance
ggplot(lambda_true_df) +
  geom_line(aes(x=y, y=free_var, color="free")) +
  geom_line(aes(x=y, y=fixed_var, color="fixed")) +
  expand_limits(x=0, y=0)

# Free shrinks more (just a little bit though)
ggplot(lambda_true_df) +
  geom_line(aes(x=y, y=free_mean, color="free")) +
  geom_line(aes(x=y, y=fixed_mean, color="fixed")) +
  geom_abline(aes(slope=1, intercept=0)) +
  expand_limits(x=0, y=0)

ggplot(lambda_true_df) +
  geom_point(aes(x=y, y=free_var - fixed_var, color="MCMC")) +
  geom_point(aes(x=y, y=diag(lr_cov_corr_exact), color="LR")) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0))

# # This should be an MCMC estimate of dmom / dalpha
# ggplot(lambda_stats_correction) +
#   geom_point(aes(x=lr_diff, y=lr_diff_exact)) +
#   geom_abline(aes(slope=1, intercept=0)) +
#   expand_limits(x=0, y=0)

# Compare exact vs MCMC estimates
lambda_sd_check  <-
  dcast(lambda_stats, g ~ method, value.var="sd") %>%
  inner_join(lambda_true_df, by="g")

ggplot(lambda_sd_check) +
  geom_point(aes(x=sqrt(free_var), y=free, color="free")) +
  geom_point(aes(x=sqrt(fixed_var), y=fixed, color="fixed")) +
  geom_abline(aes(slope=1, intercept=0)) +
  expand_limits(x=0, y=0)

lambda_mean_check  <-
  mutate(lambda_stats, upper=mean + 2 * sd, lower=mean - 2 * sd) %>%
  melt(id.var=c("g", "method")) %>%
  filter(variable != "sd") %>%
  dcast(g ~ method + variable, value.var="value") %>%
  inner_join(select(lambda_true_df, g, true_lambda, y), by="g")

ggplot(lambda_mean_check) +
  geom_point(aes(x=true_lambda, y=free_mean, color="free")) +
  geom_crossbar(aes(x=true_lambda, y=free_mean, ymin=free_lower, ymax=free_upper, color="free")) +
  geom_point(aes(x=true_lambda, y=fixed_mean, color="fixed")) +
  geom_crossbar(aes(x=true_lambda, y=fixed_mean, ymin=fixed_lower, ymax=fixed_upper, color="fixed")) +
  geom_abline(aes(slope=1, intercept=0)) +
  expand_limits(x=0, y=0)


# Comparisions
lambda_stats_lr <-
  data.frame(g=1:n_g, method="lr", sd=sqrt(diag(lr_cov))) %>%
  rbind(lambda_stats) %>%
  inner_join(lambda_true_df, by="g")

lambda_stats_lr_wide <-
  dcast(lambda_stats_lr, g ~ method, value.var="sd") %>%
  inner_join(lambda_true_df, by="g")


ggplot(dcast(lambda_stats, g ~ method, value.var="sd") %>%
  inner_join(lambda_true_df, by="g")) +
  geom_point(aes(x=y, y=fixed))

# L1 error
with(lambda_stats_lr_wide, sum(abs(lr - free)))
with(lambda_stats_lr_wide, sum(abs(fixed - free)))

ggplot(lambda_stats_lr_wide) +
  geom_point(aes(x=free, y=fixed, color="fixed")) +
  geom_point(aes(x=free, y=lr, color="lr")) +
  geom_abline(aes(intercept=0, slope=1)) +
  ylab("Approximate sd")

ggplot(lambda_stats_lr_wide) +
  geom_point(aes(x=y, y=free, color="free")) +
  geom_point(aes(x=y, y=fixed, color="fixed")) +
  geom_point(aes(x=y, y=lr, color="linear response")) +
  ylab("Approximate sd")


ggplot(lambda_stats_lr_wide) +
  geom_point(aes(x=true_sd, y=free, color="free")) +
  geom_point(aes(x=true_sd, y=fixed, color="fixed")) +
  geom_point(aes(x=true_sd, y=lr, color="linear response")) +
  geom_abline(aes(intercept=0, slope=1)) +
  ylab("Approximate sd")

ggplot(lambda_stats_lr_wide) +
  # geom_point(aes(x=y, y=free, color="free")) +
  geom_point(aes(x=y, y=fixed, color="fixed")) +
  ylab("Approximate sd")


# Sanity check
max(abs(lambda_means - filter(lambda_stats, method == "fixed")$mean))
max(abs(sqrt(diag(sample_cov)) - filter(lambda_stats, method == "fixed")$sd))

}


