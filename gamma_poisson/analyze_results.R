library(ggplot2)
library(dplyr)
library(reshape2)
library(rstan)

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

ggplot(lambda_summary) +
    geom_point(aes(x=fixed, y=true_lambda, color="true")) +
    geom_point(aes(x=fixed, y=free, color="free")) +
    geom_abline(aes(intercept=0, slope=1))

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


##############################
# Mean and variance

lambda_stats <-
  select(lambda_grouped, g, draw, method, lambda) %>%
  group_by(g, method) %>%
  summarize(sd=sd(lambda), mean=mean(lambda)) %>%
  inner_join(lambda_true_df, by="g")
  
# Free has more variance
ggplot(dcast(lambda_stats, g ~ method, value.var="sd")) +
  geom_point(aes(x=fixed, y=free)) +
  geom_abline(aes(slope=1, intercept=0))

# Fixed shrinks more
ggplot(dcast(lambda_stats, g + y ~ method, value.var="mean")) +
  geom_point(aes(x=y, y=y - free, color="y - free")) +
  geom_point(aes(x=y, y=y - fixed, color="y - fixed")) +
  geom_hline(aes(yintercept=0))



############################
# Hyperparameter posteriors

prior_df <- data.frame(gamma=stan_results$lambda_free$prior_gamma,
                       beta=stan_results$lambda_free$prior_beta)
beta_est <- stan_results$stan_dat$prior_beta
gamma_est <- stan_results$stan_dat$prior_gamma
stan_results$stan_dat$prior_gamma_mean

gamma_prior_sd <- sqrt(stan_results$stan_dat$prior_gamma_var)
beta_prior_sd <- sqrt(stan_results$stan_dat$prior_beta_var)

# Sanity check
median(lambda_free$gamma_alpha / lambda_free$gamma_beta) - gamma_est
median(lambda_free$beta_alpha / lambda_free$beta_beta) - beta_est

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

M_at <- matrix(NaN, 2, n_g)
M_at[1, ] <- 1 / (1 + beta_est) - 1 / beta_est
M_at[2, ] <- gamma_est / (beta_est ^ 2) - (gamma_est + t(y)) / ((1 + beta_est) ^ 2)

dalpha_dt <- -1 * solve(M_aa, M_at)

# Means within a draw
lambda_moments <-
  ungroup(lambda_fixed_df)  %>%
  group_by(draw) %>%
  summarize(e=mean(lambda), e_log=mean(log(lambda))) %>%
  arrange(draw)

ell_alpha <- matrix(NaN, max(lambda_moments$draw), 2)
ell_alpha[, 1] <- log(beta_est) - digamma(gamma_est) + lambda_moments$e_log
ell_alpha[, 2] <- gamma_est / beta_est - lambda_moments$e
ell_alpha <- ell_alpha * n_g

(gamma_est + dalpha_dt[1, 1]) / (beta_est + dalpha_dt[2, 1]) - gamma_est / beta_est

# Lambda draws as a matrix
lambda_draws <- as.matrix(
  select(lambda_fixed_df, lambda, g, draw) %>%
    dcast(g ~ draw, value.var="lambda") %>%
    select(-g))

lambda_draws <- lambda_draws - rowMeans(lambda_draws)

cov_corr <- lambda_draws %*% ell_alpha %*% dalpha_dt / n_draws
diag(cov_corr)

lambda_means <- rowMeans(lambda_draws)
sample_cov <- lambda_draws %*% t(lambda_draws) / (n_draws - 1) -
              lambda_means %*% t(lambda_means) * n_draws / (n_draws - 1)

lr_cov <- sample_cov + cov_corr

# The change is not that large!
plot(diag(sample_cov), diag(lr_cov)); abline(0, 1)

# Sanity check
max(abs(lambda_means - filter(lambda_stats, method == "fixed")$mean))
max(abs(sqrt(diag(sample_cov)) - filter(lambda_stats, method == "fixed")$sd))
