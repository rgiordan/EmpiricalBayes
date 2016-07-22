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



#########
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

# There's actually quite a lot of dispersion even with a highly informative prior
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


