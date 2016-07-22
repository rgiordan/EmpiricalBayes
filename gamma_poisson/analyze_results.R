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

group_by(rmse, method, measure) %>% summarize(mean=mean(value))
dcast(rmse, g ~ measure + method)