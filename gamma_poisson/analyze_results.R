library(ggplot2)
library(dplyr)
library(reshape2)
library(rstan)
library(Matrix)
library(boot) # For inverse logit

project_directory <- file.path(Sys.getenv("GIT_REPO_LOC"), "EmpiricalBayes/gamma_poisson")
data_directory <- file.path(project_directory, "data/")
stan_directory <- file.path(project_directory, "stan/")

stan_draws_file <- file.path(data_directory, "gamma_poisson_mcmc_draws.Rdata")
stan_results <- environment()
load(stan_draws_file, env=stan_results)