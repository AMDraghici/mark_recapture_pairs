## Load packages
library(parallel)
library(dclone)
library(tidyverse)
library(readxl)
library(lubridate)
library(rjags)
library(nimble)
library(coda)
library(ggmcmc)

## MCMC parameters
niter <- 1e4
nburnin <- 5e3
nchains <- 5
nthin <- 10

## Load scripts
`%+%` <- function(a, b) paste0(a, b)
src_dir <- getwd()#"/home/sbonner/Students/Statistics/A_Draghici/Research/mark_recapture_pair_swap"
source(file.path(src_dir,"Scripts","pair_swap_mod_nimble.R"))
source(file.path(src_dir,"Scripts","jolly_seber_mod_nimble.R"))
source(file.path(src_dir,"Scripts","00_fn_sim_pair_data.R"))

# TESTING SIMULATED DATA METHOD -------------------------------------------------------------------------------------------------
# Set number of occasions and animals
k = 15
n = 100

# Seeds for Testing
set.seed(1)
# set.seed(1e5)

# Parameter Grid 
param_list <- list(
  n            = n, # Number of Animals
  k            = k, # Occasions
  lf           = 15, # Data Augmentation for Females (M_F)
  lm           = 15, # Data Augmentation for Males (M_M)
  prop.female  = 0.5, # Proportion of simulated individuals to be female
  delta        = rep(0.8, k), # Probability that mating is attempted
  phi.f        = rep(0.9, k), # Marginal Prob of Female Survival
  phi.m        = rep(0.9, k), # Marginal Prob of Male Survival
  gam          = rep(0.0, k), # Correlation in Survival Prob of Mates
  p.f          = rep(0.8, k), # Marginal Prob of Female Recapture
  p.m          = rep(0.8, k), # Marginal Prob of Male Recapture
  rho          = rep(0, k), # Correlation in male survival rates
  betas        = list(beta0 = 90, beta1 = 0), # inv.logit(Beta0 + Beta1 * hij) = Prob of reforming a pair from t-1 after hij times together
  rand_init    = F, # Randomize Initial Entry (just leave as F)
  init         = sample(1, n, TRUE), # Initial Entry into population for individual n
  show_unmated = T, # Include unmated observations in attempt to mate step
  data_aug     = T  # Add individuals to data augmentation
)

# Generate One set of Data
ps_data <- sim_dat(param_list) # pair-swap data
check_sim_output(ps_data)
js_data <- format_to_js(ps_data) 

## Compile model PS------------------------------------------------------------------------
nimble_params <- c("PF","PM","rho","PhiF","PhiM","gamma","beta0", "beta1",
                   "delta","eps","gl","gu","ru","rl","NF","NM")

fit <- run_pair_swap_nimble_parallel(data = ps_data, 
                                     params = nimble_params,
                                     niter = niter, 
                                     nthin = nthin, 
                                     nburnin = nburnin,
                                     ncores = nchains)

samples <- fit$samples
inits <- fit$inits
seeds <- fit$seeds

## Compile model JS-----------------------------------------------------------------------

js_params <- c("PF","PM", "PhiF","PhiM", "xi", "NF", "NM")

fit_js <- run_js_nimble_parallel(data = js_data, 
                              params = js_params,
                              niter = niter, 
                              nthin = nthin, 
                              nburnin = nburnin,
                              ncores = nchains)

samples_js <- fit_js$samples
inits_js <- fit_js$inits
seeds_js <- fit_js$seed


## Summary
summ <- summary(samples_js)
round(cbind(summ[[1]][,"Mean"],summ[[2]][,c("2.5%","97.5%")])[1:6,],3)

## Convergence diagnostics
gelman.diag(samples[,c("NF","NM","PhiF","PhiM","PF","PM", "rho","beta0","delta")])

## Effective sample size
ess <- round(effectiveSize(samples)/(nchains * niter),2)
ess
## Traceplots.
ggs_traceplot(ggs(samples,c("rho")))  + 
  geom_hline(yintercept = param_list$rho[1], col = "red", size = 1.5)
ggs_traceplot(ggs(samples,c("gamma"))) + 
  geom_hline(yintercept = param_list$gam[1], col = "red", size = 1.5)

ggs_traceplot(ggs(samples,c("P"))) # + 
# geom_hline(yintercept = 0.7, col = "red", size = 0.5)
ggs_traceplot(ggs(samples,c("delta"))) + 
  geom_hline(yintercept = param_list$delta[1], col = "red", size = 1.5)
ggs_traceplot(ggs(samples,c("beta")))

ggs_traceplot(ggs(samples,c("N")))

