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
niter <- 10000
nburnin <- 5000
nchains <- 2
thin <- 1

## Load scripts
`%+%` <- function(a, b) paste0(a, b)
src_dir <- getwd()#"/home/sbonner/Students/Statistics/A_Draghici/Research/mark_recapture_pair_swap"
source(file.path(src_dir,"Scripts","13_pair_swap_mod_nimble.R"))
source(file.path(src_dir,"Scripts","11_jolly_seber_mod_nimble.R"))
source(file.path(src_dir,"Scripts","00_fn_sim_pair_data.R"))


# Functions

# 
# ## Load data
# out_dir <- getwd() %+% "/Simulation/Error_Check_Sim/no_repartner_no_corr_bad/"
# i <- 101
# ps_run_i <- readRDS(out_dir %+% "ps_run_" %+% i %+% ".rds")
# ps_data <- ps_run_i$data
# 
# js_run_i <- readRDS(out_dir %+% "js_run_" %+% i %+% ".rds")
# 


# THINGS TO CHECK
# CAN WE SET DELTA TO A VALUE THAT IS COMPATIBLE WITH THE INITIAL VALUE?
# Does the NaN = 0 trick work? 
# Should we get rid of the mating/not mating mechanism altogether? if so how do we deal with initialization bug...maybe set no partner as a viable option?
#-> that could actually work...try that after the NaN trick


# TESTING SIMULATED DATA METHOD -------------------------------------------------------------------------------------------------
# Set number of occasions and animals
k = 10
n = 50

# Seeds for Testing
set.seed(42)
# set.seed(1e5)

# Parameter Grid 
param_list <- list(
  n            = n, # Number of Animals
  k            = k, # Occasions
  lf           = 20, # Data Augmentation for Females (M_F)
  lm           = 20, # Data Augmentation for Males (M_M)
  prop.female  = 0.5, # Proportion of simulated individuals to be female
  delta        = rep(0.85, k), # Probability that mating is attempted
  phi.f        = rep(0.9, k), # Marginal Prob of Female Survival
  phi.m        = rep(0.9, k), # Marginal Prob of Male Survival
  gam          = rep(0, k), # Correlation in Survival Prob of Mates
  p.f          = rep(0.9, k), # Marginal Prob of Female Recapture
  p.m          = rep(0.9, k), # Marginal Prob of Male Recapture
  rho          = rep(0, k), # Correlation in male survival rates
  betas        = list(beta0 = 90, beta1 = 0), # inv.logit(Beta0 + Beta1 * hij) = Prob of reforming a pair from t-1 after hij times together
  rand_init    = F, # Randomize Initial Entry (just leave as F)
  init         = sample(1, n, TRUE), # Initial Entry into population for individual n
  show_unmated = T, # Include unmated observations in attempt to mate step
  data_aug     = T,  # Add individuals to data augmentation
  full_repartnership = F
)

# Generate One set of Data
ps_data <- sim_dat(param_list) # pair-swap data
js_data <- format_to_js(ps_data) 

## Compile model PS------------------------------------------------------------------------
nimble_params <- c("PF","PM","rho",
                   "rho_raw","gamma_raw",
                   "PhiF","PhiM","gamma","beta0", "beta1",
                   "delta","eps","gl","gu","ru","rl","NF","NM")

CpsMCMC_List <- compile_pair_swap_nimble(ps_data,nimble_params)

# Generate inital values
inits <- lapply(1:nchains, function(i) generate_nimble_init_pairs(ps_data))

## Generate samples
samples <- run_nimble(CpsMCMC_List$CpsMCMC,
                      niter = niter,
                      nburnin = nburnin,
                      thin = thin,
                      inits = inits,
                      nchains = nchains)

## Compile model JS-----------------------------------------------------------------------

js_params <- c("PF","PM", "PhiF","PhiM", "xi", "NF", "NM")

CpsMCMC_JS_List <- compile_jolly_seber_nimble(js_data,js_params)

# Generate inital values
inits_js <- lapply(1:nchains, function(i) generate_init_js(js_data))

## Generate samples
samples_js <- run_nimble(CpsMCMC_JS_List$CjsMCMC,
                         niter = niter,
                         nburnin = nburnin,
                         thin = thin,
                         inits = inits_js,
                         nchains = nchains)


## Summary
summ <- summary(samples)
round(cbind(summ[[1]][,"Mean"],summ[[2]][,c("2.5%","97.5%")])[1:6,],3)

## Convergence diagnostics
gelman.diag(samples[,c("PhiF","PhiM","PF","PM")])

## Effective sample size
ess <- round(effectiveSize(samples)/(nchains * niter),2)
ess
## Save results
# saveRDS(samples,"test_pair_swap_mcmc_bad.rds")

## Traceplots
ggs_traceplot(ggs(samples,c("rho")))  + 
  geom_hline(yintercept = param_list$rho[1], col = "red", size = 1.5)
ggs_traceplot(ggs(samples,c("gamma"))) + 
  geom_hline(yintercept = param_list$gam[1], col = "red", size = 1.5)

ggs_traceplot(ggs(samples_js,c("P")))  + 
  geom_hline(yintercept = 0.9, col = "red", size = 1.5)
ggs_traceplot(ggs(samples,c("delta"))) + 
  geom_hline(yintercept = param_list$delta[1], col = "red", size = 1.5)
ggs_traceplot(ggs(samples,c("beta")))

# Delta is too low
# Arepartner seems okay (although beta is often too high)
# gamma/rho are too high
# need 

