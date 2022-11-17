# Load packages
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
niter <- 1e5
nburnin <- niter/4
nchains <- 1
nthin <- 10

## Load scripts
`%+%` <- function(a, b) paste0(a, b)
src_dir <- getwd()#"/home/sbonner/Students/Statistics/A_Draghici/Research/mark_recapture_pair_swap"
source(file.path(src_dir,"Scripts","cormack_jolly_seber_mod_nimble.R"))
source(file.path(src_dir,"Scripts","fn_sim_pair_data3.R"))
source(file.path(src_dir, "Scripts", "fn_correlation_estimators.R"))

# TESTING SIMULATED DATA METHOD -------------------------------------------------------------------------------------------------
PF <- 0.45
PM <- 0.45
PhiF <- 0.8
PhiM <- 0.8
gam_true <- 0.25
rho_true <- 0.4
PFM <- compute_jbin_cjs(PF,PM,rho_true)$prob.mf 
PhiMF <- compute_jbin_cjs(PhiF,PhiM,gam_true)$prob.mf
prob_prod <- PFM * PhiMF
n_pop <- 500
k <- 30

# Parameter Grid 
param_list <- list(
  n            = n_pop, # Number of Animals
  k            = k, # Occasions
  prop.female  = 0.45, # Proportion of simulated individuals to be female
  delta        = rep(1, k), # Probability that mating is attempted
  phi.f        = rep(PhiF, k), # Marginal Prob of Female Survival
  phi.m        = rep(PhiM, k), # Marginal Prob of Male Survival
  gam          = rep(gam_true, k), # Correlation in Survival Prob of Mates
  p.f          = rep(PF, k), # Marginal Prob of Female Recapture
  p.m          = rep(PM, k), # Marginal Prob of Male Recapture
  rho          = rep(rho_true, k), # Correlation in male survival rates
  betas        = list(beta0 = 1000, beta1 = 1000), # inv.logit(Beta0 + Beta1 * hij) = Prob of reforming a pair from t-1 after hij times together
  rand_init    = F, # Randomize Initial Entry (just leave as F)
  init         = sample(1, n_pop, TRUE), # Initial Entry into population for individual n
  show_unmated = T # Include unmated observations in attempt to mate step
)

# Generate One set of Data
ps_data <- sim_dat(param_list) # pair-swap data
cjs_data <- format_to_cjs(ps_data)

# Run model CJS---------------------------------------------------------------------------------------------
cjs_params <- c("PF","PM", "PhiF","PhiM")
x <- Sys.time()
fit_cjs <- run_cjs_nimble_parallel(data    = cjs_data,
                                   params  = cjs_params,
                                   niter   = niter,
                                   nthin   = nthin,
                                   nburnin = nburnin,
                                   ncores  = nchains)

samples_cjs <- fit_cjs$samples
inits_cjs <- fit_cjs$inits
seeds_cjs <- fit_cjs$seed
y <- Sys.time()

difftime(y,x,units = "mins")

pred_probs <- summary(samples_cjs)$statistics[,1]

# Compute Recapture Correlation Estimate---------------------------------------------------------------------
x <- Sys.time()
rho <- compute_recapture_correlation(ps_data = ps_data, 
                                     PF      = pred_probs[1],
                                     PM      = pred_probs[2])
names(rho) <- "Mean"

rho_bs <- compute_bootstrap_estimates_recapture_correlation(ps_data = ps_data,
                                                            iter    = 10000,
                                                            PF      = pred_probs[1],
                                                            PM      = pred_probs[2])

se_rho <- sd(rho_bs)
names(se_rho) <- "SE"
quantiles_rho <- quantile(rho_bs, c(0.025, 0.5, 0.75, 0.975))

summ_rho <- c(rho, se_rho, quantiles_rho)
y <- Sys.time()
difftime(y,x,units = "mins")

# Compute Survival Correlation Estimate---------------------------------------------------------------------
x <- Sys.time()
gamma <- compute_survival_correlation(ps_data = ps_data,
                                      PFM     = compute_jbin_cjs(pred_probs[1], pred_probs[2], rho)$prob.mf,
                                      PhiF    = pred_probs[3],
                                      PhiM    = pred_probs[4])
names(gamma) <- "Mean"

gamma_bs <- compute_bootstrap_estimates_survival_correlation(ps_data = ps_data,
                                                             iter    = 10000,
                                                             recapture_correlation = 0,
                                                             PF      = pred_probs[1],
                                                             PM      = pred_probs[2],
                                                             PhiF    = pred_probs[3],
                                                             PhiM    = pred_probs[4])
se_gamma <- sd(gamma_bs)
names(se_gamma) <- "SE"
quantiles_gamma <- quantile(gamma_bs, c(0.025, 0.5, 0.75, 0.975))

summ_gamma <- c(gamma, se_gamma, quantiles_gamma)
y <- Sys.time()
difftime(y,x,units = "mins")
# Return Results----------------------------------------------------------------------------------------------
summ_corr <- as.data.frame(rbind(summ_rho, summ_gamma))
summ_corr$Parameter <- c("rho","gamma")
summ_corr <- summ_corr[,c("Parameter","Mean", "SE", "2.5%","50%","75%","97.5%")]
rownames(summ_corr) <- NULL

results <- list(fit_cjs   = fit_cjs,
                summ_corr = summ_corr,
                rho_bs    = unname(rho_bs),
                gamma_bs  = unname(gamma_bs))
