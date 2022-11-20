## Load scripts ------------------------------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
src_dir <- getwd()
source(file.path(src_dir, "Scripts", "fn_generic.R"))
source(file.path(src_dir, "Scripts", "fn_sim_pair_data.R"))
source(file.path(src_dir, "Scripts", "fn_correlation_estimators.R"))

# Load packages
libs <- c("tidyverse","RMark", "nimble", "readxl", "lubridate", "parallel", "coda")
load_packages(libs, FALSE)
source(file.path(src_dir, "Scripts", "fn_cormack_jolly_seber_mod_nimble.R"))
source(file.path(src_dir, "Scripts", "fn_pair_swap_mod_nimble.R"))
source(file.path(src_dir, "Scripts", "fn_process_hduck_data.R"))
# Simulate Data ------------------------------------------------------------------------------------------------
x0 <- Sys.time()
PM       <- 0.7
PF       <- 0.7
PhiF     <- 0.72
PhiM     <- 0.75
gam_true <- 0.8
rho_true <- 0.43
n_pop    <- 312
k        <- 30

set.seed(42)

# Parameter Grid 
param_list <- list(
  n            = n_pop,              # Number of Animals
  k            = k,                  # Occasions
  prop.female  = 0.5,               # Proportion of simulated individuals to be female
  delta        = rep(1, k),          # Probability that mating is attempted
  phi.f        = rep(PhiF, k),       # Marginal Prob of Female Survival
  phi.m        = rep(PhiM, k),       # Marginal Prob of Male Survival
  gam          = rep(gam_true, k),   # Correlation in Survival Prob of Mates
  p.f          = rep(PF, k),         # Marginal Prob of Female Recapture
  p.m          = rep(PM, k),         # Marginal Prob of Male Recapture
  rho          = rep(rho_true, k),   # Correlation in male survival rates
  betas        = list(beta0 = 1000, 
                      beta1 = 1000), # logit pair reform params, beta1 is history coef, assume always repartner
  rand_init    = F,                  # Randomize Initial Entry (just leave as F)
  init         = NULL,       # Initial Entry into population for individual n
  show_unmated = T                   # Include unmated observations in attempt to mate step 
)

# Generate One set of Data
ps_data <- sim_dat(param_list) # pair-swap data
cjs_data <- format_to_cjs(ps_data)
y <- Sys.time()
difftime(y,x0,units = "mins")

# Run CJS Model MARK -----------------------------------------------------------------------------------------
x <- Sys.time()
cjs_out <- run_cjs_model_mark(cjs_data = cjs_data,
                              PhiF     = PhiF,
                              PhiM     = PhiM,
                              PF       = PF,
                              PM       = PM)

pred_probs <- cjs_out$Est
y <- Sys.time()
difftime(y,x,units = "mins")
#-------------------------------------------------------------------------------------------------------------

# Compute Recapture Correlation Estimate----------------------------------------------------------------------
x <- Sys.time()
rho <- compute_recapture_correlation(ps_data = ps_data, 
                                     PF      = pred_probs[4],
                                     PM      = pred_probs[3])
names(rho) <- "Mean"

# Bootstrap To Estimate SE 
rho_bs <- compute_bootstrap_estimates_recapture_correlation(ps_data = ps_data,
                                                            iter    = 10000,
                                                            PF      = pred_probs[4],
                                                            PM      = pred_probs[3])

# Collect Results
se_rho <- sd(rho_bs)
names(se_rho) <- "SE"
quantiles_rho <- quantile(rho_bs, c(0.025, 0.5, 0.75, 0.975))
summ_rho <- c(rho, se_rho, quantiles_rho)
y <- Sys.time()
difftime(y,x,units = "mins")
#-------------------------------------------------------------------------------------------------------------

# Compute Survival Correlation Estimate-----------------------------------------------------------------------
x <- Sys.time()
gamma <- compute_survival_correlation(ps_data = ps_data,
                                      PFM     = compute_jbin_cjs(prob.f = pred_probs[4],
                                                                 prob.m = pred_probs[3],
                                                                 cor    = rho)$prob.mf,
                                      PhiF    = pred_probs[2],
                                      PhiM    = pred_probs[1])
names(gamma) <- "Mean"

# Bootstrap To Estimate SE 
gamma_bs <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                             iter                  = 10000,
                                                             recapture_correlation = rho,
                                                             PF                    = pred_probs[4],
                                                             PM                    = pred_probs[3],
                                                             PhiF                  = pred_probs[2],
                                                             PhiM                  = pred_probs[1])

# Collect Results
se_gamma <- sd(gamma_bs)
names(se_gamma) <- "SE"
quantiles_gamma <- quantile(gamma_bs, c(0.025, 0.5, 0.75, 0.975))

summ_gamma <- c(gamma, se_gamma, quantiles_gamma)
y <- Sys.time()
difftime(y,x,units = "mins")
#-------------------------------------------------------------------------------------------------------------

# Return Results----------------------------------------------------------------------------------------------

# Gather Correlation Results
summ_corr <- as.data.frame(rbind(summ_rho, summ_gamma))
summ_corr$Parameter <- c("rho","gamma")
summ_corr <- summ_corr[,c("Parameter","Mean", "SE", "2.5%","50%","75%","97.5%")]
rownames(summ_corr) <- NULL
#-------------------------------------------------------------------------------------------------------------

# Run Nimble Model with Imputed Fates ------------------------------------------------------------------------

## Compile model PS-------------------------------------------------------------------------------------------

## MCMC parameters
niter <- 1e4
nburnin <- niter/4
nchains <- 5
nthin <- 10
nimble_params <- c("PF","PM","PhiF","PhiM",
                   "gl","gu","gamma",
                   "ru","rl","rho")


# ACCOUNT FOR RECRUITMENT....
x <- Sys.time()
fit <- run_pair_swap_nimble_parallel(data    = ps_data, 
                                     params  = nimble_params,
                                     niter   = niter, 
                                     nthin   = nthin, 
                                     nburnin = nburnin,
                                     ncores  = nchains)

samples <- fit$samples
inits <- fit$inits
seeds <- fit$seed
y <- Sys.time()

difftime(y,x,units = "hours")

## Compile model with TRUTH PS--------------------------------------------------------------------------------

## MCMC parameters
ps_data2 <- ps_data
ps_data2$apairs_f_imputed <- ps_data$pairs_f
ps_data2$apairs_m_imputed <- ps_data$pairs_m

# ACCOUNT FOR RECRUITMENT....
x <- Sys.time()
fit2 <- run_pair_swap_nimble_parallel(data    = ps_data2, 
                                      params  = nimble_params,
                                      niter   = niter, 
                                      nthin   = nthin, 
                                      nburnin = nburnin,
                                      ncores  = nchains)

samples2 <- fit2$samples
inits2 <- fit2$inits
seeds2 <- fit2$seed
y <- Sys.time()

difftime(y,x,units = "hours")


# Run MCMC CJS Model ---------------------------------------------------------------------------------------

## Compile model CJS----------------------------------------------------------------------------------------
nimble_cjs_params <- c("PF","PM","PhiF","PhiM")

# ACCOUNT FOR RECRUITMENT....
x <- Sys.time()

fit_cjs <- run_cjs_nimble_parallel(data    = cjs_data, 
                                   params  = nimble_cjs_params,
                                   niter   = niter, 
                                   nthin   = nthin, 
                                   nburnin = nburnin,
                                   ncores  = nchains)

samples_cjs <- fit_cjs$samples
inits_cjs   <- fit_cjs$inits
seeds_cjs   <- fit_cjs$seed
y <- Sys.time()

difftime(y,x,units = "hours")
# ---------------------------------------------------------------------------------------------------------

# Review Results ------------------------------------------------------------------------------------------

# Base Model Run
cjs_out
summ_corr

# MCMC Summary
summary(samples_cjs)
summary(samples)
summary(samples2)

# Convergence Diagnostics
gelman.diag(samples)
gelman.diag(samples2)
gelman.diag(samples_cjs)

# Plot Densities (Imputed)
ggmcmc::ggs(samples) %>% ggmcmc::ggs_density("P")
ggmcmc::ggs(samples) %>% ggmcmc::ggs_density("gamma")
ggmcmc::ggs(samples) %>% ggmcmc::ggs_density("rho")

# Plot Densities (Known)
ggmcmc::ggs(samples2) %>% ggmcmc::ggs_density("P")
ggmcmc::ggs(samples2) %>% ggmcmc::ggs_density("gamma")
ggmcmc::ggs(samples2) %>% ggmcmc::ggs_density("rho")
# ---------------------------------------------------------------------------------------------------------



