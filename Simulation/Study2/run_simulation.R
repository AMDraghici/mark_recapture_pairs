## Load Custom Scripts ---------------------------------------------------------------------------------------------
`%+%`      <- function(a, b) paste0(a, b)
src_dir    <-"/home/mdraghic/projects/def-sbonner/mdraghic/mark_recapture_pair_swap/"
out_dir    <- src_dir %+% "/Simulation/Study2/Output/"
source(file.path(src_dir, "Scripts", "fn_generic.R"))
source(file.path(src_dir, "Scripts", "fn_sim_pair_data.R"))
source(file.path(src_dir, "Scripts", "fn_correlation_estimators.R"))

# Load packages ---------------------------------------------------------------------------------------------------
libs <- c("tidyverse","RMark", "nimble","lubridate", "parallel", "coda")
load_packages(libs, FALSE)

# Load Scripts with Dependencies ----------------------------------------------------------------------------------
source(file.path(src_dir, "Scripts", "fn_cormack_jolly_seber_mod_nimble.R"))
source(file.path(src_dir, "Scripts", "fn_pair_swap_mod_nimble.R"))
#------------------------------------------------------------------------------------------------------------------

# Read command line arguments--------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
cat("Running iteration #" %+% iter %+% "...", "\n")
#------------------------------------------------------------------------------------------------------------------

# Simulate Data ---------------------------------------------------------------------------------------------------
x <- Sys.time()

# True Parameter Settings
PM                 <- 0.7
PF                 <- 0.7
PhiF               <- 0.72
PhiM               <- 0.75
gam_true           <- 0.8
rho_true           <- 0.43
n_pop              <- 312
k                  <- 30
ru_true            <- compute_jbin_param_cjs(PF,PM)$cor_lower_bound
rl_true            <- compute_jbin_param_cjs(PF,PM)$cor_upper_bound
gu_true            <- compute_jbin_param_cjs(PhiF,PhiM)$cor_upper_bound
gl_true            <- compute_jbin_param_cjs(PhiF,PhiM)$cor_upper_bound
true_params        <- c(PF,PM, PhiF, PhiM, gam_true, gl_true, gu_true, rho_true, rl_true, ru_true)
param_names        <- c("PF","PM", "PhiF", "PhiM", "gamma", "gl", "gu", "rho", "rl", "ru")
true_param_df      <- data.frame(Truth = true_params, Parameter = param_names)


# Store Random Seed for Reproducibility 
random_seed <- .Random.seed

# Parameter Grid 
param_list <- list(
  n            = n_pop,              # Number of Animals
  k            = k,                  # Occasions
  prop.female  = 0.5,                # Proportion of simulated individuals to be female
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
  init         = NULL,               # Initial Entry into population for individual n
  show_unmated = T                   # Include unmated observations in attempt to mate step 
)

# Generate One set of pair-swap data
ps_data <- sim_dat(param_list) 
# Convert to standard CJS data 
cjs_data <- format_to_cjs(ps_data)

y <- Sys.time()
difftime(y,x,units = "mins")
# ---------------------------------------------------------------------------------------------------------------

# Run CJS Model MARK --------------------------------------------------------------------------------------------
x <- Sys.time()
cjs_out <- run_cjs_model_mark(cjs_data = cjs_data) %>% 
           left_join(true_param_df, by = "Parameter") %>% 
           mutate(Bias    = Truth - Est,
                  In95    = 1*(Truth <= UB & Truth >= LB),
                  Range95 = UB - LB,
                  iter    = iter) 


phim_mark <- cjs_out %>% filter(Parameter == "PhiM") %>% pull(Est)
phif_mark <- cjs_out %>% filter(Parameter == "PhiF") %>% pull(Est)
pm_mark   <- cjs_out %>% filter(Parameter == "PM") %>% pull(Est)
pf_mark   <- cjs_out %>% filter(Parameter == "PF") %>% pull(Est)

y <- Sys.time()
difftime(y,x,units = "mins")
#---------------------------------------------------------------------------------------------------------------

# Compute Recapture Correlation Estimate------------------------------------------------------------------------
x <- Sys.time()
rho <- compute_recapture_correlation(ps_data = ps_data, 
                                     PF      = pf_mark,
                                     PM      = pm_mark)
names(rho) <- "Est"

# Bootstrap To Estimate SE 
rho_bs <- compute_bootstrap_estimates_recapture_correlation(ps_data = ps_data,
                                                            iter    = 10000,
                                                            PF      = pf_mark,
                                                            PM      = pm_mark)

# Collect Results
se_rho <- sd(rho_bs)
names(se_rho) <- "SE"
mean_rho_bs <- mean(rho_bs)
names(mean_rho_bs) <- "Est_Btstrp"
quantiles_rho <- quantile(rho_bs, c(0.025, 0.25, 0.5, 0.75, 0.975))
summ_rho <- c(rho, mean_rho_bs, se_rho, quantiles_rho)
y <- Sys.time()
difftime(y,x,units = "mins")
#-------------------------------------------------------------------------------------------------------------

# Compute Survival Correlation Estimate-----------------------------------------------------------------------
x <- Sys.time()
gamma <- compute_survival_correlation(ps_data = ps_data,
                                      PFM     = compute_jbin_cjs(prob.f = pf_mark,
                                                                 prob.m = pm_mark,
                                                                 cor    = rho)$prob.mf,
                                      PhiF    = phif_mark,
                                      PhiM    = phim_mark)
names(gamma) <- "Est"

# Bootstrap To Estimate SE 
gamma_bs <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                             iter                  = 10000,
                                                             recapture_correlation = rho,
                                                             PF                    = pf_mark,
                                                             PM                    = pm_mark,
                                                             PhiF                  = phif_mark,
                                                             PhiM                  = phim_mark)

# Collect Results
se_gamma <- sd(gamma_bs)
names(se_gamma) <- "SE"
mean_gamma_bs <- mean(gamma_bs)
names(mean_gamma_bs) <- "Est_Btstrp"
quantiles_gamma <- quantile(gamma_bs, c(0.025, 0.25, 0.5, 0.75, 0.975))
summ_gamma <- c(gamma,mean_gamma_bs, se_gamma, quantiles_gamma)
y <- Sys.time()
difftime(y,x,units = "mins")
#-------------------------------------------------------------------------------------------------------------

# Return Results----------------------------------------------------------------------------------------------

# Gather Correlation Results
summ_corr <- as.data.frame(rbind(summ_rho, summ_gamma))
summ_corr$Parameter <- c("rho","gamma")
summ_corr <- summ_corr[,c("Parameter", "Est", "Est_Btstrp", "SE", "2.5%", "25%", "50%","75%","97.5%")]
rownames(summ_corr) <- NULL

# Add Study Statistics 
summ_corr <- summ_corr %>% 
             left_join(true_param_df, by = "Parameter") %>% 
             mutate(Bias    = Truth - Est,
                    In95    = 1*(Truth <= `97.5%` & Truth >= `2.5%`),
                    In50    = 1*(Truth <= `75%`  & Truth  >= `25%`),
                    Range95 = `97.5%` - `2.5%`,
                    Range50 = `75%` - `25%`,
                    iter    = iter) 

#-------------------------------------------------------------------------------------------------------------

# Run Nimble Model with Imputed Pairs ------------------------------------------------------------------------

## MCMC parameters
niter   <- 1e4
nburnin <- niter/4
nchains <- 1
nthin   <- 10
nimble_params <- c("PF","PM", "PhiF","PhiM", "gl","gu","gamma", "ru","rl","rho")

# Run Nimble Pipeline
x <- Sys.time()
fit <- run_pair_swap_nimble_parallel(data    = ps_data, 
                                     params  = nimble_params,
                                     niter   = niter, 
                                     nthin   = nthin, 
                                     nburnin = nburnin,
                                     ncores  = nchains)

# Store Results, add sim study statistics
summ_posterior_ps_imputed <- gather_posterior_summary(fit$samples, 
                                                      nchains = nchains) %>% 
                             left_join(true_param_df, by = "Parameter") %>% 
                             mutate(Bias    = Truth - Mean,
                                    In95    = 1*(Truth <= `97.5%` & Truth >= `2.5%`),
                                    In50    = 1*(Truth <= `75%`  & Truth  >= `25%`),
                                    Range95 = `97.5%` - `2.5%`,
                                    Range50 = `75%` - `25%`,
                                    iter    = iter) 
y <- Sys.time()

difftime(y,x,units = "hours")
# ------------------------------------------------------------------------------------------------------------


# Run Nimble Model with KNOWN Pairs --------------------------------------------------------------------------

## Add known pairs (from simulated data) to mcmc_data object
ps_data2 <- ps_data
ps_data2$apairs_f_imputed <- ps_data$pairs_f
ps_data2$apairs_m_imputed <- ps_data$pairs_m

# Run Nimble Pipeline
x <- Sys.time()
fit2 <- run_pair_swap_nimble_parallel(data    = ps_data2, 
                                      params  = nimble_params,
                                      niter   = niter, 
                                      nthin   = nthin, 
                                      nburnin = nburnin,
                                      ncores  = nchains)

# Store Results, add sim study statistics
summ_posterior_ps_known <- gather_posterior_summary(fit2$samples, 
                                                    nchains = nchains) %>% 
                           left_join(true_param_df, by = "Parameter") %>% 
                           mutate(Bias    = Truth - Mean,
                                  In95    = 1*(Truth <= `97.5%` & Truth >= `2.5%`),
                                  In50    = 1*(Truth <= `75%`  & Truth  >= `25%`),
                                  Range95 = `97.5%` - `2.5%`,
                                  Range50 = `75%` - `25%`,
                                  iter    = iter) 
y <- Sys.time()
difftime(y,x,units = "hours")
# ------------------------------------------------------------------------------------------------------------


# Run Nimble Standard CJS Model ------------------------------------------------------------------------------
nimble_cjs_params <- c("PF","PM","PhiF","PhiM")


# Run Nimble Pipeline
x <- Sys.time()
fit_cjs <- run_cjs_nimble_parallel(data    = cjs_data, 
                                   params  = nimble_cjs_params,
                                   niter   = niter, 
                                   nthin   = nthin, 
                                   nburnin = nburnin,
                                   ncores  = nchains)

# Store Results, add sim study statistics
summ_posterior_cjs <- gather_posterior_summary(fit_cjs$samples, 
                                                    nchains = nchains) %>% 
                      left_join(true_param_df, by = "Parameter") %>% 
                      mutate(Bias    = Truth - Mean,
                             In95    = 1*(Truth <= `97.5%` & Truth >= `2.5%`),
                             In50    = 1*(Truth <= `75%`  & Truth  >= `25%`),
                             Range95 = `97.5%` - `2.5%`,
                             Range50 = `75%` - `25%`,
                             iter    = iter) 
y <- Sys.time()
difftime(y,x,units = "hours")
# ---------------------------------------------------------------------------------------------------------

# Review Results ------------------------------------------------------------------------------------------

results <- list()

# Store Results in List
results$random_seed               <- random_seed
results$mark_out                  <- cjs_out
results$summ_corr                 <- summ_corr
results$summ_posterior_cjs        <- summ_posterior_cjs
results$summ_posterior_ps_imputed <- summ_posterior_ps_imputed
results$summ_posterior_ps_known   <- summ_posterior_ps_known

saveRDS(results, out_dir %+% "results_sim2_" %+% i %+% ".rds")
# ---------------------------------------------------------------------------------------------------------