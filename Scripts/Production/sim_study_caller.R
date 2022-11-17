## Load scripts ------------------------------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
src_dir <- getwd()
source(file.path(src_dir, "Scripts", "Production", "fn_generic.R"))
source(file.path(src_dir, "Scripts", "Production", "fn_sim_pair_data.R"))
source(file.path(src_dir, "Scripts", "Production", "fn_correlation_estimators.R"))

# Load packages
libs <- c("tidyverse","RMark")
load_packages(libs, FALSE)

# Simulate Data ------------------------------------------------------------------------------------------------
x0 <- Sys.time()
PM       <- 0.85
PF       <- 0.85
PhiF     <- 0.8
PhiM     <- 0.8
gam_true <- 0.25
rho_true <- 0.4
n_pop    <- 100
k        <- 10

# Parameter Grid 
param_list <- list(
  n            = n_pop,              # Number of Animals
  k            = k,                  # Occasions
  prop.female  = 0.45,               # Proportion of simulated individuals to be female
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
  init         = rep(1,n_pop),       # Initial Entry into population for individual n
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
                                     PF      = pred_probs[3],
                                     PM      = pred_probs[4])
names(rho) <- "Mean"

# Bootstrap To Estimate SE 
rho_bs <- compute_bootstrap_estimates_recapture_correlation(ps_data = ps_data,
                                                            iter    = 10000,
                                                            PF      = pred_probs[3],
                                                            PM      = pred_probs[4])

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
                                      PFM     = compute_jbin_cjs(pred_probs[1], pred_probs[2], rho)$prob.mf,
                                      PhiF    = pred_probs[1],
                                      PhiM    = pred_probs[2])
names(gamma) <- "Mean"

# Bootstrap To Estimate SE 
gamma_bs <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                             iter                  = 10000,
                                                             recapture_correlation = NULL,
                                                             PF                    = pred_probs[3],
                                                             PM                    = pred_probs[4],
                                                             PhiF                  = pred_probs[1],
                                                             PhiM                  = pred_probs[2])

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

results <- list(ps_data   = ps_data,
                cjs_data  = cjs_data,
                cjs_out   = cjs_out,
                summ_corr = summ_corr,
                rho_bs    = unname(rho_bs),
                gamma_bs  = unname(gamma_bs))

saveRDS(results, out_dir %+% "results.rds")

y0 <- Sys.time()
difftime(y0,x0,units = "mins")
#-------------------------------------------------------------------------------------------------------------