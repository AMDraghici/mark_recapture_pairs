## Load scripts ------------------------------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
src_dir <- getwd()
source(file.path(src_dir, "Scripts", "Production", "fn_generic.R"))
source(file.path(src_dir, "Scripts", "Production", "fn_sim_pair_data.R"))
source(file.path(src_dir, "Scripts", "Production", "fn_correlation_estimators.R"))

# Load packages
libs <- c("tidyverse","RMark", "nimble", "readxl", "lubridate", "parallel", "coda")
load_packages(libs, FALSE)

source(file.path(src_dir, "Scripts", "Production", "fn_cormack_jolly_seber_mod_nimble.R"))
source(file.path(src_dir, "Scripts", "Production", "fn_pair_swap_mod_nimble.R"))
source(file.path(src_dir, "Scripts", "Production", "fn_process_hduck_data.R"))


# Pull/Process Data ---------------------------------------------------------------------------------------------------------------------------
dat_dir <- src_dir %+% "/Data/RE__Harlequin_duck_data/"
cap.data <- gather_hq_data(dat_dir) %>% 
  build_cr_df() %>% 
  populate_missing_mate_data() %>%
  populate_missing_mate_data() %>% 
  add_implied_states() %>%
  add_last_capture() %>% 
  clean_filtered() 

drop_id_yr_filter <- cap.data %>%
  group_by(animal_id) %>%
  summarize(num = sum(recapture_individual)) %>% 
  filter(num < 1) %>%
  pull(animal_id)

cap.data <- cap.data %>%
  filter(!(animal_id %in% drop_id_yr_filter)) %>%
  assign_ids_bysex()

ps_data  <- build_nimble_data(cap.data, F)
cjs_data <- format_to_cjs(ps_data)

# Run CJS Model MARK -----------------------------------------------------------------------------------------
x <- Sys.time()
cjs_out <- run_cjs_model_mark(cjs_data = cjs_data,
                              PhiF     = 0,
                              PhiM     = 0,
                              PF       = 0,
                              PM       = 0)

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
                                      PFM     = compute_jbin_cjs(pred_probs[3], pred_probs[4], rho)$prob.mf,
                                      PhiF    = pred_probs[1],
                                      PhiM    = pred_probs[2])
names(gamma) <- "Mean"

# Bootstrap To Estimate SE 
gamma_bs <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                             iter                  = 10000,
                                                             recapture_correlation = rho,
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

# Run MCMC CJS Model ---------------------------------------------------------------------------------------

## Compile model CJS----------------------------------------------------------------------------------------
nimble_params <- c("PF","PM","PhiF","PhiM")

# ACCOUNT FOR RECRUITMENT....
x <- Sys.time()

fit_cjs <- run_cjs_nimble_parallel(data    = cjs_data, 
                                   params  = nimble_params,
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

