## Load scripts ------------------------------------------------------------------------------------------------------------------
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


# Pull/Process Data -------------------------------------------------------------------------------------------------------------
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

ps_data  <- build_nimble_data(cap.data)
cjs_data <- format_to_cjs(ps_data)

# Run CJS Model MARK -----------------------------------------------------------------------------------------------------------
x <- Sys.time()
cjs_out <- run_cjs_model_mark(cjs_data = cjs_data)
phim_mark <- cjs_out %>% filter(Parameter == "PhiM") %>% pull(Est)
phif_mark <- cjs_out %>% filter(Parameter == "PhiF") %>% pull(Est)
pm_mark   <- cjs_out %>% filter(Parameter == "PM") %>% pull(Est)
pf_mark   <- cjs_out %>% filter(Parameter == "PF") %>% pull(Est)
y <- Sys.time()
difftime(y,x,units = "mins")
#-------------------------------------------------------------------------------------------------------------------------------

# Compute Recapture Correlation Estimate----------------------------------------------------------------------------------------
x <- Sys.time()
rho <- compute_recapture_correlation(ps_data = ps_data, 
                                     PF      = pf_mark,
                                     PM      = pm_mark)
names(rho) <- "Est"

# Bootstrap To Estimate SE 
# Non-Parametric Bootstrap To Estimate SE 
rho_bs_np <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                               iter       = 1000,
                                                               PF         = pf_mark,
                                                               PM         = pm_mark,
                                                               rho        = rho,
                                                               parametric = F)
# Collect Results
mean_bstrp_rho_np        <- mean(rho_bs_np)
names(mean_bstrp_rho_np) <- "Est_Btstrp"
se_bstrp_rho_np          <- sd(rho_bs_np)
names(se_bstrp_rho_np)   <- "SE"
quantiles_rho_np         <- quantile(rho_bs_np, c(0.025, 0.25, 0.5, 0.75, 0.975))
status_np                <- 0
names(status_np)         <- "Parametric"
summ_rho_np              <- c(rho, mean_bstrp_rho_np, se_bstrp_rho_np, quantiles_rho_np, status_np)

rho_bs_sp <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                               iter       = 1000,
                                                               PF         = pf_mark,
                                                               PM         = pm_mark,
                                                               rho        = rho,
                                                               parametric = T)
# Collect Results
mean_bstrp_rho_sp         <- mean(rho_bs_sp)
names(mean_bstrp_rho_sp)  <- "Est_Btstrp"
se_bstrp_rho_sp           <- sd(rho_bs_sp)
names(se_bstrp_rho_sp)    <- "SE"
quantiles_rho_sp         <- quantile(rho_bs_sp, c(0.025, 0.25, 0.5, 0.75, 0.975))
status_sp                <- 1
names(status_sp)         <- "Parametric"
summ_rho_sp              <- c(rho, mean_bstrp_rho_sp, se_bstrp_rho_sp, quantiles_rho_sp, status_sp)

y <- Sys.time()
difftime(y,x,units = "mins")
#------------------------------------------------------------------------------------------------------------------------------

# Compute Survival Correlation Estimate----------------------------------------------------------------------------------------
x <- Sys.time()
gamma <- compute_survival_correlation(ps_data = ps_data,
                                      PFM     = compute_jbin_cjs(prob.f = pf_mark, prob.m = pm_mark, rho)$prob.mf,
                                      PhiF    = phif_mark,
                                      PhiM    = phim_mark)

# Non-Parametric Bootstrap To Estimate SE 
cat("Non-Parametric Bootstrapping to get standard error estimates of gamma...","\n")
gamma_bs_np <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                iter                  = 1000,
                                                                rho                   = rho,
                                                                PF                    = pf_mark,
                                                                PM                    = pm_mark,
                                                                gamma                 = gamma,
                                                                PhiF                  = phif_mark,
                                                                PhiM                  = phim_mark,
                                                                parametric            = F)


# Semi-Parametric Bootstrap To Estimate SE 
cat("Semi-Parametric Bootstrapping to get standard error estimates of gamma...","\n")
gamma_bs_sp <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                iter                  = 1000,
                                                                rho                   = rho,
                                                                PF                    = pf_mark,
                                                                PM                    = pm_mark,
                                                                gamma                 = gamma,
                                                                PhiF                  = phif_mark,
                                                                PhiM                  = phim_mark,
                                                                parametric            = T)

# Collect Results
names(gamma)               <- "Est"
mean_bstrp_gamma_np        <- mean(gamma_bs_np)
mean_bstrp_gamma_sp        <- mean(gamma_bs_sp)
names(mean_bstrp_gamma_np) <- "Est_Btstrp"
names(mean_bstrp_gamma_sp) <- "Est_Btstrp"
se_bstrp_gamma_np          <- sd(gamma_bs_np)
se_bstrp_gamma_sp          <- sd(gamma_bs_sp)
names(se_bstrp_gamma_np)   <- "SE"
names(se_bstrp_gamma_sp)   <- "SE"
quantiles_gamma_np         <- quantile(gamma_bs_np, c(0.025, 0.25, 0.5, 0.75, 0.975))
quantiles_gamma_sp         <- quantile(gamma_bs_sp, c(0.025, 0.25, 0.5, 0.75, 0.975))
status_gamma_np            <- 0
status_gamma_sp            <- 1
names(status_gamma_np)     <- "Parametric"
names(status_gamma_sp)     <- "Parametric"
summ_gamma_np <- c(gamma, mean_bstrp_gamma_np, se_bstrp_gamma_np, quantiles_gamma_np,status_gamma_np)
summ_gamma_sp <- c(gamma, mean_bstrp_gamma_sp, se_bstrp_gamma_sp, quantiles_gamma_sp,status_gamma_sp)
y <- Sys.time()
difftime(y,x,units = "mins")
#---------------------------------------------------------------------------------------------------------------------------

# Return Results------------------------------------------------------------------------------------------------------------

# Gather Correlation Results
summ_corr           <- as.data.frame(rbind(summ_rho_np,summ_rho_sp, summ_gamma_np, summ_gamma_sp))
summ_corr$Parameter <- c("rho","rho","gamma","gamma")
summ_corr           <- summ_corr[,c("Parameter","Est","Est_Btstrp","SE", "2.5%","25%","50%","75%","97.5%", "Parametric")]
rownames(summ_corr) <- NULL
# --------------------------------------------------------------------------------------------------------------------------