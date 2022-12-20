## Load scripts ------------------------------------------------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
src_dir <- getwd()
source(file.path(src_dir, "Scripts", "fn_generic.R"))
source(file.path(src_dir, "Scripts", "fn_sim_pair_data.R"))
source(file.path(src_dir, "Scripts", "fn_correlation_estimators.R"))
source(file.path(src_dir, "Scripts", "fn_process_hduck_data.R"))

# Load packages
libs <- c("tidyverse","RMark", "readxl", "lubridate")
load_packages(libs, FALSE)

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
# Effective Sample Size CJS Model
n_eff    <- sum(colSums(cjs_data$x[,1:(cjs_data$k-1)]))
# Run CJS Model MARK -----------------------------------------------------------------------------------------------------------
x <- Sys.time()

# Run all 4 possible cases 
cjs_list <- list()
versions <- c("B","S","R","N")

# For loop grabbing model
for(i in 1:length(versions)){
  cat(paste0("Hduck Study: computing CJS estimates for Version:", versions[i], " ..."),"\n")
  cjs_list[[i]] <- run_cjs_model_mark(cjs_data = cjs_data,
                                      title    = "mark_hduck_" %+% versions[i] %+% "_",
                                      version  = versions[i])  
}

cjs_out <- do.call(rbind, cjs_list) %>% mutate(Range95  = UB - LB) 

# Get Predicted Probs for Correlation Estimators (using full model version)
phim_mark <- cjs_out %>% filter(Version == "B" & Parameter == "PhiM") %>% pull(Est)
phif_mark <- cjs_out %>% filter(Version == "B" & Parameter == "PhiF") %>% pull(Est)
pm_mark   <- cjs_out %>% filter(Version == "B" & Parameter == "PM") %>% pull(Est)
pf_mark   <- cjs_out %>% filter(Version == "B" & Parameter == "PF") %>% pull(Est)
y <- Sys.time()
difftime(y,x,units = "mins")
#-------------------------------------------------------------------------------------------------------------------------------

# Compute Recapture Correlation Estimate----------------------------------------------------------------------

# 1. Likelihood Approach -------------------------------------------------------------------------------------
cat(paste0("Hduck Study:  - Estimating recapture correlation, rho, using likelihood approach..."),"\n")
rho_list <- compute_recapture_correlation(ps_data = ps_data, 
                                          PF      = pf_mark,
                                          PM      = pm_mark,
                                          model   = "likelihood")

rho        <- rho_list$rho
n_eff_rho  <- rho_list$n_eff_rho
names(rho) <- "Est"

# Non-Parametric Bootstrap To Estimate SE 
cat(paste0("Hduck Study:  - Non-Parametric Bootstrapping to get standard error estimates of rho (likelihood)..."),"\n")
rho_bs_np <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                               iter       = 1000,
                                                               PF         = pf_mark,
                                                               PM         = pm_mark,
                                                               rho        = rho,
                                                               use_block  = FALSE,
                                                               parametric = FALSE,
                                                               model      = "likelihood")

# Non-Parametric conditional estimator results
summ_rho_np              <- compute_btsrp_summary(rho, rho_bs_np, parameteric = 0, pearson = 0)             

# Parametric Bootstrap To Estimate SE 
cat(paste0("Hduck Study:  - Semi-Parametric Bootstrapping to get standard error estimates of rho (likelihood)..."),"\n")
rho_bs_sp <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                               iter       = 1000,
                                                               PF         = pf_mark,
                                                               PM         = pm_mark,
                                                               rho        = rho,
                                                               use_block  = FALSE,
                                                               parametric = TRUE,
                                                               model      = "likelihood")

# Semi-Parametric conditional estimator
summ_rho_sp              <- compute_btsrp_summary(rho, rho_bs_sp, parameteric = 1, pearson = 0)       

#-------------------------------------------------------------------------------------------------------------

# 2. Pearson Full Approach -----------------------------------------------------------------------------------
cat(paste0("Hduck Study:  - Estimating recapture correlation rho using full pearson ..."),"\n")
pearson_rho <- compute_recapture_correlation(ps_data = ps_data, 
                                             PF      = pf_mark,
                                             PM      = pm_mark,
                                             model   = "full_pearson")$rho
names(pearson_rho) <- "Est"

# Non-Parametric Bootstrap To Estimate SE 
cat(paste0("Hduck Study:  - Non-Parametric Bootstrapping to get standard error estimates of rho (pearson) ..."),"\n")
rho_bs_np_pearson <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                                       iter       = 1000,
                                                                       PF         = pf_mark,
                                                                       PM         = pm_mark,
                                                                       rho        = pearson_rho,
                                                                       use_block  = FALSE,
                                                                       parametric = FALSE,
                                                                       model      = "full_pearson")

# Non-Parametric conditional pearson estimator
summ_rho_np_pearson   <- compute_btsrp_summary(pearson_rho, rho_bs_np_pearson, parameteric = 0, pearson = 1)    

# Parametric Bootstrap To Estimate SE 
cat(paste0("Hduck Study:  - Semi-Parametric Bootstrapping to get standard error estimates of rho (pearson)..."),"\n")
rho_bs_sp_pearson <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                                       iter       = 1000,
                                                                       PF         = pf_mark,
                                                                       PM         = pm_mark,
                                                                       rho        = pearson_rho,
                                                                       use_block  = FALSE,
                                                                       parametric = TRUE,
                                                                       model      = "full_pearson")
# Semi-Parametric conditional pearson estimator
summ_rho_sp_pearson   <- compute_btsrp_summary(pearson_rho, rho_bs_sp_pearson, parameteric = 1, pearson = 1)   

#-------------------------------------------------------------------------------------------------------------

# 2. Pearson Conditional Approach ----------------------------------------------------------------------------
cat(paste0("Hduck Study:  - Estimating recapture correlation rho using pseudo-pearson..."),"\n")
pearson_partial_rho <- compute_recapture_correlation(ps_data = ps_data, 
                                                     PF      = pf_mark,
                                                     PM      = pm_mark,
                                                     model   = "partial_pearson")$rho
names(pearson_partial_rho) <- "Est"

# Non-Parametric Bootstrap To Estimate SE 
cat(paste0("Hduck Study:  - Non-Parametric Bootstrapping to get standard error estimates of rho (Psuedo-Pearson)..."),"\n")
rho_bs_np_pearson_partial <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                                               iter       = 1000,
                                                                               PF         = pf_mark,
                                                                               PM         = pm_mark,
                                                                               rho        = pearson_partial_rho,
                                                                               use_block  = FALSE,
                                                                               parametric = FALSE,
                                                                               model      = "partial_pearson")

# Non-Parametric conditional pearson estimator
summ_rho_np_pearson_partial   <- compute_btsrp_summary(pearson_partial_rho, rho_bs_np_pearson_partial, parameteric = 0, pearson = 2)   

# Parametric Bootstrap To Estimate SE 
cat(paste0("Hduck Study:  - Semi-Parametric Bootstrapping to get standard error estimates of rho (Psuedo-Pearson)..."),"\n")
rho_bs_sp_pearson_partial <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                                               iter       = 1000,
                                                                               PF         = pf_mark,
                                                                               PM         = pm_mark,
                                                                               rho        = pearson_partial_rho,
                                                                               use_block  = FALSE,
                                                                               parametric = TRUE,
                                                                               model      = "partial_pearson")
# Semi-Parametric conditional pearson estimator
summ_rho_sp_pearson_partial   <- compute_btsrp_summary(pearson_partial_rho, rho_bs_sp_pearson_partial, parameteric = 1, pearson = 2)   
#-------------------------------------------------------------------------------------------------------------

# Compute Survival Correlation Estimate-----------------------------------------------------------------------

# Conditional on Rho from Likelihood Approach ----------------------------------------------------------------
cat(paste0("Hduck Study:  - Estimating survival correlation gamma|rho-likelihood..."),"\n")

# If estimate of rho fails (no valid observations) pass dummy values of 10
if(rho == 10){
  gamma <- 10
  gamma_bs_np <- rep(10, 1000)
  gamma_bs_sp <- rep(10, 1000)
} else {
  
  # Estimate Gamma from observed data
  gamma_list <- compute_survival_correlation(ps_data = ps_data,
                                             PFM     = compute_jbin_cjs(prob.f = pf_mark,
                                                                        prob.m = pm_mark,
                                                                        corr   = rho)$prob.mf,
                                             PhiF    = phif_mark,
                                             PhiM    = phim_mark)
  
  # Get Gamma and Effective Sample Size
  gamma       <- gamma_list$gamma
  n_eff_gamma <- gamma_list$n_eff_gamma 
  
  # Non-Parametric Bootstrap To Estimate SE 
  cat(paste0("Hduck Study:  - Non-Parametric Bootstrapping to get standard error estimates of gamma|rho-likelihood..."),"\n")
  gamma_bs_np <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                  iter                  = 1000,
                                                                  rho                   = rho,
                                                                  PF                    = pf_mark,
                                                                  PM                    = pm_mark,
                                                                  gamma                 = gamma,
                                                                  PhiF                  = phif_mark,
                                                                  PhiM                  = phim_mark,
                                                                  use_block             = FALSE,
                                                                  parametric            = FALSE)
  
  
  # Semi-Parametric Bootstrap To Estimate SE 
  cat(paste0("Hduck Study:  - Semi-Parametric Bootstrapping to get standard error estimates of gamma|rho-likelihood..."),"\n")
  gamma_bs_sp <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                  iter                  = 1000,
                                                                  rho                   = rho,
                                                                  PF                    = pf_mark,
                                                                  PM                    = pm_mark,
                                                                  gamma                 = gamma,
                                                                  PhiF                  = phif_mark,
                                                                  PhiM                  = phim_mark,
                                                                  use_block             = FALSE,
                                                                  parametric            = TRUE)
}

# Collect Results
names(gamma) <- "Est"
summ_gamma_np   <- compute_btsrp_summary(gamma, gamma_bs_np, parameteric = 0, pearson = 0)   
summ_gamma_sp   <- compute_btsrp_summary(gamma, gamma_bs_sp, parameteric = 1, pearson = 0)   
#-------------------------------------------------------------------------------------------------------------

# Conditional on Rho from Likelihood Approach ----------------------------------------------------------------
cat(paste0("Hduck Study:  - Estimating survival correlation gamma|rho-pearson..."),"\n")

# If estimate of rho fails (no valid observations) pass dummy values of 10
if(is.na(pearson_rho)|pearson_rho == 10){
  pearson_gamma <- 10
  pearson_gamma_bs_np <- rep(10, 1000)
  pearson_gamma_bs_sp <- rep(10, 1000)
} else {
  
  # Estimate Gamma from observed data
  pearson_gamma <- compute_survival_correlation(ps_data = ps_data,
                                                PFM     = compute_jbin_cjs(prob.f = pf_mark,
                                                                           prob.m = pm_mark,
                                                                           corr   = pearson_rho)$prob.mf,
                                                PhiF    = phif_mark,
                                                PhiM    = phim_mark)$gamma
  
  # Non-Parametric Bootstrap To Estimate SE 
  cat(paste0("Hduck Study:  - Non-Parametric Bootstrapping to get standard error estimates of gamma|rho-pearson..."),"\n")
  pearson_gamma_bs_np <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                          iter                  = 1000,
                                                                          rho                   = pearson_rho,
                                                                          PF                    = pf_mark,
                                                                          PM                    = pm_mark,
                                                                          gamma                 = pearson_gamma,
                                                                          PhiF                  = phif_mark,
                                                                          PhiM                  = phim_mark,
                                                                          use_block             = FALSE,
                                                                          parametric            = FALSE)
  
  
  # Semi-Parametric Bootstrap To Estimate SE 
  cat(paste0("Hduck Study:  - Semi-Parametric Bootstrapping to get standard error estimates of gamma|rho-pearson..."),"\n")
  pearson_gamma_bs_sp <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                          iter                  = 1000,
                                                                          rho                   = pearson_rho,
                                                                          PF                    = pf_mark,
                                                                          PM                    = pm_mark,
                                                                          gamma                 = pearson_gamma,
                                                                          PhiF                  = phif_mark,
                                                                          PhiM                  = phim_mark,
                                                                          use_block             = FALSE,
                                                                          parametric            = TRUE)
}

# Collect Results
names(pearson_gamma) <- "Est"
summ_pearson_gamma_np   <- compute_btsrp_summary(pearson_gamma, pearson_gamma_bs_np, parameteric = 0, pearson = 1)   
summ_pearson_gamma_sp   <- compute_btsrp_summary(pearson_gamma, pearson_gamma_bs_sp, parameteric = 1, pearson = 1) 
#-------------------------------------------------------------------------------------------------------------

# Conditional on Rho from Likelihood Approach ----------------------------------------------------------------
cat(paste0("Hduck Study:  - Estimating survival correlation gamma|rho-psuedo-pearson..."),"\n")

# If estimate of rho fails (no valid observations) pass dummy values of 10
if(is.na(pearson_partial_rho)|pearson_partial_rho == 10){
  pearson_partial_gamma <- 10
  pearson_partial_gamma_bs_np <- rep(10, 1000)
  pearson_partial_gamma_bs_sp <- rep(10, 1000)
} else {
  
  # Estimate Gamma from observed data
  pearson_partial_gamma <- compute_survival_correlation(ps_data = ps_data,
                                                        PFM     = compute_jbin_cjs(prob.f = pf_mark,
                                                                                   prob.m = pm_mark,
                                                                                   corr   = pearson_partial_rho)$prob.mf,
                                                        PhiF    = phif_mark,
                                                        PhiM    = phim_mark)$gamma
  
  # Non-Parametric Bootstrap To Estimate SE 
  cat(paste0("Hduck Study:  - Non-Parametric Bootstrapping to get standard error estimates of gamma|rho-psuedo-pearson..."),"\n")
  pearson_partial_gamma_bs_np <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                                  iter                  = 1000,
                                                                                  rho                   = pearson_partial_rho,
                                                                                  PF                    = pf_mark,
                                                                                  PM                    = pm_mark,
                                                                                  gamma                 = pearson_partial_gamma,
                                                                                  PhiF                  = phif_mark,
                                                                                  PhiM                  = phim_mark,
                                                                                  use_block             = FALSE,
                                                                                  parametric            = FALSE)
  
  
  # Semi-Parametric Bootstrap To Estimate SE 
  cat(paste0("Hduck Study:  - Semi-Parametric Bootstrapping to get standard error estimates of gamma|rho-psuedo-pearson..."),"\n")
  pearson_partial_gamma_bs_sp <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                                  iter                  = 1000,
                                                                                  rho                   = pearson_partial_rho,
                                                                                  PF                    = pf_mark,
                                                                                  PM                    = pm_mark,
                                                                                  gamma                 = pearson_partial_gamma,
                                                                                  PhiF                  = phif_mark,
                                                                                  PhiM                  = phim_mark,
                                                                                  use_block             = FALSE,
                                                                                  parametric            = TRUE)
}

# Collect Results
names(pearson_partial_gamma) <- "Est"
summ_pearson_partial_gamma_np   <- compute_btsrp_summary(pearson_partial_gamma, pearson_partial_gamma_bs_np, parameteric = 0, pearson = 2)   
summ_pearson_partial_gamma_sp   <- compute_btsrp_summary(pearson_partial_gamma, pearson_partial_gamma_bs_sp, parameteric = 1, pearson = 2) 
#-------------------------------------------------------------------------------------------------------------

# Return Results----------------------------------------------------------------------------------------------

cat("Formatting output ...","\n")

# Gather Correlation Results
summ_corr           <- as.data.frame(rbind(summ_rho_np,
                                           summ_rho_sp,
                                           summ_rho_np_pearson,
                                           summ_rho_sp_pearson, 
                                           summ_rho_np_pearson_partial,
                                           summ_rho_sp_pearson_partial,
                                           summ_gamma_np, 
                                           summ_gamma_sp,
                                           summ_pearson_gamma_np,
                                           summ_pearson_gamma_sp,
                                           summ_pearson_partial_gamma_np,
                                           summ_pearson_partial_gamma_sp))

summ_corr$Parameter <- c(rep("rho", 6), rep("gamma",6))
summ_corr           <- summ_corr[,c("Parameter","Est","Est_Btstrp","SE", "2.5%","25%","50%","75%","97.5%", "Parametric", "Pearson")]
summ_corr           <- summ_corr %>% 
  mutate(Bias_Btstrp2_Mean   = Est   - Est_Btstrp,
         Bias_Btstrp2_Median = Est   - `50%`,
         Range95             = `97.5%` - `2.5%`,
         Range50             = `75%` - `25%`,
         Cover095            = 1*(0 <= `97.5%` & 0 >= `2.5%`),    
         Cover050            = 1*(0 <= `75%`   & 0 >= `25%`)) 

rownames(summ_corr) <- NULL

# Compute CChat Adjustment
summ_chat <- data.frame(Method     = c("Likelihood", "Pearson", "Psuedo-Pearson"),
                        CChatRho   = c(compute_proposed_chat(rho, 0), 
                                       compute_proposed_chat(pearson_rho, 0), 
                                       compute_proposed_chat(pearson_partial_rho, 0)),
                        CChatGamma = c(compute_proposed_chat(0, gamma), 
                                       compute_proposed_chat(0, pearson_gamma), 
                                       compute_proposed_chat(0, pearson_partial_gamma)),
                        CChat      = c(compute_proposed_chat(rho, gamma), 
                                       compute_proposed_chat(pearson_rho, pearson_gamma), 
                                       compute_proposed_chat(pearson_partial_rho, pearson_partial_gamma)))

# Add Correlation Chat Correction to Mark-Recapture Results
summ_cjs <- cjs_out %>% 
  mutate(CChatAdj_Likelihood      = ifelse(Parameter == "P", 
                                           compute_proposed_chat(rho, 0), 
                                           ifelse(Parameter == "Phi",
                                                  compute_proposed_chat(0, gamma), 
                                                  1)),
         CChatAdj_Pearson         = ifelse(Parameter == "P", 
                                           compute_proposed_chat(pearson_rho, 0), 
                                           ifelse(Parameter == "Phi",
                                                  compute_proposed_chat(0, pearson_gamma), 
                                                  1)),
         CChatAdj_Partial_Pearson = ifelse(Parameter == "P", 
                                           compute_proposed_chat(pearson_partial_rho, 0), 
                                           ifelse(Parameter == "Phi",
                                                  compute_proposed_chat(0, pearson_partial_gamma), 
                                                  1)),
         UBAdj_Likelihood         = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Likelihood), alpha = 0.05)[["ub"]],
         LBAdj_Likelihood         = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Likelihood), alpha = 0.05)[["lb"]],
         UBAdj_Pearson            = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Pearson), alpha = 0.05)[["ub"]],
         LBAdj_Pearson            = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Pearson), alpha = 0.05)[["lb"]],
         UBAdj_Partial_Pearson    = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Partial_Pearson), alpha = 0.05)[["ub"]],
         LBAdj_Partial_Pearson    = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Partial_Pearson), alpha = 0.05)[["lb"]],
         Range95_Likelihood       = UBAdj_Likelihood - LBAdj_Likelihood,
         Range95_Pearson          = UBAdj_Pearson - LBAdj_Pearson,
         Range95_Partial_Pearson  = UBAdj_Partial_Pearson - LBAdj_Partial_Pearson
  ) 

# Produce Likelihood Ratio Tests
summ_lrt <- compute_lrt_summ(summ_cjs = summ_cjs,
                             n_eff    = n_eff,
                             iter     = 1,
                             scenario = 1)


# Produce AICC comparisons
summ_aic <- compute_aic_summ(summ_cjs = summ_cjs,
                             n_eff    = n_eff,
                             iter     = 1,
                             scenario = 1)

# Summarize Sample Sizes
summ_n <- data.frame(n_eff                       = n_eff,
                     n_eff_rho                   = n_eff_rho,
                     n_eff_gamma                 = n_eff_gamma,
                     iter                        = 1,
                     scenario                    = 1)

# ----------------------------------------------------------------------------------------------------------------------------------------