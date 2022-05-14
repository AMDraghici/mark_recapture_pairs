# Prepare Environment -------------------------------------------------------------------------------------------------------------------------
#Libraries 
libs <- c("parallel", "dclone", "tidyverse", "readxl", "lubridate", "rjags")

## If package is not installed then do so
for (pkg in libs) {
  if (!pkg %in% rownames(installed.packages())) {
    install.packages(pkg)
  }
}

## Attach our libraries
lapply(libs, require, character.only = TRUE)

# Get directories 
`%+%` <- function(a, b) paste0(a, b)
script_dir <- getwd() %+% "/Scripts/"
dat_dir <- getwd() %+% "/Data/RE__Harlequin_duck_data/"
out_dir <- getwd() %+% "/Output/"

# Source in custom functions 
source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "01_fn_model_code.R")
source(script_dir %+% "02_fn_process_hduck_data.R")

# Data processing ----------------------------------------------------------------------------------------------------------------------------

# Read in and process HDUCK data
cap.data <- gather_hq_data(dat_dir) %>% build_cr_df() %>%  add_implied_states() %>% assign_ids_bysex()
jags_data <- build_jags_data(cap.data)
cjs_data <- format_to_cjs(jags_data)
js_data <- format_to_js(jags_data)

# Store Data
saveRDS(cap.data, out_dir %+% "raw_data_hduck.rds")
saveRDS(jags_data, out_dir %+% "jags_data_pair_swap_hduck.rds")
saveRDS(cjs_data, out_dir %+% "jags_data_std_cjs_hduck.rds")

# Run MCMC Simulation -----------------------------------------------------------------------------------------------------------------------------

# MCMC parameters
par_settings <- list('n.iter' = 1e4,
                     'n.thin' = 10,
                     'n.burn' = 1e3,
                     'n.chains' = 10,
                     'n.adapt' = 1e3)


# Jags parameters and model script

# Run standard Model

# Parameters to save and path to model 
jags_params_std <- c("pF", "pM", "phiF", "phiM")
jags_model_std <- script_dir %+% "/10_cjs_mod_standard.R"

# Run JAGS parallel 
jags_samples_std <- run_jags_parallel(jags_data = cjs_data,
                                      jags_model = jags_model_std,
                                      jags_params = jags_params_std,
                                      par_settings = par_settings,
                                      out_dir = out_dir,
                                      outname = "std_cjs_hduck")

# Run pair-swap model 

# Parameters to save and path to model 
jags_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps")
jags_model <- script_dir %+% "/11_mod_pair_swap_notime.R"

# Run JAGS parallel 
jags_samples <- run_jags_parallel(jags_data = jags_data, 
                                  jags_model = jags_model,
                                  jags_params = jags_params, 
                                  par_settings = par_settings,
                                  out_dir = out_dir,
                                  outname = "pair_swap_hduck")

#--------------------------------------------------------------------------------------------------------------------------------------------