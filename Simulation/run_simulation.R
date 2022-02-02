# Prep Code -----------------------------------------------------------------------------------------------

# Grab Libraries
library(parallel)
library(tidyverse)
library(readxl)
library(lubridate)
library(rjags)
library(nimble)

# Set up path and in-line concat fn
`%+%` <- function(a, b) paste0(a, b)
proj_dir <- "/home/mdraghic/projects/def-sbonner/mdraghic/mark_recapture_pair_swap/"
script_dir <- proj_dir %+% "/Scripts/"

# Pull Custom Code
source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "01_fn_model_code.R")
source(script_dir %+% "12_pair_swap_mod_nimble.R")
out_dir <- proj_dir %+% "/Simulation/Output/"

## Options (ECHO FOR LOGS)
options(echo = TRUE)

## Read command line arguments----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
k <- as.numeric(args[1])
pars_mat_name <- args[2]
par_index <- (k >= 1)*1 + (k > 100)*1 + (k > 200)*1+ (k > 300)*1 

cat(k,"\n")
cat(par_index, "\n")
cat(pars_mat_name,"\n")

## Load parameter matrix ---------------------------------------------------------------------------------
pars_list <- readRDS(pars_mat_name)
parameter_list <- pars_list[[par_index]]
saveRDS(parameter_list, out_dir %+% "parameter_list_" %+% k %+% ".rds")

# Generate a Dataset and Format to Pair-Swap, CJS, and JS --------------------------------------------------
jags_data <- sim_dat(parameter_list)
cjs_data  <- format_to_cjs(jags_data)
js_data   <- format_to_js(jags_data)

# Store Data
saveRDS(jags_data, out_dir %+% "jags_data_test_" %+% k %+% ".rds")
saveRDS(cjs_data, out_dir %+% "cjs_data_test_" %+% k %+% ".rds")
saveRDS(js_data, out_dir %+% "js_data_test_" %+% k %+% ".rds")

# RUN MODELS ------------------------------------------------------------------------------------------------

# Number of iterations is universal to all three runs
par_settings <- list(`n.iter` = 1e5,
                     `n.thin` = 50,
                     `n.burn` = 5e4,
                     `n.chains` = 1,
                     `n.adapt` = 5e4)

# RUN 1 CJS MODEL
cjs_jags_params <- c("pF","pM", "phiF", "phiM")
cjs_jags_model <- script_dir %+% "/10_cjs_mod_standard.R"
cjs_run <- run_jags(jags_data = cjs_data,
                    jags_model = cjs_jags_model,
                    jags_params = cjs_jags_params, 
                    par_settings = par_settings, 
                    debug = F)
saveRDS(cjs_run, out_dir %+% "cjs_run_" %+% k %+% ".rds")


# RUN 2 JS MODEL ----------------------------------------------------------------------------------------------
js_jags_params <- c("pF","pM", "phiF", "phiM", "eps")
js_jags_model <- script_dir %+% "/10_js_mod_standard.R"
js_run <- run_jags(jags_data    = js_data,
                   jags_model   = js_jags_model,
                   jags_params  = js_jags_params, 
                   par_settings = par_settings, 
                   debug        = F)
saveRDS(js_run, out_dir %+% "js_run_" %+% k %+% ".rds")

# RUN 3 PAIR-SWAP JS MODEL -------------------------------------------------------------------------------------
nimble_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps", "gl", "gu", "ru", "rl")
CpsMCMC <- compile_pair_swap_nimble(jags_data, params = nimble_params)
ps_run <- run_nimble(CpsMCMC,
                     niter   = par_settings$`n.iter` + (par_settings$`n.adapt` + par_settings$`n.burn`),
                     nburnin = par_settings$`n.burn` + (par_settings$`n.adapt`), 
                     thin    = par_settings$`n.thin`)#,
                     # seed = F)
saveRDS(ps_run, out_dir %+% "ps_run_" %+% k %+% ".rds")

