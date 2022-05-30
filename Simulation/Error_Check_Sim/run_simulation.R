# Prep Code ------------------------------------------------------------------------------------------------------------------
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
source(script_dir   %+% "00_fn_sim_pair_data.R")
source(script_dir   %+% "01_fn_model_code.R")
source(script_dir   %+% "11_jolly_seber_mod_nimble.R")
source(script_dir   %+% "12_pair_swap_mod_nimble.R")
out_dir <- proj_dir %+% "/Simulation/Output/"

## Options (ECHO FOR LOGS)
options(echo = TRUE)

## Read command line arguments----------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
k <- as.numeric(args[1])
pars_mat_name <- args[2]

par_index <- (k >= 1)*1 + (k > 50)*1 + (k > 100)*1 + (k > 150)*1 

cat(k,"\n")
cat(par_index, "\n")
cat(pars_mat_name,"\n")

## Load parameter matrix ---------------------------------------------------------------------------------------------------
pars_list <- readRDS(pars_mat_name)
parameter_list <- pars_list[[par_index]]
saveRDS(parameter_list, out_dir %+% "parameter_list_" %+% k %+% ".rds")

# Generate a Pair-Swap Dataset and format into Jolly-Seber------------------------------------------------------------------
ps_data <- sim_dat(parameter_list)
js_data <- format_to_js(ps_data)

# RUN MODELS ---------------------------------------------------------------------------------------------------------------

# Jolly-Seber MODEL --------------------------------------------------------------------------------------------------------
js_params <- c("PF","PM","PhiF","PhiM", "eps","xi")

# COMPILE NIMBLE JOLLY-SEBER MODEL
start <- Sys.time()
CjsMCMC_List <- compile_jolly_seber_nimble(js_data = js_data, params = js_params)
end <- Sys.time()
print(start-end)

# RUN MODEL AND STORE RESULTS
start <- Sys.time()
js_samples <- run_nimble(CmdlMCMC = CjsMCMC_List$CjsMCMC,
                         niter    = 1e5,
                         nburnin  = 5e4,
                         thin     = 10)
end <- Sys.time()
print(start-end)

# Build output list 
js_run <- list(data        = js_data,
               inits       = CjsMCMC_List$nimble_inits,
               samples     = js_samples)

# Store kth iter 
saveRDS(js_run, out_dir %+% "js_run_" %+% k %+% ".rds")

# RUN 3 PAIR-SWAP JS MODEL -------------------------------------------------------------------------------------------------
ps_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps", "gl", "gu", "ru", "rl","NF","NM","xi")

start <- Sys.time()
CpsMCMC_List <- compile_pair_swap_nimble(ps_data = ps_data, params = ps_params)
end <- Sys.time()
print(start-end)

# RUN MODEL AND STORE RESULTS
start <- Sys.time()
ps_samples <- run_nimble(CmdlMCMC = CpsMCMC_List$CpsMCMC,
                         niter    = 1e5,
                         nburnin  = 5e4, 
                         thin     = 10)
end <- Sys.time()
print(start-end)

# Build output list 
ps_run <- list(data        = ps_data,
               inits       = CpsMCMC_List$nimble_inits,
               samples     = ps_samples)

saveRDS(ps_run, out_dir %+% "ps_run_" %+% k %+% ".rds")
#--------------------------------------------------------------------------------------------------------------------------
