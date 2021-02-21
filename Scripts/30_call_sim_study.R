library(parallel)
library(dclone)

`%+%` <- function(a, b) paste0(a, b)
script_dir <- getwd() %+% "/Scripts/"

source(script_dir %+% "00_fn_sim_pair_data_rework.R")
source(script_dir %+% "02_fn_model_code.R")
out_dir <- getwd() %+% "/Output/"

k = 5

param_list <- list(
  n = 100, 
  k = k, 
  prop.female = 0.5,
  delta = rep(0.9, k),
  phi.f = rep(0.8, k),
  phi.m = rep(0.8, k),
  gam = rep(0.5, k),
  p.f = rep(0.5, k),
  p.m = rep(0.7, k),
  rho = rep(0.6, k),
  betas = list(beta0 = 0, beta1 = 10),
  rand_sex = F,
  rand_init = F,
  init = rep(1,100)
)

# # Pull individual dataset
jags_data <- do.call(simulate_cr_data, param_list)
cjs_data <- format_to_cjs(jags_data)
# 
# 
# # Multiple Datasets using parallel
# sim_cr_dat(parameter_list = param_list, iterations =  2)


# Run JAGS
#set.seed(42)
jags_data <- sim_cr_dat(parameter_list = param_list, iterations =  100)

## MCMC parameters  
par_settings <- list('n.iter' = 20, 
                     'n.thin' = 1,
                     'n.burn' = 10,
                     'n.chains' = 2,
                     'n.adapt' = 10)

## Jags parameters and model script

# Run Full Model + No Groups
jags_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps")
jags_model <- script_dir %+% "/10_mod_pair_swap.R"


for(i in 1:100){
  dat <- jags_data[[i]]
  ## Run jags in parallel and save results
  jags_samples <- run_jags_parallel(dat, 
                                    jags_model,
                                    jags_params, 
                                    par_settings,
                                    out_dir,
                                    outname = "TESTING_MODEL")
  
}

# Error in checkForRemoteErrors(val) : 
#   2 nodes produced errors; first error: Error in node amating_f[49,4]
# Node inconsistent with parents
# Data isn't being passed through Hidden state correctly - correct this next 
# Initializing graph nodes and adapting chains with 10 iterations ... 
# Error in checkForRemoteErrors(val) : 
#   2 nodes produced errors; first error: Error in node amating_f[18,3]
# Node inconsistent with parents


# To do
# Build logic for first encounter 
# Fix data generation (hidden states)
# Work out math for how transition probabilities are handled (and get their unconditional form)
# work out logic for hidden covariate
# work out logic for penalities based on position to ensure that the mariginal distributions are correctly balanced! (if possible)
# Clean up code to share 
# Add standard CJS Model
# Code that generates standard cjs data (should be archived somewhere)
# Sim study functions (graphs and stuff)
# Clean up the HH duck data (and create general purpose functions too)
# Run that study and compare to standard models
# For both studies will likely need to take a learning detour and figure out sharcnet 
# Writing/conferences/chapter 2 vs journal article
