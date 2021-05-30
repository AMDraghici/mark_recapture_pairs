library(parallel)
library(dclone)
library(tidyverse)
library(readxl)
library(lubridate)

`%+%` <- function(a, b) paste0(a, b)
script_dir <- getwd() %+% "/Scripts/"
dat_dir <- getwd() %+% "/Data/FW__Harlequin_Ducks/"

source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "02_fn_model_code.R")
source(script_dir %+% "03_fn_process_hduck_data.R")
out_dir <- getwd() %+% "/Output/"

## HDUCK Data
# 
# cap.data <- gather_hq_data(dat_dir) %>% build_cr_df() %>%  add_implied_states() %>% assign_ids_bysex()
# jags_data <- build_jags_data(cap.data)
# SIM DATA

k = 5
n = 50

param_list <- list(
  n = n, 
  k = k, 
  prop.female = 0.5,
  delta = rep(0.9, k),
  phi.f = rep(0.8, k),
  phi.m = rep(0.8, k),
  gam = rep(0.4, k),
  p.f = rep(0.75, k),
  p.m = rep(0.75, k),
  rho = rep(0.70, k),
  betas = list(beta0 = 0.0, beta1 = 3),
  rand_sex = F,
  rand_init = F,
  init = rep(1,n)
)

# # Pull individual dataset
set.seed(42)
jags_data <- do.call(simulate_cr_data, param_list)
cjs_data <- format_to_cjs(jags_data)
# 
# 
# # Multiple Datasets using parallel
# sim_cr_dat(parameter_list = param_list, iterations =  2)

# # Run JAGS
# jags_data <- sim_cr_dat(parameter_list = param_list, iterations =  100)
# 
# ## MCMC parameters  
# par_settings <- list('n.iter' = 10000, 
#                      'n.thin' = 10,
#                      'n.burn' = 2000,
#                      'n.chains' = 4,
#                      'n.adapt' = 2000)
# Vaillancourt 
# ## Jags parameters and model script
# 
# # Run standard Model
# 
# jags_params <- c("pF", "pM", "phiF", "phiM")
# jags_model <- script_dir %+% "/14_cjs_mod.R"
# 
# 
# jags_samples <- run_jags_parallel(cjs_data, 
#                                   jags_model,
#                                   jags_params, 
#                                   par_settings,
#                                   out_dir,
#                                   outname = "T1_CJS_STD")


# Run Full Model + No Groups
## MCMC parameters  
par_settings <- list('n.iter' = 1e2, 
                     'n.thin' = 1,
                     'n.burn' = 1e2,
                     'n.chains' = 2,
                     'n.adapt' = 1e2)


jags_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps")
jags_model <- script_dir %+% "/10_mod_pair_swap_notime.R"


jags_samples2 <- run_jags_parallel(jags_data, 
                                  jags_model,
                                  jags_params, 
                                  par_settings,
                                  out_dir,
                                  outname = "TESTING_MODEL2")


gather_posterior_summary(jags_samples2) #%>% 
  # add_true_values(param_list) %>% 
  # plot_caterpillar(params = jags_params) +
  # geom_point(aes(x = Parameter, y = true), size = 3, alpha = 0.75, color = "darkblue")

# To do
# Update Model Code
# - Check that the pair-swap sampling makes sense 
# - Once pair-swap is reasonable, add histories and repartner logic 
# - Then make sure Hduck data is generated correctly 
# - If simulation is reasonable look at
# - Data Augementation -> need enough spots for all birds (double check when doing model code) +
#  - -> sometimes birds observed with unseen mate (MEET w/ Simon to Discuss)
# - Speedup by optimizing code
# - Speedup by NIMBLIZING Code
# - Write base R code to generate simulation study datasets
# - Write code to sample iteratively for simulation study (MEET w/ Simon First to Discuss)
# - Code up many time step option (maybe not needed for this study)?

# To do in a while
# Choose good parameters for simulation study + objective 
# Deploy simulation study to sharcnet 
# Write code to build statistics highlighting key aspects of the study
# Once code is fast enough and tested run hduck data through sharcnet 
# Maybe do the chat (or  log-likelihood) thing from chapter 1 to suggest correlation exists (do it after running hduck)
# Pull together results
# Write 
