library(parallel)
library(dclone)
library(tidyverse)
library(readxl)
library(lubridate)

#setwd("C:/Users/Alex/Documents/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/")
`%+%` <- function(a, b) paste0(a, b)
script_dir <- getwd() %+% "/Scripts/"
dat_dir <- getwd() %+% "/Data/FW__Harlequin_Ducks/"

source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "02_fn_model_code.R")
source(script_dir %+% "03_fn_process_hduck_data.R")
out_dir <- getwd() %+% "/Output/"

#HDUCK Data

cap.data <- gather_hq_data(dat_dir) %>% build_cr_df() %>%  add_implied_states() %>% assign_ids_bysex()
cap.data <- cap.data %>% filter(initial_entry < 28)
jags_data <- build_jags_data(cap.data)

# Still need to remove some impossible cases
# - If initial is after current remove from psi matrix (set 1 to 0)
# - should area be considered
# - some pairs need to be manually updated due to inconsistency wrt mate band and animal id in capture record
# - whats a logical "dead" age 
# - some known fate information in the dataset to consider as well (in columns)


#SIM DATA

k = 5
n = 10

param_list <- list(
  n = n, 
  k = k, 
  prop.female = 0.5,
  delta = rep(0.9, k), 
  phi.f = rep(0.8, k),
  phi.m = rep(0.8, k),
  gam = rep(0.6, k),
  p.f = rep(0.75, k),
  p.m = rep(0.75, k),
  rho = rep(0.6, k),
  betas = list(beta0 = 0.0, beta1 = 0.5),
  rand_sex = F,
  rand_init = F,
  init = rep(1,n) 
)

# Pull individual dataset
set.seed(100)
jags_data <- do.call(simulate_cr_data, param_list)
cjs_data <- format_to_cjs(jags_data)

# 
# jags_data$psi <- jags_data$apairs
# jags_data$psi[is.na(jags_data$psi)] <- 1
# jags_data$psi <- jags_data$psi[1:jags_data$nf+1,1:jags_data$nm+1,]
# jags_data$psi <- jags_data$psi[,,1:jags_data$k+1]
# 
# psi <- array(NA,dim = c(jags_data$nf,jags_data$nm+1,jags_data$k))
# 
# for(i in 1:jags_data$k){
#   psi[,,i] <- cbind(jags_data$psi[,,i],rep(0,jags_data$nm))
# }
# 
# jags_data$psi <- psi
# 
# jags_data$apairs_f <- cbind(rep((jags_data$nf+1),nrow(jags_data$apairs_f)),jags_data$apairs_f)
# 
# # 
# 
# # Multiple Datasets using parallel
# sim_cr_dat(parameter_list = param_list, iterations =  2)

# # Run JAGS
# jags_data <- sim_cr_dat(parameter_list = param_list, iterations =  100)
# 
# ## MCMC parameters  
par_settings <- list('n.iter' = 1e4, 
                     'n.thin' = 10,
                     'n.burn' = 1e3,
                     'n.chains' = 2,
                     'n.adapt' = 1e3)

# Vaillancourt 
# ## Jags parameters and model script
# 
# Run standard Model

jags_params <- c("pF", "pM", "phiF", "phiM")
jags_model <- script_dir %+% "/11_cjs_mod_standard.R"


jags_samples <- run_jags_parallel(cjs_data,
                                  jags_model,
                                  jags_params,
                                  par_settings,
                                  out_dir,
                                  outname = "T1_CJS_STD")


# TEST DATA WITH PROGRAM MARK TO SEE RESULTS 

# Run Full Model + No Groups
## MCMC parameters  
par_settings <- list('n.iter' = 1e4, 
                     'n.thin' = 10,
                     'n.burn' = 1e3,
                     'n.chains' = 2,
                     'n.adapt' = 1e3)


jags_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps")
jags_model <- script_dir %+% "/11_mod_pair_swap_notime.R"


jags_samples2 <- run_jags_parallel(jags_data, 
                                  jags_model,
                                  jags_params, 
                                  par_settings,
                                  out_dir,
                                  outname = "TESTING_MODEL2")


gather_posterior_summary(jags_samples) #%>% 
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



# REPARTNER PROCESS and AREPARTNER PROCESS NEED TO BE FIXED (DONE!)
# NIMBLE WONT COMPILE...EIGENVALUE ERRORS 
# JAGS RUNS GREAT
# LARGE DATASETS MAY BE A PROBLEM

# TO DO
# 1. NEED CODE REVIEW FROM SIMON AFTER I DOCUMENT ALL THE STEPS
# 2. DO THE HDUCK CONVERSION 