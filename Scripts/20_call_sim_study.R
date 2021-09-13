library(parallel)
library(dclone)
library(tidyverse)
library(readxl)
library(lubridate)
library(rjags)

#setwd("C:/Users/Alex/Documents/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/")
`%+%` <- function(a, b) paste0(a, b)
script_dir <- getwd() %+% "/Scripts/"
dat_dir <- getwd() %+% "/Data/RE__Harlequin_duck_data/"

source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "01_fn_model_code.R")
source(script_dir %+% "02_fn_process_hduck_data.R")
out_dir <- getwd() %+% "/Output/"

# #HDUCK Data
# 
# cap.data <- gather_hq_data(dat_dir) %>% build_cr_df() %>%  add_implied_states() %>% assign_ids_bysex()
# cap.data <- cap.data #%>% filter(initial_entry < 28)
# jags_data <- build_jags_data(cap.data)
# cjs_data <- format_to_cjs(jags_data)


# for(i in 1:length(jags_data)){
#   print(names(jags_data[i]))
#   if(is.null(dim(jags_data[[i]]))){
#     print(jags_data[[i]])
#   } else {
#     print(dim(jags_data[[i]]))
#   }
# }

#SIM DATA

k = 8
n = 100

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
  betas = list(beta0 = 1.0, beta1 = 1.5),
  rand_sex = F,
  rand_init = F,
  init = rep(1,n)#sample(k-1, n, TRUE)
)

# # # Pull individual dataset
# # #set.seed(42)
jags_data <- sim_dat(param_list)
cjs_data <- format_to_cjs(jags_data)

# Multiple Datasets using parallel
data_list <- sim_cr_dat(parameter_list = param_list, iterations =  100)
shuffled_list <- replicate_shuffled_data(jags_data, 4)


# # Run JAGS
# jags_data <- sim_cr_dat(parameter_list = param_list, iterations =  100)
# 
# ## MCMC parameters  
# par_settings <- list('n.iter' = 1e4, 
#                      'n.thin' = 10,
#                      'n.burn' = 1e3,
#                      'n.chains' = 2,
#                     'n.adapt' = 1e3)

# Vaillancourt 
# ## Jags parameters and model script
# 
# # Run standard Model
# 
jags_params <- c("pF", "pM", "phiF", "phiM")
jags_model <- script_dir %+% "/10_cjs_mod_standard.R"
 
 
jags_samples <- run_jags_parallel(cjs_data,
                                  jags_model,
                                  jags_params,
                                  par_settings,
                                  out_dir,
                                  save = F, 
                                  outname = "T1_CJS_STD")


# TEST DATA WITH PROGRAM MARK TO SEE RESULTS 

# Run Full Model + No Groups
## MCMC parameters  
par_settings <- list('n.iter' = 10000, 
                     'n.thin' = 10,
                     'n.burn' = 1000,
                     'n.chains' = 1,
                     'n.adapt' = 1000)


jags_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps")
jags_model <- script_dir %+% "/11_mod_pair_swap_notime.R"

x <- run_jags(jags_data = jags_data, 
              jags_model  = jags_model,
              jags_params = jags_params,
              par_settings = par_settings,
              debug = F)


jags_samples2 <- run_jags_parallel(jags_data, 
                                   jags_model,
                                   jags_params, 
                                   par_settings,
                                   out_dir,
                                   outname = "TESTING_MODEL2")


# x <- run_jags(jags_data,
#          jags_model,
#          jags_params, 
#          par_settings)

gather_posterior_summary(jags_samples2) #%>% 
# add_true_values(param_list) %>% 
# plot_caterpillar(params = jags_params) +
# geom_point(aes(x = Parameter, y = true), size = 3, alpha = 0.75, color = "darkblue")

# To do

## ADD RECRUITMENT LOGIC FOR SIMULATION STUDY

# Update Model Code
# - Check that the pair-swap sampling makes sense (DONE)
# - Once pair-swap is reasonable, add histories and repartner logic (DONE) 
# - Then make sure Hduck data is generated correctly (DONE)
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

#k=10,n=100, beta = 3, 12k iter + 2 cores = 70 hours


#SIM DATA

k = 5
n = 50

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
  betas = list(beta0 = 1.0, beta1 = 1.5),
  rand_sex = F,
  rand_init = F,
  init  = rep(1,n)#sample(k-1, n, TRUE)
)

par_settings <- list('n.iter' = 10, 
                     'n.thin' = 1,
                     'n.burn' = 10,
                     'n.chains' = 1,
                     'n.adapt' = 10)

jags_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps")# "psi_raw", "psi_cond", "psi_cond2", "male_taken_jt")
jags_model <- script_dir %+% "/11_mod_pair_swap_notime.R"
jags_data <- sim_dat(param_list)


x <- run_jags(jags_data = jags_data, 
              jags_model  = jags_model,
              jags_params = jags_params,
              par_settings = par_settings,
              debug = F)


jags_data_list <- replicate_shuffled_data(jags_data, 100)
jags_data_list <- sim_cr_dat(param_list, iterations = 10, ncores = 5)

x <- run_jags_simulation_parallel(jags_data_list = jags_data_list,
                                  jags_model  = jags_model,
                                  jags_params = jags_params,
                                  par_settings = par_settings,
                                  out_dir = out_dir,
                                  save = F,
                                  ncores = 5)


# x <- list()
# for(i in 1:20){
#   x[[i]] <- run_jags(jags_data = jags_data, #jags_data_list[[i]],
#                      jags_model  = jags_model, 
#                      jags_params = jags_params, 
#                      par_settings = par_settings,
#                      debug = T)
#   
# }
# 
# 
# jags_data <-  readRDS(getwd() %+% "/jags_data_out_debug.rds") #jags_data_list[[i]]
# jags_init_out <- readRDS(getwd() %+% "/jags_init_out_debug.rds")
# saveRDS(jags_init_out,getwd() %+% "/jags_init_out_debug.rds")
# jags_init <- jags_init_out$jags_inits
# jags_debug <- jags_init_out$jags_debug
# 
# 
# x <- run_jags(jags_data = jags_data,
#               jags_model  = jags_model, 
#               jags_params = jags_params, 
#               par_settings = par_settings,
#               debug = T)
# 
# # apairs_f[10,6]
# # Cannot normalize density 
# 
# 
# for(i in 1:5){
#   print("Time i: "  %+% i %+% "--------------------------------------")
#   for(j in 1:25){
#     print("Female j: "  %+% j)
#     print(which(jags_debug$psi_cond2[j,,i]==1))
#   }
# }
# 
# 
# 
# 
# dim(jags_data$psi)
# 
# jags_data$psi[,nm+2,] <- jags_data$psi[,nm+1,]
