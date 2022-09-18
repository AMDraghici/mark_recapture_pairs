# Load packages
library(parallel)
library(dclone)
library(tidyverse)
library(readxl)
library(lubridate)
library(rjags)
library(nimble)
library(coda)
library(ggmcmc)

## MCMC parameters
niter <- 1e5
nburnin <- niter/2
nchains <- 5
nthin <- 20

## Load scripts
`%+%` <- function(a, b) paste0(a, b)
setwd("C:/Users/Alex/Documents/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/")
src_dir <- getwd()#"/home/sbonner/Students/Statistics/A_Draghici/Research/mark_recapture_pair_swap"
source(file.path(src_dir,"Scripts","jolly_seber_mod_nimble.R"))
source(file.path(src_dir,"Scripts","pair_swap_mod_nimble4.R"))
source(file.path(src_dir,"Scripts","fn_sim_pair_data2.R"))
# source(file.path(src_dir,"Scripts","fn_process_hduck_data.R"))

# TESTING SIMULATED DATA METHOD -------------------------------------------------------------------------------------------------
# Set number of occasions and animals
k = 30
n = 300

# Seeds for Testing
set.seed(2*pi)
# set.seed(1e5)

# Parameter Grid 
param_list <- list(
  n            = n, # Number of Animals
  k            = k, # Occasions
  lf           = 30, # Data Augmentation for Females (M_F)
  lm           = 30, # Data Augmentation for Males (M_M)
  prop.female  = 0.5, # Proportion of simulated individuals to be female
  delta        = rep(1, k), # Probability that mating is attempted
  phi.f        = rep(0.9, k), # Marginal Prob of Female Survival
  phi.m        = rep(0.9, k), # Marginal Prob of Male Survival
  gam          = rep(0.9, k), # Correlation in Survival Prob of Mates
  p.f          = rep(0.8, k), # Marginal Prob of Female Recapture
  p.m          = rep(0.8, k), # Marginal Prob of Male Recapture
  rho          = rep(0.8, k), # Correlation in male survival rates
  betas        = list(beta0 = 1000, beta1 = 0.1), # inv.logit(Beta0 + Beta1 * hij) = Prob of reforming a pair from t-1 after hij times together
  rand_init    = F, # Randomize Initial Entry (just leave as F)
  init         = sample(1, n, TRUE), # Initial Entry into population for individual n
  show_unmated = T, # Include unmated observations in attempt to mate step
  data_aug     = T  # Add individuals to data augmentation
)

# Generate One set of Data
ps_data <- sim_dat(param_list) # pair-swap data
# check_sim_output(ps_data)
js_data <- format_to_js(ps_data) 

# 0.021858, 0.021858 0.015658 0.016329 
# 0.18,0.16,0.16

## Compile model PS------------------------------------------------------------------------
nimble_params <- c("PF","PM","rho","PhiF","PhiM","gamma",
                   "eps","gl","gu","ru","rl","NF","NM","xi","Phi00","Phif0","Phim0","Phifm","P00","Pf0","Pm0","Pfm")

# ACCOUNT FOR RECRUITMENT....
x <- Sys.time()
fit <- run_pair_swap_nimble_parallel(data = ps_data, 
                                     params = nimble_params,
                                     niter = niter, 
                                     nthin = nthin, 
                                     nburnin = nburnin,
                                     ncores = nchains)

samples <- fit$samples
inits <- fit$inits
seeds <- fit$seed
y <- Sys.time()

difftime(y,x,units = "hours")

# inits = generate_nimble_init_pairs(ps_data)
# samples <- run_nimble(CmdlMCMC,
#                       niter = niter,
#                       nburnin = niter/2,
#                       nchains = nchains,
#                       thin  = nthin,
#                       inits = inits,
#                       seed = F)

## Compile model JS-----------------------------------------------------------------------

js_params <- c("PF","PM", "PhiF","PhiM", "xi", "NF", "NM")
x <- Sys.time()
fit_js <- run_js_nimble_parallel(data = js_data, 
                              params = js_params,
                              niter = niter, 
                              nthin = nthin, 
                              nburnin = nburnin,
                              ncores = nchains)

samples_js <- fit_js$samples
inits_js <- fit_js$inits
seeds_js <- fit_js$seed
y <- Sys.time()

difftime(y,x,units = "hours")

## Summary

summ <- summary(samples)
round(cbind(summ[[1]][,"Mean"],summ[[2]][,c("2.5%","97.5%")])[1:14,],3)

summ2 <- summary(samples_js)
round(cbind(summ2[[1]][,"Mean"],summ2[[2]][,c("2.5%","97.5%")])[1:6,],3)


## Convergence diagnostics
gelman.diag(samples[,c("NF","NM","PhiF","PhiM","PF","PM", "rho","gamma")])
gelman.diag(samples_js[,c("NF","NM","PhiF","PhiM","PF","PM")])

## Effective sample size
ess <- round(effectiveSize(samples)/(nchains * (niter-nburnin)/nthin),2)
ess_js <- round(effectiveSize(samples_js)/(nchains * (niter-nburnin)/nthin),2)
ess
ess_js

chain <- 4
## Traceplots.
ggs(samples) %>% 
  # filter(Chain == chain) %>%
  ggs_density("gamma") #+ 
  
# geom_hline(yintercept = param_list$gam[1], col = "red", size = 1.5)

ggs(samples) %>% 
  # filter(Chain == chain) %>%
  ggs_traceplot("PhiF") #
suchain <- 4
ggs(samples) %>%
  
  # filter(Chain == chain) %>%
  ggs_traceplot("PM") # + ylim(c(0,1)) 
# # geom_hline(yintercept = 0.7, col = "red", size = 0.5)
# ggs(samples)%>% 
#   # filter(Chain == chain)  %>%
#   ggs_traceplot("delta")  #+ 
#   
# # geom_hline(yintercept = param_list$delta[1], col = "red", size = 1.5)
# ggs(samples) %>% 
#   # filter(Chain == chain) %>% 
#   ggs_traceplot("beta")

ggs(samples) %>% 
  # filter(Chain == chain) %>%
  ggs_traceplot("N")

ggs(samples)%>%
  # filter(Chain == chain)  %>%
  ggs_traceplot("delta")

hduck_run <- list(ps_samples = samples,
                  ps_data    = ps_data,
                  ps_inits   = inits,
                  js_samples = samples_js,
                  js_data    = js_data,
                  js_inits   = inits_js)
saveRDS(hduck_run, "hduck_run.rds")

x <- readRDS("hduck_run.rds")

#HDUCK Data
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

ps_data <- build_nimble_data(cap.data, data_aug = T, lf = 100, lm = 100)
