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
nburnin <- niter/4
nchains <- 1
nthin <- 10

## Load scripts
`%+%` <- function(a, b) paste0(a, b)
# setwd("C:/Users/Alex/Documents/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/")
src_dir <- getwd()#"/home/sbonner/Students/Statistics/A_Draghici/Research/mark_recapture_pair_swap"
# source(file.path(src_dir,"Scripts","jolly_seber_mod_nimble.R"))
source(file.path(src_dir,"Scripts","pair_swap_mod_nimble9.R"))
source(file.path(src_dir,"Scripts","cormack_jolly_seber_mod_nimble.R"))
source(file.path(src_dir,"Scripts","fn_sim_pair_data3.R"))
# source(file.path(src_dir,"Scripts","fn_process_hduck_data.R"))

# TESTING SIMULATED DATA METHOD -------------------------------------------------------------------------------------------------
# Set number of occasions and animals
k = 30
n = 300

# Seeds for Testing
set.seed(pi)
# set.seed(1e5)

# Parameter Grid 
param_list <- list(
  n            = n, # Number of Animals
  k            = k, # Occasions
  prop.female  = 0.45, # Proportion of simulated individuals to be female
  delta        = rep(1, k), # Probability that mating is attempted
  phi.f        = rep(0.8, k), # Marginal Prob of Female Survival
  phi.m        = rep(0.8, k), # Marginal Prob of Male Survival
  gam          = rep(-0.1, k), # Correlation in Survival Prob of Mates
  p.f          = rep(0.7, k), # Marginal Prob of Female Recapture
  p.m          = rep(0.7, k), # Marginal Prob of Male Recapture
  rho          = rep(0.8, k), # Correlation in male survival rates
  betas        = list(beta0 = 1e3, beta1 = 0), # inv.logit(Beta0 + Beta1 * hij) = Prob of reforming a pair from t-1 after hij times together
  rand_init    = F, # Randomize Initial Entry (just leave as F)
  init         = sample(1:29, n, TRUE), # Initial Entry into population for individual n
  show_unmated = T # Include unmated observations in attempt to mate step
)

# Generate One set of Data
ps_data <- sim_dat(param_list) # pair-swap data
cjs_data <- format_to_cjs(ps_data)
x <- lapply(1:k, function(t) rowSums(ps_data$psi[,1:(ps_data$nm),t]))


## Compile model PS------------------------------------------------------------------------
nimble_params <- c("PF","PM","PhiF","PhiM","beta0","beta1",
                   "gl","gu","gamma",
                   "ru","rl","rho")

# ACCOUNT FOR RECRUITMENT....
x <- Sys.time()
fit <- run_pair_swap_nimble_parallel(data    = ps_data, 
                                     params  = nimble_params,
                                     niter   = niter, 
                                     nthin   = nthin, 
                                     nburnin = nburnin,
                                     ncores  = nchains)

samples <- fit$samples
inits <- fit$inits
seeds <- fit$seed
y <- Sys.time()

difftime(y,x,units = "hours")

## Compile model CJS-----------------------------------------------------------------------

cjs_params <- c("PF","PM", "PhiF","PhiM")
x <- Sys.time()
fit_cjs <- run_cjs_nimble_parallel(data = cjs_data,
                                  params = cjs_params,
                                  niter = niter,
                                  nthin = nthin,
                                  nburnin = nburnin,
                                  ncores = nchains)

samples_cjs <- fit_cjs$samples
inits_cjs <- fit_cjs$inits
seeds_cjs <- fit_cjs$seed
y <- Sys.time()

difftime(y,x,units = "hours")

## Summary

summ <- summary(samples)
round(cbind(summ[[1]][,"Mean"],summ[[2]][,c("2.5%","97.5%")]),3)
# 
summ2 <- summary(samples_cjs)
round(cbind(summ2[[1]][,"Mean"],summ2[[2]][,c("2.5%","97.5%")]),3)


## Convergence diagnostics
gelman.diag(samples)
gelman.diag(samples_cjs)

## Effective sample size
ess <- round(effectiveSize(samples)/(nchains * (niter-nburnin)/nthin),2)
ess_cjs <- round(effectiveSize(samples_cjs)/(nchains * (niter-nburnin)/nthin),2)
ess
ess_cjs

chain <- 3
## Traceplots.
p1 <- ggs(samples) %>% 
  # filter(Chain == chain) %>%
  ggs_density("gamma") +
  geom_vline(xintercept = param_list$gam[1], col = "red")   + xlim(c(-1,1))

p2 <- ggs(samples) %>% 
  # filter(Chain == chain) %>%
  ggs_density("rho") +
  geom_vline(xintercept = param_list$rho[1], col = "red") + xlim(c(-1,1))
gridExtra::grid.arrange(p1,p2,nrow = 2)

  

# geom_hline(yintercept = param_list$gam[1], col = "red", size = 1.5)

p1 <- ggs(samples) %>% 
  # filter(Chain == chain) %>%
  # ggs_traceplot("PF") +
  ggs_density("PF") +
  geom_vline(xintercept = param_list$p.f[1], col = "red") + xlim(c(0,1))
p2 <- ggs(samples) %>%
  # filter(Chain == chain) %>%
  ggs_density("PM") +
  geom_vline(xintercept = param_list$p.m[1], col = "red") + xlim(c(0,1))

gridExtra::grid.arrange(p1,p2,nrow=2)

p1 <- ggs(samples) %>% 
  # filter(Chain == chain) %>%
  ggs_density("PhiF") +
  geom_vline(xintercept = param_list$phi.f[1], col = "red") + xlim(c(0,1))
p2 <- ggs(samples) %>%
  # filter(Chain == chain) %>%
  ggs_density("PhiM") +
  geom_vline(xintercept = param_list$phi.m[1], col = "red") + xlim(c(0,1))

gridExtra::grid.arrange(p1,p2,nrow=2)


p1 <- ggs(samples) %>% 
  # filter(Chain == chain) %>%
  ggs_density("beta0")# +
  # geom_vline(xintercept = param_list$phi.f[1], col = "red") + xlim(c(0,1))
p2 <- ggs(samples) %>%
  # filter(Chain == chain) %>%
  ggs_density("beta1") #+
  # geom_vline(xintercept = param_list$phi.m[1], col = "red") + xlim(c(0,1))

gridExtra::grid.arrange(p1,p2,nrow=2)

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
