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
niter <- 3e4
nburnin <- niter/2
nchains <- 5
nthin <- 5

## Load scripts
`%+%` <- function(a, b) paste0(a, b)
src_dir <- getwd()#"/home/sbonner/Students/Statistics/A_Draghici/Research/mark_recapture_pair_swap"
source(file.path(src_dir,"Scripts","jolly_seber_mod_nimble.R"))
source(file.path(src_dir,"Scripts","pair_swap_mod_nimble3.R"))
source(file.path(src_dir,"Scripts","fn_sim_pair_data.R"))

# TESTING SIMULATED DATA METHOD -------------------------------------------------------------------------------------------------
# Set number of occasions and animals
k = 30
n = 300

# Seeds for Testing
set.seed(30)
# set.seed(1e5)

# Parameter Grid 
param_list <- list(
  n            = n, # Number of Animals
  k            = k, # Occasions
  lf           = 35, # Data Augmentation for Females (M_F)
  lm           = 35, # Data Augmentation for Males (M_M)
  prop.female  = 0.5, # Proportion of simulated individuals to be female
  delta        = rep(1, k), # Probability that mating is attempted
  phi.f        = rep(0.9, k), # Marginal Prob of Female Survival
  phi.m        = rep(0.9, k), # Marginal Prob of Male Survival
  gam          = rep(0, k), # Correlation in Survival Prob of Mates
  p.f          = rep(0.8, k), # Marginal Prob of Female Recapture
  p.m          = rep(0.8, k), # Marginal Prob of Male Recapture
  rho          = rep(0.8, k), # Correlation in male survival rates
  betas        = list(beta0 = 90, beta1 = 0), # inv.logit(Beta0 + Beta1 * hij) = Prob of reforming a pair from t-1 after hij times together
  rand_init    = F, # Randomize Initial Entry (just leave as F)
  init         = sample(1, n, TRUE), # Initial Entry into population for individual n
  show_unmated = T, # Include unmated observations in attempt to mate step
  data_aug     = T  # Add individuals to data augmentation
)

# Generate One set of Data
ps_data <- sim_dat(param_list) # pair-swap data
# check_sim_output(ps_data)
js_data <- format_to_js(ps_data) 

# # Known settings
# ps_data$zf <- c(rep(1,ps_data$true_pop_f),rep(0, ps_data$lf))
# ps_data$zm <- c(rep(1,ps_data$true_pop_m),rep(0, ps_data$lm))

#
# ps_data$amating_f <- rbind(ps_data$mating_f[1:ps_data$true_pop_f,],matrix(0,nrow = ps_data$lf+1, ncol = ps_data$k))
# ps_data$amating_m <- rbind(ps_data$mating_m[1:ps_data$true_pop_m,],matrix(0,nrow = ps_data$lm+1, ncol = ps_data$k))
# ps_data$arepartner <- rbind(ps_data$repartner[,2:ps_data$k], matrix(0,nrow = ps_data$lf, ncol = ps_data$k-1))
# ps_data$recruit_f <- rbind(ps_data$recruit_f_true[1:ps_data$true_pop_f,],matrix(1,nrow = ps_data$lf, ncol = ps_data$k))
# ps_data$recruit_m <- rbind(ps_data$recruit_m_true[1:ps_data$true_pop_m,],matrix(1,nrow = ps_data$lm, ncol = ps_data$k))
# ps_data$af        <- rbind(ps_data$sf[1:ps_data$true_pop_f,],matrix(1,nrow = ps_data$lf, ncol = ps_data$k),rep(0,ps_data$k))
# ps_data$am        <- rbind(ps_data$sm[1:ps_data$true_pop_m,],matrix(1,nrow = ps_data$lm, ncol = ps_data$k),rep(0,ps_data$k))


# RHO FAILS WHEN EVERYTHING BUT AMATING/ZF (AND AREREPARTNER REMOVED) is on
# COULD DA be the problem? 

# DA Seems fine...
# Amating_F/amating_m included corrupts the result...

# TRY WITHOUT PSI (not including throws error at initial value generation - must check this)

# Turn on amating only - why does it fail?

# for(i in 1:nrow(ps_data$amating_f)){
#   x <- max(which(1 == ps_data$amating_f[i,]))
#   if(abs(x)==Inf) next
#   ps_data$af[i,1:x] <- 1
#   rm(x)
# }
# 
# 
# for(i in 1:nrow(ps_data$amating_m)){
#   x <- max(which(1 == ps_data$amating_m[i,]))
#   if(abs(x)==Inf) next
#   ps_data$am[i,1:x] <- 1
#   rm(x)
# }
# for(i in 1:ps_data$true_pop_f){
#   ps_data$amating_f[i,] <- ifelse(ps_data$recruit_f[i,]==0,NA,ps_data$amating_f[i,])
# }
#
# for(i in 1:ps_data$true_pop_m){
#   ps_data$amating_m[i,] <- ifelse(ps_data$recruit_m[i,]==0,NA,ps_data$amating_m[i,])
# }
#


# # ADDING KNOWN PAIRS
# psi_known <- array(0, dim = dim(ps_data$psi))
# psi_known[,dim(psi_known)[2],] <- 0
# 
# for(t in 1:ps_data$k){
#   for(i in 1:ps_data$true_pop_f){
#     for(j in 1:ps_data$true_pop_m){
#       if(ps_data$pairs_f[i,t]==j){
#         psi_known[i,j,t] <- 1
#       }
#     }
#   }
# }
# 
# 
# # psi_known <- ps_data$psi
# # 
# # for(t in 1:ps_data$k){
# #   for(i in 1:ps_data$nf){
# #     for(j in 1:ps_data$nm){
# #       if((is.na(ps_data$zf[i])|is.na(ps_data$zm[j])) & is.na(ps_data$apf[i,t])){
# #         psi_known[i,j,t] <- 1
# #       }
# #     }
# #   }
# # }
# 
# ps_data$psi <- psi_known

## Compile model PS------------------------------------------------------------------------
nimble_params <- c("PF","PM","rho","PhiF","PhiM","gamma",#"beta0", #"beta1",
                   "eps","gl","gu","ru","rl","NF","NM","xi")

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
round(cbind(summ[[1]][,"Mean"],summ[[2]][,c("2.5%","97.5%")])[1:6,],3)

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

chain <- 1
## Traceplots.
ggs(samples) %>% 
  filter(Chain == chain) %>%
  ggs_traceplot("gamma") #+ 
  
# geom_hline(yintercept = param_list$gam[1], col = "red", size = 1.5)

ggs(samples) %>% 
  # filter(Chain == chain) %>%
  ggs_density("rho") #
chain <- 4
ggs(samples_js) %>%
  # filter(Chain == chain) %>%
  ggs_density("P")# + ylim(c(0,1)) # + 
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
  ggs_traceplot("xi")

# Store Results
saveRDS(samples, "samples_simple_80gam_0rho_300_30.rds")
saveRDS(inits, "inits_simple_80gam_0rho_300_30.rds")
saveRDS(ps_data, "ps_data_simple_80gam_0rho_300_30.rds")
saveRDS(seeds, "seeds_simple_80gam_0rho_300_30.rds")


saveRDS(samples_js, "samples_js_80gam_0rho_simple_300_30.rds")
saveRDS(inits_js, "inits_js_80gam_0rho_simple_300_30.rds")
saveRDS(js_data, "js_data_80gam_0rho_simple_300_30.rds")
saveRDS(seeds_js, "seeds_js_80gam_0rho_simple_300_30.rds")
