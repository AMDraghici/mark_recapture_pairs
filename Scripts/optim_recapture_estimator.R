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

## Load scripts
`%+%` <- function(a, b) paste0(a, b)
# setwd("C:/Users/Alex/Documents/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/")
src_dir <- getwd()#"/home/sbonner/Students/Statistics/A_Draghici/Research/mark_recapture_pair_swap"
# source(file.path(src_dir,"Scripts","jolly_seber_mod_nimble.R"))
# source(file.path(src_dir,"Scripts","pair_swap_mod_nimble9.R"))
source(file.path(src_dir,"Scripts","cormack_jolly_seber_mod_nimble.R"))
source(file.path(src_dir,"Scripts","fn_sim_pair_data3.R"))

# TESTING SIMULATED DATA METHOD -------------------------------------------------------------------------------------------------
# Set number of occasions and animals
k = 30
n = 150

# Seeds for Testing
# set.seed(pi)
# set.seed(1e5)

# Parameter Grid 
param_list <- list(
  n            = n, # Number of Animals
  k            = k, # Occasions
  prop.female  = 0.5, # Proportion of simulated individuals to be female
  delta        = rep(1, k), # Probability that mating is attempted
  phi.f        = rep(0.8, k), # Marginal Prob of Female Survival
  phi.m        = rep(0.8, k), # Marginal Prob of Male Survival
  gam          = rep(0, k), # Correlation in Survival Prob of Mates
  p.f          = rep(0.75, k), # Marginal Prob of Female Recapture
  p.m          = rep(0.75, k), # Marginal Prob of Male Recapture
  rho          = rep(0.5, k), # Correlation in male survival rates
  betas        = list(beta0 = 1e3, beta1 = 0), # inv.logit(Beta0 + Beta1 * hij) = Prob of reforming a pair from t-1 after hij times together
  rand_init    = F, # Randomize Initial Entry (just leave as F)
  init         = sample(1, n, TRUE), # Initial Entry into population for individual n
  show_unmated = T # Include unmated observations in attempt to mate step
)

# Generate One set of Data
ps_data <- readRDS("~/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/Simulation_Data_Gamma.rds")
# cjs_data <- format_to_cjs(ps_data)
# x <- lapply(1:k, function(t) rowSums(ps_data$psi[,1:(ps_data$nm),t]))

recap_f <- ps_data$recap_f 
recap_m <- ps_data$recap_m
apairs_f <- ps_data$apairs_f
first_capture_f = ps_data$first_capture_f
PF <- PM <- param_list$p.f[1]
nf <- ps_data$nf
nm <- ps_data$nm 
k  <- ps_data$k    


compute_recap_cor <- function(x, PF, PM){
  sigF <- sqrt(PF * (1-PF))
  sigM <- sqrt(PM * (1-PM))
  
  rl <- compute_jbin_param_cjs(PF,PM)$cor_lower_bound
  ru <- compute_jbin_param_cjs(PF,PM)$cor_upper_bound
  
  y <- (x - (PF * PM))/(sigF * sigM)
  
  return(pmax(pmin(ru, y),rl))
}

plot(x = seq(0,1,by = 0.01), compute_recap_cor(seq(0,1,by = 0.01), 0.75,0.75), type = "l")


compute_correlation <- function(ps_data, PF, PM, known_pairs){
  
  N <- 0
  n <- 0
  if(known_pairs){
    apairs_f <- ps_data$pairs_f 
  } else{
    apairs_f <- ps_data$apairs_f 
  }
  
  nm <- ps_data$nm
  nf <- ps_data$nf
  k <- ps_data$k
  first_capture_f <- ps_data$first_capture_f
  recap_f <- ps_data$recap_f
  recap_m <- ps_data$recap_m
  
  for(i in 1:nf){
  
    
    for(t in (first_capture_f[i]+1):k){
      
      if(!is.na(apairs_f[i,t]) & !is.na(apairs_f[i,t-1])){
        
        if(apairs_f[i,t] == (nm+1)) next
        
        if(apairs_f[i,t] == apairs_f[i,t-1]){
          
          if(recap_f[i,t] == 1 & recap_m[apairs_f[i,t],t] == 1){
            n <- n + 1
            N <- N +1
          } else {
            N <- N+1
          } 
        } 
      }
    }
  }
  
  obs <- compute_recap_cor(n/N,PF, PM)
  
  return(list(obs = obs,
              n = n,
              N = N))
}

# compute_correlation(ps_data, PF, PM, T)

# 
rho_corr_known <- lapply(1:length(ps_data), function(i) compute_correlation(ps_data[[i]],
                                                                              PF,
                                                                              PM,
                                                                              known_pairs = T))


rho_corr_unknown <- lapply(1:length(ps_data), function(i) compute_correlation(ps_data[[i]], 
                                                                              PF,
                                                                              PM,
                                                                              known_pairs = F))


rho_corr <- rho_corr_known

n  <- sapply(1:length(ps_data), function(i) rho_corr[[i]]$n)
N  <- sapply(1:length(ps_data), function(i) rho_corr[[i]]$N)
rho <- sapply(1:length(ps_data), function(i) rho_corr[[i]]$obs)

# Construct Boundary Plot
y <- compute_recap_cor(seq(0,1,by = 0.01), PF, PM)
plot(y = y , x =seq(0,1,by = 0.01), type = "l" , ylab = "rho_Hat",
     xlab = "MLE of Bernoulli Construction", 
     main = "Predicted Values of rho from Simulation Study")
abline(v = qbinom(c(0.05,0.95),N, compute_jbin_cjs(0.75,0.75,0)$prob.mf)/N, col = "red")
abline(h = 0, lty = 2)
points(x = n/N, y = rho)

# 95% CI for Y (MLE)
lb <- qbinom(c(0.025),N, compute_jbin_cjs(0.75,0.75,0)$prob.mf)/N
ub <- qbinom(c(0.975),N, compute_jbin_cjs(0.75,0.75,0)$prob.mf)/N
1 - mean((n/N) >= ub|(n/N) <= lb)

mean(rho)
plot(density(rho))

