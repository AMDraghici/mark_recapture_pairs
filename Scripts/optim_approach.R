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
# source(file.path(src_dir,"Scripts","fn_process_hduck_data.R"))

recapture_correlations <- list()
for(i in 1:1000){
  
  # TESTING SIMULATED DATA METHOD -------------------------------------------------------------------------------------------------
  # Set number of occasions and animals
  k = 30
  n = 400
  
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
  ps_data <- sim_dat(param_list) # pair-swap data
  cjs_data <- format_to_cjs(ps_data)
  x <- lapply(1:k, function(t) rowSums(ps_data$psi[,1:(ps_data$nm),t]))
  
  recap_f <- ps_data$recap_f 
  recap_m <- ps_data$recap_m
  apairs_f <- ps_data$apairs_f
  af <- ps_data$af
  am <- ps_data$am
  first_capture_f = ps_data$first_capture_f
  PhiF <- PhiM <- param_list$phi.f[1]
  PF <- PM <- param_list$p.f[1]
  # 
  # partial_likelihood <- function(pars,
  #                                PhiF,
  #                                PhiM,
  #                                PF,
  #                                PM,
  #                                first_capture_f,
  #                                recap_f,
  #                                recap_m,
  #                                apairs_f){
  # 
  #   # gamma <- pars[1]
  #   rho <- pars[1]
  # 
  #   surv_dist <- compute_jbin_cjs(PhiF, PhiM, 0)
  #   recap_dist <- compute_jbin_cjs(PF, PM, rho)
  # 
  #   obs_dist <- c(surv_dist$prob.mf * recap_dist$prob.00 +
  #                   surv_dist$prob.f0 * (1-PF) +
  #                   surv_dist$prob.m0 * (1-PM) +
  #                   surv_dist$prob.00,
  #                 surv_dist$prob.f0 * PF + surv_dist$prob.mf * recap_dist$prob.f0,
  #                 surv_dist$prob.m0 * PM + surv_dist$prob.mf * recap_dist$prob.m0,
  #                 surv_dist$prob.mf * recap_dist$prob.mf)
  # 
  #   if(round(sum(obs_dist),1) != 1) browser()
  # 
  #   lli <- log(obs_dist)
  # 
  #   nm <- dim(recap_m)[1]-1
  #   nf <- dim(recap_f)[1]-1
  #   k  <- dim(recap_f)[2]
  #   ll <- 0
  # 
  #   for(i in 1:nf){
  #     for(t in (first_capture_f[i]+1):k){
  #       if(!is.na(apairs_f[i,t])){
  #         if(apairs_f[i,t] != (nm+1)){
  #           index <- 1 + recap_f[i,t] + 2 * recap_m[apairs_f[i,t],t]
  #           ll <- ll + lli[index]
  #         }
  #       }
  #     }
  #   }
  # 
  #   return(-ll)
  # }
  # # 
  # # 
  # partial_likelihood <- function(pars, 
  #                                PhiF,
  #                                PhiM,
  #                                PF,
  #                                PM,
  #                                first_capture_f,
  #                                recap_f,
  #                                recap_m,
  #                                af,
  #                                am,
  #                                apairs_f){
  #   
  #   gamma <- pars[1]
  #   rho <- pars[2]
  #   
  #   surv_dist <- compute_jbin_cjs(PhiF, PhiM, gamma)
  #   recap_dist <- compute_jbin_cjs(PF, PM, rho)
  # 
  #   surv_cond <- c(surv_dist$prob.mf/PhiM,surv_dist$prob.f0/(1-PhiM))
  #   recap_cond <- c(recap_dist$prob.mf/PM, recap_dist$prob.f0/(1-PM))
  #   
  #  
  #   llsi <- log(surv_cond)
  #   llri <- log(recap_cond)
  #   
  #   nm <- dim(recap_m)[1]-1
  #   nf <- dim(recap_f)[1]-1
  #   k  <- dim(recap_f)[2]
  #   ll <- 0
  #   
  #   for(i in 1:nf){
  #     for(t in (first_capture_f[i]:(k-1))){
  #       if(!is.na(apairs_f[i,t])){
  #         if(apairs_f[i,t] != (nm+1)){
  #           if(is.na(am[apairs_f[i,t],t])){
  #             next
  #           }
  #           if(is.na(af[i,t])){
  #             next
  #           }
  #           ll <- ll + llsi[am[apairs_f[i,t],t] + 1]
  #         }
  #       }
  #     }
  #   }
  #   
  #   
  #   for(i in 1:nf){
  #     for(t in (first_capture_f[i]+1):k){
  #       if(!is.na(apairs_f[i,t])){
  #         if(apairs_f[i,t] != (nm+1)){
  #           ll <- ll + llri[recap_m[apairs_f[i,t],t] + 1]
  #         }
  #       }
  #     }
  #   }
  #   
  #   return(-ll)
  # }
  # 
  
  
  
  partial_likelihood <- function(pars,
                                 PF,
                                 PM,
                                 first_capture_f,
                                 recap_f,
                                 recap_m,
                                 apairs_f){
    rho <- pars[1]
    recap_dist <- compute_jbin_cjs(PF, PM, rho)
    
    obs_dist <- c(recap_dist$prob.00,
                  recap_dist$prob.f0,
                  recap_dist$prob.m0,
                  recap_dist$prob.mf)
    
    if(round(sum(obs_dist),1) != 1) browser()
    
    lli <- log(obs_dist)
    
    nm <- dim(recap_m)[1]-1
    nf <- dim(recap_f)[1]-1
    k  <- dim(recap_f)[2]
    ll <- 0
    
    for(i in 1:nf){
      for(t in (first_capture_f[i]+1):k){
        if(!is.na(apairs_f[i,t])){
          if(apairs_f[i,t] != (nm+1)){
            if(t != k){
              
              if(!is.na(apairs_f[i,t-1])){
                if((apairs_f[i,t-1] == apairs_f[i,t])){
                  index <- 1 + recap_f[i,t] + 2 * recap_m[apairs_f[i,t],t]
                  ll <- ll + lli[index] 
                }
              } else {
                if(!is.na(apairs_f[i,t+1])){
                  if((apairs_f[i,t]==apairs_f[i,t+1])){
                    index <- 1 + recap_f[i,t] + 2 * recap_m[apairs_f[i,t],t]
                    ll <- ll + lli[index] 
                  }
                }
              }
            }
            
            if(t == k){
              if(!is.na(apairs_f[i,t-1])){
                if((apairs_f[i,t-1] == apairs_f[i,t])){
                  index <- 1 + recap_f[i,t] + 2 * recap_m[apairs_f[i,t],t]
                  ll <- ll + lli[index]
                }
              }
            }
          }
        }
      }
    }
    
    return(-ll)
  }
  
  
  # gl <- compute_jbin_param_cjs(PhiF,PhiM)$cor_lower_bound
  # gu <- compute_jbin_param_cjs(PhiF,PhiM)$cor_upper_bound
  
  rl <- compute_jbin_param_cjs(PF,PM)$cor_lower_bound
  ru <- compute_jbin_param_cjs(PF,PM)$cor_upper_bound
  
  
  lower = c(rl)
  upper = c(ru)
  
  # 
  # optim(par             = runif(2,min = lower,max = upper),
  #       fn              = partial_likelihood,
  #       PhiM            = PhiM,
  #       PhiF            = PhiF,
  #       PF              = PF,
  #       PM              = PM,
  #       first_capture_f = first_capture_f,
  #       recap_f         = recap_f,
  #       recap_m         = recap_m,
  #       apairs_f        = apairs_f,
  #       control         = list(factr = 1e5),
  #       method          = "L-BFGS-B",
  #       lower           = lower,
  #       upper           = upper)
  
  recapture_correlations[[i]] <- nlminb(start          = runif(1,min = lower,max = upper),
                                        objective       = partial_likelihood,
                                        PF              = PF,
                                        PM              = PM,
                                        first_capture_f = ps_data$first_capture_f,
                                        recap_f         = recap_f,
                                        recap_m         = recap_m,
                                        apairs_f        = apairs_f,
                                        lower           = lower,
                                        upper           = upper)$par
}


mesh <- expand.grid(rho =  seq(rl,ru,by = 0.01))

nll <- sapply(1:nrow(mesh), function(i) partial_likelihood(pars =  c(mesh$rho[i]),
                                                           PF              = PF,
                                                           PM              = PM,
                                                           first_capture_f = first_capture_f,
                                                           recap_f         = recap_f,
                                                           recap_m         = recap_m,
                                                           apairs_f        = apairs_f) )


mesh$nll <- nll

mesh %>% 
  ggplot(aes(x = rho, y = -nll)) +
  geom_line() + 
  theme(legend.position = "none")

mesh %>% ggplot(aes(rho,gamma, fill = as.numeric(nll)))+ 
  geom_tile()

plotly::plot_ly(x = mesh$gamma, y = mesh$rho, z = exp(-(mesh$nll)))
plotly::plot_ly(x = mesh$gamma, y = mesh$rho, z = (-(mesh$nll)))


# P^MF PHI^MF

compute_surv_cor <- function(x, PMF, PHIF, PHIM){
  sigF <- sqrt(PHIF * (1-PHIF))
  sigM <- sqrt(PHIM * (1-PHIM))
  
  gl <- compute_jbin_param_cjs(PHIF,PHIM)$cor_lower_bound
  gu <- compute_jbin_param_cjs(PHIF,PHIM)$cor_upper_bound
  
  y <- ((x/PMF) - (PHIF * PHIM))/(sigF * sigM)
  
  return(pmax(pmin(gu, y),gl))
}


# TESTING SIMULATED DATA METHOD -------------------------------------------------------------------------------------------------
# Set number of occasions and animals
k = 30
n = 400

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
ps_data <- sim_dat(param_list) # 


N <- 0
n <- 0
apairs_f <- ps_data$apairs_f
nm <- ps_data$nm
nf <- ps_data$nf
k <- ps_data$k
first_capture_f <- ps_data$first_capture_f
recap_f <- ps_data$recap_f
recap_m <- ps_data$recap_m

for(i in 1:nf){
  for(j in (first_capture_f[i]+1):k){
    if(!is.na(apairs_f[i,j-1])){
      
      if(apairs_f[i,j-1] == (nm+1)) next
    
      if(!is.na(apairs_f[i,j])){
        if(apairs_f[i,j] == apairs_f[i,j-1]){
          if(recap_f[i,j] == 1 & recap_m[apairs_f[i,j],j] == 1){
            n <- n + 1
            N <- N +1
          } else {
            N <- N+1
          } 
        } else {
          N <- N + 1
        }
      } else {
        N <- N + 1
      }
    }
  }
}

obs <- compute_surv_cor(n/N,0.65625, 0.8, 0.8)



y <- compute_surv_cor(seq(0,1,by = 0.01),0.65625, 0.8, 0.8)
plot(y = y , x =seq(0,1,by = 0.01), type = "l" )
abline(v = qbinom(c(0.05,0.95),N, 0.8*0.8 * 0.65625)/N, col = "red")
abline(v = n/N, col = "blue")
abline(h = obs, col = "blue")

