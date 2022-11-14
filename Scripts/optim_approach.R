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
src_dir <- getwd()
source(file.path(src_dir,"Scripts","cormack_jolly_seber_mod_nimble.R"))
source(file.path(src_dir,"Scripts","fn_sim_pair_data3.R"))


partial_likelihood <- function(pars,
                               PF,
                               PM,
                               nf,
                               nm,
                               k,
                               first_capture_f,
                               recap_f,
                               recap_m,
                               af,
                               am,
                               apairs_f){
  rho <- pars
  recap_dist <- compute_jbin_cjs(prob.f = PF, 
                                 prob.m = PM, 
                                 corr   = rho)
  
  obs_dist <- c(recap_dist$prob.mf,
                recap_dist$prob.f0,
                recap_dist$prob.m0,
                recap_dist$prob.00)
  
  if(round(sum(obs_dist),1) != 1) browser()
  
  lli <- log(obs_dist)
  ll <- 0
  
  for(i in 1:nf){
    for(t in (first_capture_f[i]+1):(k-1)){
      if(!is.na(apairs_f[i,t]) & !is.na(apairs_f[i,t-1]) & !is.na(apairs_f[i,t+1])){
        if(apairs_f[i,t] != (nm+1)){
          if(!is.na(af[i,t]) & !is.na(am[apairs_f[i,t],t])){
            if(af[i,t] == 1 & am[apairs_f[i,t],t] == 1){
              if(apairs_f[i,t] == apairs_f[i,t-1] & apairs_f[i,t] == apairs_f[i,t+1]){
                if(recap_f[i,t] == 1 & recap_m[apairs_f[i,t],t] == 1){
                  ll <- ll + lli[1]
                  next
                }
                
                if(recap_f[i,t] == 1 & recap_m[apairs_f[i,t],t] == 0){
                  ll <- ll + lli[2]
                  next
                }
                
                if(recap_f[i,t] == 0 & recap_m[apairs_f[i,t],t] == 1){
                  ll <- ll + lli[3]
                  next
                }
                
                if(recap_f[i,t] == 0 & recap_m[apairs_f[i,t],t] == 0){
                  ll <- ll + lli[4]
                  next
                }
              }
            }
          }
        } 
      }
    }
  }
  
  return(-ll)
}

# TESTING SIMULATED DATA METHOD -------------------------------------------------------------------------------------------------
# Set number of occasions and animals
PF <- 0.45
PM <- 0.45
PhiF <- 0.8
PhiM <- 0.8
gam_true <- 0.25
rho_true <- 0.4
PFM <- compute_jbin_cjs(PF,PM,rho_true)$prob.mf 
PhiMF <- compute_jbin_cjs(PhiF,PhiM,gam_true)$prob.mf
prob_prod <- PFM * PhiMF
n_pop <- 350
k <- 30

# Parameter Grid 
param_list <- list(
  n            = n_pop, # Number of Animals
  k            = k, # Occasions
  prop.female  = 0.5, # Proportion of simulated individuals to be female
  delta        = rep(1, k), # Probability that mating is attempted
  phi.f        = rep(PhiF, k), # Marginal Prob of Female Survival
  phi.m        = rep(PhiM, k), # Marginal Prob of Male Survival
  gam          = rep(gam_true, k), # Correlation in Survival Prob of Mates
  p.f          = rep(PF, k), # Marginal Prob of Female Recapture
  p.m          = rep(PM, k), # Marginal Prob of Male Recapture
  rho          = rep(rho_true, k), # Correlation in male survival rates
  betas        = list(beta0 = 1000, beta1 = 1000), # inv.logit(Beta0 + Beta1 * hij) = Prob of reforming a pair from t-1 after hij times together
  rand_init    = F, # Randomize Initial Entry (just leave as F)
  init         = sample(1, n_pop, TRUE), # Initial Entry into population for individual n
  show_unmated = T # Include unmated observations in attempt to mate step
)

# ps_data_list <- lapply(1:10, function(x) sim_dat(param_list))
# saveRDS(ps_data_list, "~/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/Simulation_Data_Gamma2.rds")
ps_data_list<- readRDS("~/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/Simulation_Data_Gamma2.rds")

recapture_correlations <- list()
for(i in 1:length(ps_data_list)){
  
  # Generate One set of Data
  ps_data <- ps_data_list[[i]] 
  
  rl <- compute_jbin_param_cjs(PF,PM)$cor_lower_bound
  ru <- compute_jbin_param_cjs(PF,PM)$cor_upper_bound
  
  lower = c(rl)
  upper = c(ru)
  
  recapture_correlations[[i]] <- nlminb(start          = runif(1, min = lower,max = upper),
                                        objective       = partial_likelihood,
                                        PF              = param_list$p.f[1],
                                        PM              = param_list$p.m[1],
                                        nf              = ps_data$nf,
                                        k               = ps_data$k,
                                        nm              = ps_data$nm,
                                        first_capture_f = ps_data$first_capture_f,
                                        recap_f         = ps_data$recap_f,
                                        recap_m         = ps_data$recap_m,
                                        af              = ps_data$af,
                                        am              = ps_data$am,
                                        apairs_f        = ps_data$apairs_f,
                                        lower           = lower,
                                        upper           = upper)$par
}



plot(density(unlist(recapture_correlations)), main = "Recapture Correlation Estimates", xlab = "Rho_Hat", ylab = "PDF")
abline(v = mean(unlist(recapture_correlations)), lty = 2, col = "black")
abline(v = param_list$rho[1], col = "green")



# mesh <- expand.grid(rho =  seq(rl,ru,by = 0.01))
# 
# nll <- sapply(1:nrow(mesh), function(i) partial_likelihood(pars =          c(mesh$rho[i]),
#                                                            PF              = PF,
#                                                            PM              = PM,
#                                                            nf              = nf,
#                                                            k               = k,
#                                                            nm              = nm,
#                                                            first_capture_f = ps_data$first_capture_f,
#                                                            recap_f         = recap_f,
#                                                            recap_m         = recap_m,
#                                                            af              = af,
#                                                            am              = am,
#                                                            apairs_f        = apairs_f) )
# 
# 
# mesh$nll <- nll
# 
# mesh %>% 
#   ggplot(aes(x = rho, y = -nll)) +
#   geom_line() + 
#   theme(legend.position = "none")
# 
# mesh %>% ggplot(aes(rho,gamma, fill = as.numeric(nll)))+ 
#   geom_tile()
# 
# plotly::plot_ly(x = mesh$gamma, y = mesh$rho, z = exp(-(mesh$nll)))
# plotly::plot_ly(x = mesh$gamma, y = mesh$rho, z = (-(mesh$nll)))
# 

PFM_list <- list()

for(i in 1:length(recapture_correlations)){
  PFM_list[[i]] <- compute_jbin_cjs(PF,PM,recapture_correlations[[i]])$prob.mf 
}

PFM <- unlist(PFM_list)

# SURVIVAL CORRELATION
# TESTING SIMULATED DATA METHOD -------------------------------------------------------------------------------------------------

compute_surv_cor <- function(x, PMF, PHIF, PHIM){
  sigF <- sqrt(PHIF * (1-PHIF))
  sigM <- sqrt(PHIM * (1-PHIM))
  
  gl <- compute_jbin_param_cjs(PHIF,PHIM)$cor_lower_bound
  gu <- compute_jbin_param_cjs(PHIF,PHIM)$cor_upper_bound
  
  y <- ((x/PMF) - (PHIF * PHIM))/(sigF * sigM)
  
  return(pmax(pmin(gu, y),gl))
}


compute_correlation <- function(ps_data, PFM, PhiF, PhiM, known_pairs){
  
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
    for(j in (first_capture_f[i]+1):k){
      
      if(!is.na(apairs_f[i,j-1])){
        
        if(apairs_f[i,j-1] == (nm+1)) next
        
        if(recap_f[i,j-1] == 1 & recap_m[apairs_f[i,j-1],j-1] == 1){
          if(!is.na(apairs_f[i,j])){
            
            if((apairs_f[i,j] == apairs_f[i,j-1])){
              if(recap_f[i,j] == 1 & recap_m[apairs_f[i,j],j] == 1){
                n <- n + 1
                N <- N +1
              } else {
                N <- N+1
              } 
            } else {
              N <- N + 1
            }
          }else {
            N <- N + 1
          } 
        }
      }
    }
  }
  
  obs <- compute_surv_cor(n/N,PFM, PhiF, PhiM)
  
  return(list(obs = obs,
              n = n,
              N = N))
}



gamma_corr <- lapply(1:length(ps_data_list), function(i) compute_correlation(ps_data_list[[i]], 
                                                                             PFM[i],
                                                                             PhiF,
                                                                             PhiM,
                                                                             known_pairs = F))



n  <- sapply(1:length(ps_data_list), function(i) gamma_corr[[i]]$n)
N  <- sapply(1:length(ps_data_list), function(i) gamma_corr[[i]]$N)
gam <- sapply(1:length(ps_data_list), function(i) gamma_corr[[i]]$obs)

PFM_true <- compute_jbin_cjs(PF,PM,rho_true)$prob.mf 

# Construct Boundary Plot
y <- compute_surv_cor(seq(0,1,by = 0.01),PFM_true, PhiF, PhiM)
plot(y = y , x =seq(0,1,by = 0.01), type = "l" , ylab = "Gamma_Hat", xlab = "MLE of Bernoulli Construction", main = "Predicted Values of Gamma from Simulation Study")
abline(v = qbinom(c(0.05,0.95),N, PFM_true * PhiMF)/N, col = "red")
abline(h = 0.5, lty = 2)
points(x = n/N, y = gam)

# 95% CI for Y (MLE)
lb <- qbinom(c(0.025),N, PFM_true * PhiMF)/N
ub <- qbinom(c(0.975),N, PFM_true * PhiMF)/N
1 - mean((n/N) >= ub|(n/N) <= lb)

# 95% CI for gamma (convert by passing Y into survival function - might need to account for transformation w/ jacobian)
gam_lb <- compute_surv_cor(lb, PFM_true, PhiF, PhiM)
gam_ub <- compute_surv_cor(ub, PFM_true, PhiF, PhiM)

# Percent coverage by confidence intervals
pct_coverage <- 1 - mean(gam <= gam_lb | gam >= gam_ub)
pct_coverage

# Prediction of gamma 
pred_gam <- mean(gam) 
pred_gam 

# Density Plot of Gamma
plot(density(gam), main = "Density of Simulated Values of Gamma", xlab = "Gamma_Hat", ylab = "PDF")
abline(v = pred_gam[1], lty = 2)
abline(v = gam_true, lty = 2, col = "red")