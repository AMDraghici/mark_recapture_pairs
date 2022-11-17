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


partial_likelihood2 <- function(pars,
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
          # if(!is.na(af[i,t]) & !is.na(am[apairs_f[i,t],t])){
            # if(af[i,t] == 1 & am[apairs_f[i,t],t] == 1){
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
      # }
    # }
  }
  
  return(-ll)
}

compute_recapture_correlation <- function(x, PF, PM){
    sigF <- sqrt(PF * (1-PF))
    sigM <- sqrt(PM * (1-PM))
    
    rl <- compute_jbin_param_cjs(PF,PM)$cor_lower_bound
    ru <- compute_jbin_param_cjs(PF,PM)$cor_upper_bound
    
    y <- (x - (PF * PM))/(sigF * sigM)
    
    return(pmax(pmin(ru, y),rl))
  }
  
  
compute_rho <- function(PF,
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
 
  
  n1 <- 0
  n2 <- 0
  n3 <- 0
  n4 <- 0 
  N  <- 0
  
  for(i in 1:nf){
    for(t in (first_capture_f[i]+1):(k-1)){
      if(!is.na(apairs_f[i,t]) & !is.na(apairs_f[i,t-1]) & !is.na(apairs_f[i,t+1])){
        if(apairs_f[i,t] != (nm+1)){
          if(apairs_f[i,t] == apairs_f[i,t-1] & apairs_f[i,t] == apairs_f[i,t+1]){
            if(recap_f[i,t] == 1 & recap_m[apairs_f[i,t],t] == 1){
              n1 <- n1+1
              N <- N+1
            }
            
            if(recap_f[i,t] == 1 & recap_m[apairs_f[i,t],t] == 0){
              n2 <- n2+1
              N <- N+1
            }
            
            if(recap_f[i,t] == 0 & recap_m[apairs_f[i,t],t] == 1){
              n3 <- n3+1
              N <- N+1
            }
            
            if(recap_f[i,t] == 0 & recap_m[apairs_f[i,t],t] == 0){
              n4 <- n4+1
              N <- N+1
            }
          }
        }
      }
    } 
  }
  
  rho <- compute_recapture_correlation(n1/N, PF, PM)
  
  return(list(rho = rho,
              n1   = n1,
              n2   = n2,
              n3   = n3,
              n4   = n4,
              N   = N))
}


# TESTING SIMULATED DATA METHOD -------------------------------------------------------------------------------------------------
# Set number of occasions and animals
PF <- 0.75
PM <- 0.75
PhiF <- 0.8
PhiM <- 0.8
gam_true <- 0.29
rho_true <- 0.4
PFM <- compute_jbin_cjs(PF,PM,rho_true)$prob.mf 
PhiMF <- compute_jbin_cjs(PhiF,PhiM,gam_true)$prob.mf
prob_prod <- PFM * PhiMF
n_pop <- 300
k <- 20
set.seed(1)

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

ps_data_list <- lapply(1:25, function(x) sim_dat(param_list))
# saveRDS(ps_data_list, "~/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/Simulation_Data_Gamma4.rds")
ps_data_list<- readRDS("~/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/Simulation_Data_Gamma4.rds")

recapture_correlations <- list()
recapture_correlations_2 <- list()
for(i in 1:length(ps_data_list)){
  
  # Generate One set of Data
  ps_data <- ps_data_list[[i]] 
  
  rl <- compute_jbin_param_cjs(PF,PM)$cor_lower_bound
  ru <- compute_jbin_param_cjs(PF,PM)$cor_upper_bound
  
  lower = c(rl)
  upper = c(ru)
  
  recapture_correlations[[i]] <- nlminb(start           = runif(1, min = lower,max = upper),
                                        objective       = partial_likelihood2,
                                        PF              = param_list$p.f[1],
                                        PM              = param_list$p.m[1],
                                        nf              = ps_data$nf,
                                        k               = ps_data$k,
                                        nm              = ps_data$nm,
                                        first_capture_f = ps_data$first_capture_f,
                                        recap_f         = ps_data$recap_f,
                                        recap_m         = ps_data$recap_m,
                                        apairs_f        = ps_data$apairs_f,
                                        lower           = lower,
                                        upper           = upper)$par
  
  
  recapture_correlations_2[[i]] <- nlminb(start          = runif(1, min = lower,max = upper),
                                          objective       = partial_likelihood,
                                          PF              = param_list$p.f[1],
                                          PM              = param_list$p.m[1],
                                          nf              = ps_data$nf,
                                          k               = ps_data$k,
                                          nm              = ps_data$nm,
                                          first_capture_f = ps_data$first_capture_f,
                                          recap_f         = ps_data$recap_f,
                                          recap_m         = ps_data$recap_m,
                                          apairs_f        = ps_data$apairs_f,
                                          lower           = lower,
                                          upper           = upper)$par
  
}



plot(density(unlist(recapture_correlations)), main = "Recapture Correlation Estimates", xlab = "Rho_Hat", ylab = "PDF", xlim = c(-0.1,1), ylim = c(0,6))
abline(v = mean(unlist(recapture_correlations)), lty = 2, col = "black")
abline(v = param_list$rho[1], col = "green")

lines(density(unlist(recapture_correlations_2)), col = "red", lty = 2)
abline(v = mean(unlist(recapture_correlations_2)), lty = 2, col = "red")



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

for(i in 1:length(x)){
  PFM_list[[i]] <- compute_jbin_cjs(PF,PM,x[[i]])$prob.mf 
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

partial_likelihood_gam <- function(pars, ps_data, PFM, PhiF, PhiM){
  
  gamma <- pars
  surv_dist <- compute_jbin_cjs(prob.f = PhiF, 
                                prob.m = PhiM, 
                                corr   = gamma)
  
  prob <- surv_dist$prob.mf * PFM
  lli <- log(c(prob,1-prob))
  ll <- 0
  
  # Grab data
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
        
        if(recap_f[i,j-1] == 1 & recap_m[apairs_f[i,j-1],j-1] == 1){
          if(!is.na(apairs_f[i,j])){
            
            if((apairs_f[i,j] == apairs_f[i,j-1])){
              if(recap_f[i,j] == 1 & recap_m[apairs_f[i,j],j] == 1){
                ll <- ll + lli[1]
              } else {
                ll <- ll + lli[2]
              } 
            } else {
              ll <- ll + lli[2]
            }
          }else {
            ll <- ll + lli[2]
          } 
        }
      }
    }
  }
  
  return(-ll)
}


lower <- compute_jbin_param_cjs(PhiF,PhiM)$cor_lower_bound
upper <- compute_jbin_param_cjs(PhiF,PhiM)$cor_upper_bound

nlminb(start          = runif(1, min = lower,max = upper),
       objective       = partial_likelihood_gam,
       ps_data         = ps_data_list[[10]],
       PFM             = PFM[10],
       PhiF            = PhiF, 
       PhiM            = PhiM,
       lower           = lower,
       upper           = upper)#$par

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


## BOOTSTRAP GAMMA AND RHO ESTIMATIORS
x <- Sys.time()
ps_data <- ps_data_list[[1000]]

bootstrap_dataset_recapture <- function(ps_data, size = NULL){
  
  # Pull relevant data
  recap_f         <- ps_data$recap_f
  recap_m         <- ps_data$recap_m
  apairs_f        <- ps_data$apairs_f
  first_capture_f <- ps_data$first_capture_f   
  nf              <- ps_data$nf
  nm              <- ps_data$nm
  k               <- ps_data$k
  
  # Flag females with pair histories that are viable 
  viable_f <- vector(length = nf)
  
  # History needs to interact with recapture likelihood at least once
  for(i in 1:nf){
    for(t in (first_capture_f[i]+1):(k-1)){
      if(!is.na(apairs_f[i,t]) & !is.na(apairs_f[i,t-1]) & !is.na(apairs_f[i,t+1])){
        if(apairs_f[i,t] != (nm+1)){
          if(apairs_f[i,t] == apairs_f[i,t-1] & apairs_f[i,t] == apairs_f[i,t+1]){
            viable_f[i] <- TRUE
          }
        }
      }
    } 
  }
  
  if(is.null(size)){
    size <- sum(1 * viable_f)
  }
  
  # Construct objects for subsampling 
  recap_f_new <- matrix(NA, nrow = size, ncol = k)
  recap_m_new <- recap_m
  apairs_f_new <- matrix(NA, nrow = size, ncol = k)
  first_capture_f_new <- matrix(NA, nrow = k)
  
  index_viable_f <- which(viable_f)
  
  # sample without replacement and build new recapture dataset
  index_sample                      <- sample(index_viable_f, replace = T, size = size)
  recap_f_new[1:size, 1:k]          <- recap_f[index_sample, 1:k]
  apairs_f_new[1:size, 1:k]         <- apairs_f[index_sample, 1:k]
  first_capture_f_new[1:size]       <- first_capture_f[index_sample]
  
  return(list(recap_f         = recap_f_new,
              recap_m         = recap_m,
              apairs_f        = apairs_f_new,
              first_capture_f = first_capture_f_new,
              nf              = size,
              nm              = nm, 
              k               = k))
}

generate_bootstrap_replicates <- function(ps_data, iter, size = NULL){
  replicates <- lapply(1:iter, function(x) bootstrap_dataset_recapture(ps_data, size))
}

bootstrap_datasets <- generate_bootstrap_replicates(ps_data, 1000)

recapture_correlations_bootstrap <- list()
for(i in 1:length(bootstrap_datasets)){
  
  # Generate One set of Data
  ps_data_boot <- bootstrap_datasets[[i]] 
  
  rl <- compute_jbin_param_cjs(PF,PM)$cor_lower_bound
  ru <- compute_jbin_param_cjs(PF,PM)$cor_upper_bound
  
  lower = c(rl)
  upper = c(ru)
  
  recapture_correlations_bootstrap[[i]] <- nlminb(start          = runif(1, min = lower,max = upper),
                                                  objective       = partial_likelihood,
                                                  PF              = param_list$p.f[1],
                                                  PM              = param_list$p.m[1],
                                                  nf              = ps_data_boot$nf,
                                                  k               = ps_data_boot$k,
                                                  nm              = ps_data_boot$nm,
                                                  first_capture_f = ps_data_boot$first_capture_f,
                                                  recap_f         = ps_data_boot$recap_f,
                                                  recap_m         = ps_data_boot$recap_m,
                                                  apairs_f        = ps_data_boot$apairs_f,
                                                  lower           = lower,
                                                  upper           = upper)$par
}

y <- Sys.time()
difftime(x,y, units = "mins")

plot(density(unlist(recapture_correlations)), main = "Recapture Correlation Estimates", ylim = c(0,10), xlab = "Rho_Hat", ylab = "PDF")
abline(v = mean(unlist(recapture_correlations)), lty = 2, col = "black")
abline(v = param_list$rho[1], col = "green")



lines(density(unlist(recapture_correlations_bootstrap)), main = "Recapture Correlation Estimates", col = "red", xlab = "Rho_Hat", ylab = "PDF")
abline(v = mean(unlist(recapture_correlations_bootstrap)), lty = 2, col = "red")
abline(v = recapture_correlations_bootstrap[[5]], col = "blue")

mean(unlist(recapture_correlations_bootstrap))
sd(unlist(recapture_correlations_bootstrap)) * sqrt(1000)

mean(unlist(recapture_correlations))
sd(unlist(recapture_correlations)) 

quantile(unlist(recapture_correlations_bootstrap), c(0.025,0.25,0.75,0.975))

ps_data <- ps_data_list[[2]]


bootstrap_dataset_survival <- function(ps_data, size = NULL){
  
  # Pull relevant data
  recap_f         <- ps_data$recap_f
  recap_m         <- ps_data$recap_m
  apairs_f        <- ps_data$apairs_f
  first_capture_f <- ps_data$first_capture_f   
  nf              <- ps_data$nf
  nm              <- ps_data$nm
  k               <- ps_data$k
  
  # Flag females with pair histories that are viable 
  viable_f <- vector(length = nf)
  
  for(i in 1:nf){
    for(j in (first_capture_f[i]+1):k){
      if(!is.na(apairs_f[i,j-1])){
        if(apairs_f[i,j-1] == (nm+1)) next
        if(recap_f[i,j-1] == 1 & recap_m[apairs_f[i,j-1],j-1] == 1){
          viable_f[i] <- T
        }
      }
    }
  }
  
  
  if(is.null(size)){
    size <- sum(1 * viable_f)
  }
  
  # Construct objects for subsampling 
  recap_f_new <- matrix(NA, nrow = size, ncol = k)
  recap_m_new <- recap_m
  apairs_f_new <- matrix(NA, nrow = size, ncol = k)
  first_capture_f_new <- matrix(NA, nrow = )
  
  index_viable_f <- which(viable_f)
  
  # sample without replacement and build new recapture dataset
  index_sample                      <- sample(index_viable_f, replace = T, size = size)
  recap_f_new[1:size, 1:k]          <- recap_f[index_sample, 1:k]
  apairs_f_new[1:size, 1:k]         <- apairs_f[index_sample, 1:k]
  first_capture_f_new[1:size]       <- first_capture_f[index_sample]
  
  return(list(recap_f         = recap_f_new,
              recap_m         = recap_m,
              apairs_f        = apairs_f_new,
              first_capture_f = first_capture_f_new,
              nf              = size,
              nm              = nm, 
              k               = k))
}

generate_bootstrap_replicates_surv <- function(ps_data, iter, size = NULL){
  replicates <- lapply(1:iter, function(x) bootstrap_dataset_survival(ps_data, size))
}

bs_replicate_100 <- generate_bootstrap_replicates_surv(ps_data_list[[100]], 1000)


gamma_corr_bootstrap <- lapply(1:length(bs_replicate_100), function(i) compute_correlation(bs_replicate_100[[i]], 
                                                                                             PFM,
                                                                                             0.8,
                                                                                             0.8,
                                                                                             known_pairs = F))


n_bootstrap  <- sapply(1:length(bootstrap_datasets), function(i) gamma_corr_bootstrap[[i]]$n)
N_bootstrap  <- sapply(1:length(bootstrap_datasets), function(i) gamma_corr_bootstrap[[i]]$N)
gam_bootstrap <- sapply(1:length(bs_replicate_100), function(i) gamma_corr_bootstrap[[i]]$obs)

plot(density((gam)))
abline(v = mean(gam))
abline(v = gam_true, col = "red")
lines(density(gam_bootstrap), col = "purple")
abline(v = mean(gam_bootstrap), col = "purple",lty = 2)
abline(v = gam[1], col = "green",lty = 2)
sd(gam)
sd(gam_bootstrap)

p <- gamma_corr[[2]]$n/gamma_corr[[2]]$N

sqrt((p * (1-p))/(PFM[2]^2 * gamma_corr[[2]]$N  * (PhiF * (1-PhiF)) * (PhiM * (1-PhiM))))

