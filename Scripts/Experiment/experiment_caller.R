library(parallel)
library(dclone)
library(tidyverse)
library(readxl)
library(lubridate)

`%+%` <- function(a, b) paste0(a, b)
script_dir <- getwd() %+% "/Scripts/"
source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "02_fn_model_code.R")
out_dir <- getwd() %+% "/Output/"

n <- 1000
k <- 28

psi <- array(1, dim = c(n,n,k))

for(t in 1:k){
  for(i in 1:n){
    psi[i,1:n,t] <- 1/n 
  }
}

pairs <- array(NA, dim = c(n,n,k))

for(t in 1:k){
  for(i in 1:n){
    if(i == 1){
      pairs[i,1:n,t] <- rmultinom(1,1,psi[i,1:n,t])
    } else if(i == 2){
      pairs[i,1:n,t] <- rmultinom(1,1,psi[i,1:n,t]*(1-pairs[1,,t]))
    } else {
      pairs[i,1:n,t] <- rmultinom(1,1,psi[i,1:n,t]*(1-colSums(pairs[1:(i-1),,t])))
    }
    
  }
}

pairs2 <- apply(pairs,1,function(x) which(x == 1))
pairs2 <- t(pairs2)

jags_data <- list(pairs = pairs,
                  pairs2 = pairs2,
                  psi = psi,
                  n = n,
                  k = k)

jags_params <- c("psi2")
jags_model1 <- script_dir %+% "experiment_binom.R"
jags_model2<- script_dir %+% "experiment_multinom.R"
jags_model3<- script_dir %+% "Experiment/experiment_cat.R"

par_settings <- list('n.iter' = 1e2, 
                     'n.thin' = 1,
                     'n.burn' = 1e2,
                     'n.chains' = 2,
                     'n.adapt' = 1e2)
# 
# jags_samples1 <- run_jags_parallel(jags_data, 
#                                    jags_model1,
#                                    jags_params, 
#                                    par_settings,
#                                    out_dir,
#                                    outname = "experiment_binom")


jags_samples2 <- run_jags_parallel(jags_data, 
                                   jags_model3,
                                   jags_params, 
                                   par_settings,
                                   out_dir,
                                   outname = "experiment_cat")
