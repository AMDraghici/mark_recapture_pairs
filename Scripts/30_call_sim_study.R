library(parallel)
library(dclone)

`%+%` <- function(a, b) paste0(a, b)
script_dir <- getwd() %+% "/Scripts/"

source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "02_fn_model_code.R")
out_dir <- getwd() %+% "/Output/"

k = 4

param_list <- list(
  n = 10, 
  k = k, 
  prop.female = 0.5,
  delta = rep(0.9, k),
  phi.f = rep(0.8, k),
  phi.m = rep(0.8, k),
  gam = rep(1, k),
  p.f = rep(0.8, k),
  p.m = rep(0.8, k),
  rho = rep(0.8, k),
  betas = list(beta0 = 0, beta1 = 100),
  rand_sex = F,
  rand_init = F,
  init = rep(1,10)
)

# # Pull individual dataset
# data <- do.call(simulate_cr_data, param_list)
# cjs_data <- format_to_cjs(data)
# 
# 
# # Multiple Datasets using parallel
# sim_cr_dat(parameter_list = param_list, iterations =  2)


# Run JAGS
set.seed(42)
jags_data <- sim_cr_dat(parameter_list = param_list, iterations =  2)[[1]]

## MCMC parameters  
par_settings <- list('n.iter' = 500, 
                     'n.thin' = 1,
                     'n.burn' = 500,
                     'n.chains' = 2,
                     'n.adapt' = 500) 

## Jags parameters and model script

# Run Full Model + No Groups
jags_params <- c("delta", "PhiF","PhiM","gamma","beta0","beta1")
jags_model <- script_dir %+% "/12_testing.R"

## Run jags in parallel and save results
jags_samples <- run_jags_parallel(jags_data, 
                                  jags_model,
                                  jags_params, 
                                  par_settings,
                                  out_dir,
                                  outname = "TESTING_MODEL")

x <- as.matrix(jags_samples[[1]])[1,]
x <- as.data.frame(x)
rownames_to_column(x) %>% filter(rowname %in% paste0("apf[",as.character(1:10),",2]"))
