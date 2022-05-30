# Prep ----------------------------------------------------

# Call R Libs
library(parallel)
library(tidyverse)
library(readxl)
library(lubridate)
library(rjags)

# In-line string concat 
`%+%` <- function(a, b) paste0(a, b)

# Paths
proj_dir <- "/home/mdraghic/projects/def-sbonner/mdraghic/mark_recapture_pair_swap/"
script_dir <- proj_dir %+% "/Scripts/"
out_dir <- proj_dir %+% "/Simulation/"

# Source code
source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "01_fn_model_code.R")

# Simulation Settings -------------------------------------
k = 8
n = 100

# Baseline Parameters that do not vary
to_vary <- list("phi.f" = c(0.8,0.8),
                "phi.m"= c(0.8,0.8),
                "p.f" = c(0.9,0.9),
                "p.m" = c(0.9,0.9),
                "gam" = c(0,0),
                "rho" = c(0,0),
                "n" = c(n)) 
parameter_grid <-list()

ik <- 1
for(i in 1:2){
  for(j in 1:2){
    parameter_grid[[ik]] <- list(n = n,
                                 k = k,
                                 lf = 20,
                                 lm = 20,
                                 prop.female = 0.5,
                                 delta = rep(0.9, k),
                                 phi.f = rep(to_vary$phi.f[i], k),
                                 phi.m = rep(to_vary$phi.m[i], k),
                                 gam = rep(to_vary$gam[i], k),
                                 p.f = rep(to_vary$p.f[i], k),
                                 p.m = rep(to_vary$p.m[i], k),
                                 rho = rep(to_vary$rho[i], k),
                                 betas = list(beta0 = 90, beta1 = 0),
                                 rand_init = F,
                                 init = sample(1, n, TRUE),
                                 show_unmated = ifelse(j == 1, T, F),
                                 data_aug = T)
    ik <- ik+1
  }
}


saveRDS(parameter_grid, out_dir %+% "param_grid.rds")

