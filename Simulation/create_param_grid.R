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
out_dir <- proj_dir %+% "Simulation/"

# Source code
source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "01_fn_model_code.R")

# Simulation Settings -------------------------------------

k = 15 
n = 200

#set.seed(42)
parameter_list <- list(
  n = n,
  k = k,
  prop.female = 0.5,
  delta = rep(0.9, k),
  phi.f = rep(0.9, k),
  phi.m = rep(0.9, k),
  gam = rep(0.7, k),
  p.f = rep(0.8, k),
  p.m = rep(0.8, k),
  rho = rep(0.7, k),
  betas = list(beta0 = 0.5, beta1 = 1.0),
  rand_init = F,
  init = sample(x = k-1, size = n, replace = T)
)

# Baseline Parameters that do not vary
to_vary <- list("phi.f" = seq(0.9,0.9,by=0.1),
                "phi.m"= seq(0.9,0.9,by=0.1),
                "p.f" = seq(0.8,0.8,by=0.1),
                "p.m" = seq(0.8,0.8,by=0.1),
                "n" = c(n)) 

# Add parameters that do vary 
to_vary$gam <- c(0.2, 0.7) #Survival correlation 
to_vary$rho <- c(0.2, 0.7)  # recapture correlation

# Grid of model settings for simulation study 
grid_data <- generate_grid(parameter_list,to_vary)
saveRDS(grid_data$parameter_grid, out_dir %+% "param_grid.rds")

