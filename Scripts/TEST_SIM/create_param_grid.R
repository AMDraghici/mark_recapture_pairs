library(parallel)
library(dclone)
library(tidyverse)
library(readxl)
library(lubridate)
library(rjags)

#setwd("C:/Users/Alex/Documents/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/")
`%+%` <- function(a, b) paste0(a, b)
script_dir <- getwd() %+% "/Scripts/"

source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "01_fn_model_code.R")
source(script_dir %+% "02_fn_process_hduck_data.R")
out_dir <- getwd() %+% "/Output/"


#SIM DATA

k = 5 
n = 10

parameter_list <- list(
  n = n,
  k = k,
  prop.female = 0.5,
  delta = rep(0.9, k),
  phi.f = rep(0.7, k),
  phi.m = rep(0.7, k),
  gam = rep(0.6, k),
  p.f = rep(0.8, k),
  p.m = rep(0.8, k),
  rho = rep(0.5, k),
  betas = list(beta0 = 1.0, beta1 = 1.5),
  rand_sex = F,
  rand_init = F,
  init = rep(1,n)
)

# Baseline Parameters that do not vary
to_vary <- list("phi.f" = seq(0.7,0.7,by=0.1),
                "phi.m"= seq(0.7,0.7,by=0.1),
                "p.f" = seq(0.8,0.8,by=0.1),
                "p.m" = seq(0.8,0.8,by=0.1),
                "n" = c(n)) 

# Add parameters that do vary 
to_vary$gam <- seq(0, 1.0, 0.5) #Survival correlation 
to_vary$rho <- seq(0, 1.0, 0.5) # recapture correlation

# Grid of model settings for simulation study 
grid_data <- generate_grid(parameter_list,to_vary)
saveRDS(grid_data$parameter_grid, getwd() %+% "/scripts/TEST_SIM/param_grid.rds")

