# Prepare Environment -------------------------------------------------------------------------------------------------------------------------
#Libraries 
libs <- c("parallel", "dclone", "tidyverse", "readxl", "lubridate", "rjags")

## If package is not installed then do so
for (pkg in libs) {
  if (!pkg %in% rownames(installed.packages())) {
    install.packages(pkg)
  }
}

## Attach our libraries
lapply(libs, require, character.only = TRUE)

# Get directories 
`%+%` <- function(a, b) paste0(a, b)
script_dir <- getwd() %+% "/Scripts/"
dat_dir <- getwd() %+% "/Data/RE__Harlequin_duck_data/"
out_dir <- getwd() %+% "/Output/"

# Source in custom functions 
source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "02_fn_model_code.R")
source(script_dir %+% "03_fn_process_hduck_data.R")

# Prepare Data Settings -------------------------------------------------------------------------------------------------------------------------------------

# Simulation Study data settings 
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

# Prepare Model Settings --------------------------------------------------------------------------------------------------------------------------

par_settings <- list('n.iter' = 1e1, 
                     'n.thin' = 1e0,
                     'n.burn' = 1e1,
                     'n.chains' = 1,
                     'n.adapt' = 1e1)

jags_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps")
jags_model <- script_dir %+% "/11_mod_pair_swap_notime.R"

# Run the simulation ------------------------------------------------------------------------------------------------------------------------------

# Iterations per setting
iterations <- 2
# Cores to use
ncores <- 3 

run_simulation(grid_data = grid_data,
               jags_model = jags_model, 
               jags_params = jags_params,
               par_settings = par_settings,
               out_dir = out_dir,
               iterations = iterations,
               ncores = ncores)