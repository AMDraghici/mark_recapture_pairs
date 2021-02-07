`%+%` <- function(a, b) paste0(a, b)
script_dir <- getwd() %+% "/Scripts/"

source(script_dir %+% "00_fn_sim_pair_data.R")

k = 4

parameter_list <- list(
  n = 10, 
  k = k, 
  prop.female = 0.5,
  delta = rep(0.9, k),
  phi.f = rep(0.8, k),
  phi.m = rep(0.8, k),
  gam = rep(0.5, k),
  p.f = rep(0.5, k),
  p.m = rep(0.5, k),
  rho = rep(0.8, k),
  betas = list(beta0 = 0, beta1 = 100),
  rand_sex = F,
  rand_init = F,
  init = rep(1,50)
)

data <- do.call(simulate_cr_data, parameter_list)
cjs_data <- format_to_cjs(data)


