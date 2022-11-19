## Load Custom Scripts ---------------------------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
src_dir <- getwd()
out_dir <- getwd()
source(file.path(src_dir, "Scripts", "Production", "fn_generic.R"))
source(file.path(src_dir, "Scripts", "Production", "fn_sim_pair_data.R"))
source(file.path(src_dir, "Scripts", "Production", "fn_correlation_estimators.R"))

# Load packages ---------------------------------------------------------------------------------------------------
libs <- c("tidyverse","RMark")
load_packages(libs, FALSE)
# -----------------------------------------------------------------------------------------------------------------

# Execute Simulation ----------------------------------------------------------------------------------------------
scenario_grid <- get_scenarios()

i <- 56

x <- Sys.time()
results <- execute_simulation(niter      = 1e3,
                              scenario   = scenario_grid[i,"scenario"],
                              PM         = scenario_grid[i,"PM"],
                              PF         = scenario_grid[i,"PF"],
                              PhiF       = scenario_grid[i,"PhiF"],
                              PhiM       = scenario_grid[i,"PhiM"],
                              gam_true   = scenario_grid[i,"gam_true"],
                              rho_true   = scenario_grid[i,"rho_true"],
                              n_pop      = scenario_grid[i,"n_obs"],
                              k          = scenario_grid[i,"k"],
                              init       = NULL)

y <- Sys.time()
difftime(y,x,units = "mins")

# ----------------------------------------------------------------------------------------------------------------

# Return Results -------------------------------------------------------------------------------------------------
saveRDS(results, out_dir %+% "/results_" %+% i %+% ".rds")
# ----------------------------------------------------------------------------------------------------------------

ps_data_list<- readRDS("~/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/Simulation_Data_Gamma4.rds")

PF <- 0.75
PM <- 0.75
PhiF <- 0.8
PhiM <- 0.8
gam_true <- 0.5
rho_true <- 0.5
n_pop <- 300
k <- 10
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
  init         = NULL, # Initial Entry into population for individual n
  show_unmated = T # Include unmated observations in attempt to mate step
)

ps_data_list <- lapply(1:100, function(x) sim_dat(param_list))

cjs_data_list <- lapply(1:length(ps_data_list), function(i) format_to_cjs(ps_data_list[[i]]))
mark_resuts <- lapply(1:length(cjs_data_list), function(i) run_cjs_model_mark(cjs_data_list[[i]],
                                                                              0.80,
                                                                              0.80,
                                                                              0.75,
                                                                              0.75,
                                                                              Iter = i,
                                                                              scenario = 1))

compute_recapture_correlation_simulation
mark_out <- do.call(rbind, mark_resuts)
library(data.table)
setDT(mark_out)
mark_out[,.(coverage = sum(LB <= Truth & Truth <= UB), Bias = mean(Bias)),by = Parameter]

