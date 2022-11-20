## Load Custom Scripts ---------------------------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
src_dir <- getwd()
out_dir <- getwd() %+% "/Output/"
source(file.path(src_dir, "Scripts", "fn_generic.R"))
source(file.path(src_dir, "Scripts", "fn_sim_pair_data.R"))
source(file.path(src_dir, "Scripts", "fn_correlation_estimators.R"))

# Load packages ---------------------------------------------------------------------------------------------------
libs <- c("tidyverse","RMark")
load_packages(libs, FALSE)
# -----------------------------------------------------------------------------------------------------------------

# Execute Simulation ----------------------------------------------------------------------------------------------
scenario_grid <- get_scenarios()

i <- 20

x <- Sys.time()
results <- execute_simulation(niter      = 3,
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