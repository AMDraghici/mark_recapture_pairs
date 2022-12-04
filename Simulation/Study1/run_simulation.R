## Load Custom Scripts ---------------------------------------------------------------------------------------------
`%+%`      <- function(a, b) paste0(a, b)
src_dir    <- "/home/mdraghic/projects/def-sbonner/mdraghic/mark_recapture_pair_swap/"
out_dir    <- src_dir %+% "Simulation/Study1/Output/"
source(file.path(src_dir, "Scripts", "fn_generic.R"))
source(file.path(src_dir, "Scripts", "fn_sim_pair_data.R"))
source(file.path(src_dir, "Scripts", "fn_correlation_estimators.R"))

## Options (ECHO FOR LOGS)
options(echo = TRUE)

# Load packages ---------------------------------------------------------------------------------------------------
libs <- c("tidyverse","RMark")
load_packages(libs, FALSE)
# -----------------------------------------------------------------------------------------------------------------

# Read command line arguments--------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
scenario <- as.numeric(args[1])
scenario_grid <- get_scenarios()

cat("Running scenario #" %+% scenario_grid[scenario,"scenario"] %+% "...", "\n")
cat("Settings are...","\n")
head(scenario_grid[scenario,])

# Execute Simulation ----------------------------------------------------------------------------------------------
x <- Sys.time()

# Initialize Seed
runif(1)

# Run Simulation
results <- execute_simulation(niter      = 1e3,
                              scenario   = scenario_grid[scenario,"scenario"],
                              PM         = scenario_grid[scenario,"PM"],
                              PF         = scenario_grid[scenario,"PF"],
                              PhiF       = scenario_grid[scenario,"PhiF"],
                              PhiM       = scenario_grid[scenario,"PhiM"],
                              gam_true   = scenario_grid[scenario,"gam_true"],
                              rho_true   = scenario_grid[scenario,"rho_true"],
                              n_pop      = scenario_grid[scenario,"n_obs"],
                              k          = scenario_grid[scenario,"k"],
                              init       = NULL)

y <- Sys.time()
difftime(y,x,units = "mins")
# ----------------------------------------------------------------------------------------------------------------

# Return Results -------------------------------------------------------------------------------------------------
saveRDS(results, out_dir %+% "results_" %+% scenario %+% ".rds")
# ----------------------------------------------------------------------------------------------------------------