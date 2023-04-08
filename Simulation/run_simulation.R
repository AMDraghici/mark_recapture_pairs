## Load Custom Scripts ---------------------------------------------------------------------------------------------
`%+%`      <- function(a, b) paste0(a, b)
src_dir    <- "/home/mdraghic/projects/def-sbonner/mdraghic/mark_recapture_pair_swap/"
out_dir    <- src_dir %+% "Simulation/Output/"
Rcpp::sourceCpp(file.path(src_dir, "Src", "generate_pair_data.cpp"))
source(file.path(src_dir, "Scripts", "fn_generic2.R"))
source(file.path(src_dir, "Scripts", "fn_correlation_estimators.R"))
scenario_mapping <- readRDS(src_dir %+% "Simulation/scenario_mapping.rds")

## Options (ECHO FOR LOGS)
options(echo = TRUE)

# Load packages ---------------------------------------------------------------------------------------------------
libs <- c("tidyverse","marked")
load_packages(libs, FALSE)
# -----------------------------------------------------------------------------------------------------------------

# Read command line arguments--------------------------------------------------------------------------------------
args          <- commandArgs(trailingOnly = TRUE)
replicate     <- as.numeric(args[1]) 
scenario      <- scenario_mapping[replicate]
scenario_grid <- get_scenarios()

cat("Running scenario #" %+% scenario_grid[scenario,"scenario"] %+% "...", "\n")
cat("Settings are...","\n")
head(scenario_grid[scenario,])

# Execute Simulation ----------------------------------------------------------------------------------------------
x <- Sys.time()

# Initialize Seed
runif(1)

# Run Simulation
results <- execute_simulation(niter         = 100,
                              bstrp_iter    = 1000,
                              scenario      = scenario_grid[scenario,"scenario"],
                              PM            = scenario_grid[scenario,"PM"],
                              PF            = scenario_grid[scenario,"PF"],
                              PhiF          = scenario_grid[scenario,"PhiF"],
                              PhiM          = scenario_grid[scenario,"PhiM"],
                              gamma         = scenario_grid[scenario,"gam_true"],
                              rho           = scenario_grid[scenario,"rho_true"],
                              n_pop         = scenario_grid[scenario,"n_obs"],
                              k             = scenario_grid[scenario,"k"],
                              Betas         = c(scenario_grid[scenario,"Beta0"], scenario_grid[scenario,"Beta1"]),
                              Delta         = scenario_grid[scenario,"Delta"],
                              PropF         = scenario_grid[scenario,"PropF"],
                              imputed_pairs = scenario_grid[scenario,"imputed_pairs"],
                              tau           = c(rep(1/scenario_grid[scenario,"k"],scenario_grid[scenario,"k"]-1 ),0),
                              small_out     = TRUE)

y <- Sys.time()
difftime(y,x,units = "mins")
# ----------------------------------------------------------------------------------------------------------------

# Return Results -------------------------------------------------------------------------------------------------
saveRDS(results, out_dir %+% "results_" %+% replicate %+% ".rds")
# ----------------------------------------------------------------------------------------------------------------