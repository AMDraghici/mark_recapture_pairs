## Load scripts ------------------------------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
src_dir <- getwd()
source(file.path(src_dir, "Scripts", "Production", "fn_generic.R"))
source(file.path(src_dir, "Scripts", "Production", "fn_sim_pair_data.R"))
source(file.path(src_dir, "Scripts", "Production", "fn_correlation_estimators.R"))

# Load packages
libs <- c("tidyverse","RMark")
load_packages(libs, FALSE)

# Run Simualtion
PM       <- 0.85
PF       <- 0.85
PhiF     <- 0.8
PhiM     <- 0.8
gam_true <- 0.25
rho_true <- 0.4
n_pop    <- 350
k        <- 30

x <- Sys.time()
results <- execute_simulation(niter      = 10,
                              PM         = PM,
                              PF         = PF,
                              PhiF       = PhiF,
                              PhiM       = PhiM,
                              gam_true   = gam_true,
                              rho_true   = rho_true,
                              n_pop      = n_pop,
                              k          = k,
                              init       = NULL)

y <- Sys.time()
difftime(y,x,units = "mins")