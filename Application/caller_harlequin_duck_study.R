## Load scripts ------------------------------------------------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
src_dir <- getwd()
out_dir <- file.path(src_dir, "Output")
source(file.path(src_dir, "Scripts", "fn_generic.R"))
source(file.path(src_dir, "Scripts", "fn_sim_pair_data.R"))
source(file.path(src_dir, "Scripts", "fn_correlation_estimators.R"))
source(file.path(src_dir, "Scripts", "fn_process_hduck_data.R"))

# Load packages
libs <- c("tidyverse","RMark", "readxl", "lubridate")
load_packages(libs, FALSE)

# Pull/Process Data -------------------------------------------------------------------------------------------------------------
dat_dir <- src_dir %+% "/Data/RE__Harlequin_duck_data/"
cap.data <- prep_cap_data(dat_dir)
ps_data  <- format_cap_data(cap.data)
results <- execute_application(ps_data, small_out = FALSE)

# Save Results
saveRDS(cap.data, out_dir %+% "/hd_cap_data.rds")
saveRDS(ps_data,  out_dir %+%  "/hd_ps_data.rds")
saveRDS(results,  out_dir %+%  "/hd_results.rds")
