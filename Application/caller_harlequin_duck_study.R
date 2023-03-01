## Load scripts ------------------------------------------------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
src_dir <- getwd()
out_dir <- file.path(src_dir, "Output")
source(file.path(src_dir, "Scripts", "fn_generic2.R"))
source(file.path(src_dir, "Scripts", "fn_correlation_estimators.R"))
source(file.path(src_dir, "Scripts", "fn_process_hduck_data.R"))
Rcpp::sourceCpp("Src/generate_pair_data.cpp")

# Load packages
libs <- c("tidyverse", "readxl", "lubridate", "marked")
load_packages(libs, FALSE)

# Pull/Process Data -------------------------------------------------------------------------------------------------------------
dat_dir <- src_dir %+% "/Data/RE__Harlequin_duck_data/"
cap.data <- prep_cap_data(dat_dir)
ps_data  <- format_cap_data(cap.data)
results <- execute_application(ps_data, bstrp_iter = 1e4, small_out = F)

# Save Results
saveRDS(cap.data, out_dir %+% "/hd_cap_data.rds")
saveRDS(ps_data,  out_dir %+%  "/hd_ps_data.rds")
saveRDS(results,  out_dir %+%  "/hd_results.rds")
