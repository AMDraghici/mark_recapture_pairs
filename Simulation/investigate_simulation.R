# Aggregate main body simulation
## Load Custom Scripts ---------------------------------------------------------------------------------------------
`%+%`      <- function(a, b) paste0(a, b)
src_dir    <- getwd() 
# Can aggregate multiple simulation studies if desired
out_dir1   <- src_dir %+% YOUR PATH HERE
out_dir2   <- src_dir %+% YOUR PATH HERE
source(file.path(src_dir, "Scripts", "fn_generic.R"))
source(file.path(src_dir, "Scripts", "fn_sim_pair_data.R"))
source(file.path(src_dir, "Scripts", "fn_correlation_estimators.R"))

## Options (ECHO FOR LOGS)
options(echo = TRUE)

# Load packages ---------------------------------------------------------------------------------------------------
libs <- c("tidyverse","marked")
load_packages(libs, FALSE)
# -----------------------------------------------------------------------------------------------------------------

process_results <- function(path){
  # Grab Results
  files <- list.files(path)
  nruns <- length(files)
  
  # Unpack Summaries
  out_list   <- lapply(1:nruns, function(i) readRDS(path %+% files[i]))
  summ_corr  <- do.call(rbind, lapply(1:nruns, function(i) out_list[[i]]$summary_corr))
  summ_cjs   <- do.call(rbind, lapply(1:nruns, function(i) out_list[[i]]$summary_cjs))
  summ_n     <- do.call(rbind, lapply(1:nruns, function(i) out_list[[i]]$summ_n))
  
  summ_corr <- summ_corr  %>% select(-iter) 
  summ_cjs <- summ_cjs %>% select(-iter) 
  summ_n <- summ_n %>% select(-iter)
  
  return(list(summ_corr = summ_corr,
              summ_cjs = summ_cjs,
              summ_n = summ_n))
}

dirs <- c(out_dir1,
          out_dir2)

results_list <- lapply(1:length(dirs), \(i) process_results(dirs[i]))

summ_corr_list <- lapply(1:length(results_list), \(i) results_list[[i]]$summ_corr)
summ_cjs_list <- lapply(1:length(results_list), \(i) results_list[[i]]$summ_cjs)
summ_n_list <- lapply(1:length(results_list), \(i) results_list[[i]]$summ_n)

# Correction from first study
summ_corr_list[[1]] <- summ_corr_list[[1]]  %>% select(-Pval02) %>% rename(Pval02 = Pval03)

summ_corr <- do.call(rbind,summ_corr_list)
summ_cjs  <- do.call(rbind,summ_cjs_list)
summ_n    <- do.call(rbind,summ_n_list)

scenario_grid <- get_scenarios()

summ_corr <-  summ_corr %>% left_join(scenario_grid, by = c("scenario"))
summ_cjs  <-  summ_cjs %>% left_join(scenario_grid, by = c("scenario"))

# Store Results
saveRDS(summ_corr, src_dir)
saveRDS(summ_cjs,  src_dir)
saveRDS(summ_n,    src_dir)