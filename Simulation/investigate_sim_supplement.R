# Aggregate revisions
## Load Custom Scripts ---------------------------------------------------------------------------------------------
`%+%`      <- function(a, b) paste0(a, b)
src_dir    <- getwd() #"/home/mdraghic/projects/def-sbonner/mdraghic/mark_recapture_pair_swap/"
dirs   <- c(src_dir %+% "/Simulation/Run_Post_Thesis/Output/set1/",
            src_dir %+% "/Simulation/Run_Post_Thesis/Output/set2/",
            src_dir %+% "/Simulation/Run_Post_Thesis/Output/set3/",
            src_dir %+% "/Simulation/Run_Post_Thesis/Output/set4/",
            src_dir %+% "/Simulation/Run_Post_Thesis/Output/set5/")
source(file.path(src_dir, "Scripts", "fn_generic2.R"))
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
  keep_bool     <- sapply(1:nruns, function(i) ncol(out_list[[i]]$summary_corr) != 1)
  keep_index    <- 1:length(keep_bool)
  keep_index    <- keep_index[keep_bool]
  summ_corr  <- do.call(rbind, lapply(keep_index, function(i) out_list[[i]]$summary_corr))
  summ_cjs   <- do.call(rbind, lapply(keep_index, function(i) out_list[[i]]$summary_cjs))
  summ_n     <- do.call(rbind, lapply(keep_index, function(i) out_list[[i]]$summ_n))
  
  summ_corr <- summ_corr  %>% select(-iter) 
  summ_cjs <- summ_cjs %>% select(-iter) 
  summ_n <- summ_n %>% select(-iter)
  
  return(list(summ_corr = summ_corr,
              summ_cjs = summ_cjs,
              summ_n = summ_n))
}



results_list <- lapply(1:length(dirs), \(i) process_results(dirs[i]))


scenarios_1 <- sort(unique(results_list[[1]]$summ_corr$scenario))
scenarios_3 <- sort(unique(results_list[[3]]$summ_corr$scenario))


results_list[[1]]$summ_corr <- results_list[[1]]$summ_corr %>% filter((!scenario %in% scenarios_3))
results_list[[1]]$summ_cjs  <- results_list[[1]]$summ_cjs %>% filter((!scenario %in% scenarios_3))
results_list[[1]]$summ_n    <- results_list[[1]]$summ_n %>% filter((!scenario %in% scenarios_3))

results_list[[2]]$summ_corr <- results_list[[2]]$summ_corr %>% filter((!scenario %in% scenarios_1) & (!scenario %in% scenarios_3))
results_list[[2]]$summ_cjs  <- results_list[[2]]$summ_cjs %>% filter((!scenario %in% scenarios_1) & (!scenario %in% scenarios_3))
results_list[[2]]$summ_n    <- results_list[[2]]$summ_n %>% filter((!scenario %in% scenarios_1) & (!scenario %in% scenarios_3))

summ_corr_list <- lapply(1:length(results_list), \(i) results_list[[i]]$summ_corr)
summ_cjs_list  <- lapply(1:length(results_list), \(i) results_list[[i]]$summ_cjs)
summ_n_list    <- lapply(1:length(results_list), \(i) results_list[[i]]$summ_n)

summ_corr <- do.call(rbind,summ_corr_list)
summ_cjs  <- do.call(rbind,summ_cjs_list)
summ_n    <- do.call(rbind,summ_n_list)

scenario_grid <- get_scenarios_extended()

summ_corr <-  summ_corr %>% left_join(scenario_grid, by = c("scenario"))
summ_cjs  <-  summ_cjs %>% left_join(scenario_grid, by = c("scenario"))


saveRDS(summ_corr, src_dir %+% "/Output/revisions/supplement/summ_corr_revisions_final.rds")
saveRDS(summ_cjs,  src_dir %+% "/Output/revisions/supplement/summ_cjs_revisions_final.rds")
saveRDS(summ_n,    src_dir %+% "/Output/revisions/supplement/summ_n_revisions_final.rds")
