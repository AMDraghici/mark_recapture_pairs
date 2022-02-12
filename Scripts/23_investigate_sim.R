library(parallel)
library(dclone)
library(tidyverse)
library(readxl)
library(lubridate)
library(rjags)

setwd("C:/Users/Alex/Documents/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/")
`%+%` <- function(a, b) paste0(a, b)
script_dir <- getwd() %+% "/Scripts/"
dat_dir <- getwd() %+% "/Data/RE__Harlequin_duck_data/"

source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "01_fn_model_code.R")
source(script_dir %+% "02_fn_process_hduck_data.R")
out_dir <- getwd() %+% "/Simulation/Output/"


# Read in parameter data

param_sim_list <- list()
cjs_sim_list <- list()
js_sim_list <- list()
ps_sim_list <- list()
files <- list.files(out_dir)
n.failed <- 0
for(i in 1:100){
  
  # File names 
  param_file <-"parameter_list_" %+% i %+% ".rds"
  cjs_file <- "cjs_run_" %+% i %+% ".rds"
  js_file <- "js_run_" %+% i %+% ".rds"
  ps_file <- "ps_run_" %+% i %+% ".rds"
  
  mask_param <- !(param_file %in% files)
  mask_cjs <-   !(cjs_file %in% files)
  mask_js <-    !(js_file %in% files)
  mask_ps <-    !(ps_file %in% files)
  
  # Check for failed run 
  if(mask_param|mask_cjs|mask_js|mask_ps){
    print("Run " %+% i %+% " failed!")
    n.failed <- n.failed + 1
    next
  }
  
  param_sim_list[[i]] <- readRDS(out_dir %+% param_file)
  cjs_sim_list[[i]]   <- readRDS(out_dir %+% cjs_file)$jags_samples
  js_sim_list[[i]]    <- readRDS(out_dir %+% js_file)$jags_samples
  ps_sim_list[[i]]    <- readRDS(out_dir %+% ps_file)
}


ps_results  <-  process_simulation_data(ps_sim_list,param_sim_list) %>%  mutate(scenario = 1 + (iteration > 25) + (iteration > 50)  + (iteration > 75))
js_results  <-  process_simulation_data(js_sim_list,param_sim_list) %>%  mutate(scenario = 1 + (iteration > 25) + (iteration > 50)  + (iteration > 75))
cjs_results <-  process_simulation_data(cjs_sim_list,param_sim_list) %>% mutate(scenario = 1 + (iteration > 25) + (iteration > 50)  + (iteration > 75))


ps_results %>% 
  group_by(Parameter, scenario) %>% 
  summarize(coverage_50 = mean(In_50),
            coverage_95 = mean(In_95),
            avg_range_50 = mean(Range_50),
            avg_range_95 = mean(Range_95),
            avg_bias = mean(Bias),
            avg_cv = mean(coef_var)) %>% View() #%>%
 # filter(Parameter == "PhiF")


js_results %>% 
  group_by(Parameter, scenario) %>% 
  summarize(coverage_50 = mean(In_50),
            coverage_95 = mean(In_95),
            avg_range_50 = mean(Range_50),
            avg_range_95 = mean(Range_95),
            avg_bias = mean(Bias),
            avg_cv = mean(coef_var)) %>% View()

cjs_results %>% 
  group_by(Parameter, scenario) %>% 
  summarize(coverage_50 = mean(In_50),
            coverage_95 = mean(In_95),
            avg_range_50 = mean(Range_50),
            avg_range_95 = mean(Range_95),
            avg_bias = mean(Bias),
            avg_cv = mean(coef_var)) %>% View()

param_sim_list[[301]]$rho
param_sim_list[[301]]$gam

#   ggplot() + 
#  # geom_hline(yintercept = 0.95, col = "blue") + 
# #  geom_hline(yintercept = 0.5, col = "red") + 
#   geom_point(aes(x = scenario, y = avg_bias), col = "blue") #+
#  # geom_point(aes(x = scenario, y = coverage_50), col = "red")

library(ggmcmc)
cjs_sim_list[[1]] %>% ggs() %>% filter(Parameter %in% c("pM","pF","phiF","phiM")) %>% ggs_traceplot() + ylim(0,1)
js_sim_list[[1]] %>% ggs() %>% filter(Parameter %in% c("pM","pF","phiF","phiM")) %>% ggs_traceplot() + ylim(0,1)
ps_sim_list[[1]] %>% ggs() %>% filter(Parameter %in% c("PM","PF","PhiF","PhiM","rho","gamma", "delta")) %>% ggs_traceplot() + ylim(0,1)
coda.samples %>% ggs() %>% filter(Parameter %in% c("beta0","beta1")) %>% ggs_traceplot() + ylim(-5,5)
