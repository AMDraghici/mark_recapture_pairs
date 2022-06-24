# Set Up Environment
library(parallel)
library(dclone)
library(tidyverse)
library(readxl)
library(lubridate)
library(rjags)
library(ggmcmc)

setwd("YOUR PATH HERE")
`%+%` <- function(a, b) paste0(a, b)
out_dir <- getwd() %+% "/no_repartner_no_corr_good/"

# Grab Data
i <- 102
js_run_i <- readRDS(out_dir %+% "js_run_" %+% i %+% ".rds")
ps_run_i <- readRDS(out_dir %+% "ps_run_" %+% i %+% ".rds")
param_list <- readRDS(out_dir %+% "parameter_list_" %+% i %+% ".rds")

# Summaries
summary(js_run_i$samples)
summary(ps_run_i$samples)

# Some plots

# Recapture and Survival Values (JS and PS)
js_run_i$samples %>% ggs() %>% 
  filter(Parameter %in% c("PM","PF")) %>% ggs_traceplot() + geom_hline(yintercept = 0.9, col = "red")
js_run_i$samples %>% ggs() %>% 
  filter(Parameter %in% c("PhiF","PhiM")) %>% ggs_traceplot() + geom_hline(yintercept = 0.8, col = "red")

ps_run_i$samples %>% ggs() %>%
  filter(Parameter %in% c("PM","PF")) %>% ggs_traceplot() + geom_hline(yintercept = 0.9, col = "red")
ps_run_i$samples %>% ggs() %>% 
  filter(Parameter %in% c("PhiF","PhiM")) %>% ggs_traceplot() + geom_hline(yintercept = 0.8, col = "red")

# Mating Prob (desired to mate)
ps_run_i$samples %>% ggs() %>% filter(Parameter %in% c("delta")) %>% ggs_traceplot() + geom_hline(yintercept = 0.9, col = "red")

# gamma (survival corr) and rho (recapture corr) only applies to scenario 1
ps_run_i$samples %>%
  ggs() %>% filter(Parameter %in% c("rho","gamma")) %>% ggs_traceplot() + geom_hline(yintercept = 0, col = "red")

