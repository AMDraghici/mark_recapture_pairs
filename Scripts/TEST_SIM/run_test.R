library(parallel)
library(dclone)
library(tidyverse)
library(readxl)
library(lubridate)
library(rjags)


#setwd("C:/Users/Alex/Documents/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/")
`%+%` <- function(a, b) paste0(a, b)
script_dir <- getwd() %+% "/Scripts/"

source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "01_fn_model_code.R")
source(script_dir %+% "02_fn_process_hduck_data.R")
out_dir <- getwd() %+% "/Output/"

## Options
options(echo = TRUE)

## Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
k <- args[1]
pars_mat_name <- args[2]

cat(k,"\n")
cat(pars_mat_name,"\n")

## Load parameter matrix
pars_list <- readRDS(pars_mat_name)

## Run replicate
# jags_data <- sim_dat(pars_list[[k]])
jags_data_list <- sim_cr_dat(pars_list[[k]], iterations = 5, ncores = 5)