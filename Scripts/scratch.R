###Setup###
rm(list=ls())
options(scipen = 999)
##Libraries
libs <- c("tidyverse","openxlsx","data.table")

##If package is not installed then do so
for(pkg in libs){
  if(!pkg %in% rownames(installed.packages()))
    install.packages(pkg)
}

##Attach our libraries
lapply(libs,require,character.only=TRUE)

##Inline object concatenation (like python)
`%+%` <- function(a,b) paste0(a,b)

setwd("/home/alexd/Documents/Research/Chapter 2 - Dyads/")
source("/home/alexd/Documents/Research/Chapter 2 - Dyads/Code/Pair_Swap_Model/Scripts/SimulateData2.R")
##Set Working Directory and reference paths
path2dat <-getwd() %+% "/Data/"
path2out <- getwd() %+% "/Output/"
path2scripts <- getwd() %+% "/Scripts/"

###Simulate gamma-rho-CJS Data ###
# ##Sample Size (individuals - will be split into groups below)
# n <- 100
# ##Sampling Occasions
# k <- 15
# prop.female = 0.5 
# n.covariates = 2 
# beta.m = c(0.5,3)
# beta.f = c(0.5,3)
# effect.name.m = c("Intercept","History")
# effect.name.f = c("Intercept","History")
# phi.f <- rep(0.8,k)
# phi.m <- rep(0.8,k)
# gam <- rep(0.5,k)
# p.f <- rep(0.8,k)
# p.m <- rep(0.8,k)
# rho <- rep(0.5,k)
# delta <- rep(0.9,k)
k <- 5
input_list <- list(n =100,
                   k = 5,
                   prop.female = 0.5,
                   n.covariates = 2, 
                   beta.m = c(0.5,3),
                   beta.f = c(0.5,3),
                   effect.name.m = c("Intercept","History"),
                   effect.name.f = c("Intercept","History"),
                   phi.f =rep(0.8,k),
                   phi.m =rep(0.8,k),
                   gam =rep(0.5,k),
                   p.f =rep(0.8,k),
                   p.m =rep(0.8,k),
                   rho =rep(0.5,k),
                   delta =rep(0.9,k))

data <- do.call(simulate_cr_data, input_list)

system.time(do.call(simulate_cr_data, input_list))
