# Load Libraries and Custom Fns

library(parallel)
library(dclone)
library(tidyverse)
library(readxl)
library(lubridate)
library(rjags)
library(nimble)

setwd("C:/Users/Alex/Documents/Projects/Research/Chapter 2 - Dyads/Code/mark_recapture_pair_swap/")
`%+%` <- function(a, b) paste0(a, b)
script_dir <- getwd() %+% "/Scripts/"
dat_dir <- getwd() %+% "/Data/RE__Harlequin_duck_data/"

source(script_dir %+% "00_fn_sim_pair_data.R")
source(script_dir %+% "01_fn_model_code.R")
source(script_dir %+% "02_fn_process_hduck_data.R")
source(script_dir %+% "12_pair_swap_mod_nimble.R")
source(script_dir %+% "11_jolly_seber_mod_nimble.R")
out_dir <- getwd() %+% "/Output/"

# TESTING SIMULATED DATA METHOD -------------------------------------------------------------------------------------------------
# Set number of occasions and animals
k = 5
n = 20

# Seeds for Testing
# set.seed(42)
# set.seed(1e4)
# set.seed(4)
set.seed(1e5)

# Parameter Grid 
param_list <- list(
  n            = n, # Number of Animals 
  k            = k, # Occasions 
  lf           = 10, # Data Augmentation for Females (M_F)
  lm           = 10, # Data Augmentation for Males (M_M)
  prop.female  = 0.5, # Proportion of simulated individuals to be female 
  delta        = rep(0.9, k), # Probability that mating is attempted 
  phi.f        = rep(0.8, k), # Marginal Prob of Female Survival
  phi.m        = rep(0.8, k), # Marginal Prob of Male Survival
  gam          = rep(0.2, k), # Correlation in Survival Prob of Mates
  p.f          = rep(0.9, k), # Marginal Prob of Female Recapture
  p.m          = rep(0.9, k), # Marginal Prob of Male Recapture 
  rho          = rep(0.6, k), # Correlation in male survival rates 
  betas        = list(beta0 = 90, beta1 = 0), # inv.logit(Beta0 + Beta1 * hij) = Prob of reforming a pair from t-1 after hij times together
  rand_init    = F, # Randomize Initial Entry (just leave as F) 
  init         = sample(1, n, TRUE), # Initial Entry into population for individual n
  show_unmated = T, # Include unmated observations in attempt to mate step
  data_aug     = T  # Add individuals to data augmentation 
)

# Generate One set of Data
ps_data <- sim_dat(param_list) # pair-swap data
js_data <- format_to_js(ps_data) # jolly-seber data 

# Parameters of Interest for Nimble 
nimble_params <- c("PF","PM","rho",
                   "rho_raw","gamma_raw",
                   "PhiF","PhiM","gamma",
                   "delta","beta0","beta1",
                   "eps", "gl", "gu", "ru", "rl", 
                   "NF", "NM")

# PAIR SWAP 
# Compile Nimble Model for Pair-Swap 
start <- Sys.time()
CpsMCMC_List <- compile_pair_swap_nimble(ps_data, nimble_params)
end <- Sys.time()
print(start-end)

# Run Pair-Swap Model 
start <- Sys.time()
samples <- run_nimble(CpsMCMC_List$CpsMCMC,niter = 1e5,nburnin = 5e4, thin = 1)
end <- Sys.time()
print(start-end)

# Check Summary and Marginals (PF,PM,Phif,PhiM)
summary(samples)
plot_caterpillar(gather_posterior_summary(samples)) + ylim(c(0,1.0))


# JOLLY - SEBER 

# Compile nimble for Jolly-Seber Model 
start <- Sys.time()
CjsMCMC_List <- compile_jolly_seber_nimble(js_data)
end <- Sys.time()
print(start-end)

# Run Jolly-Seber Model 
start <- Sys.time()
samples2 <- run_nimble(CjsMCMC_List$CjsMCMC,niter = 1e5,nburnin = 1e4, thin = 10)
end <- Sys.time()
print(start-end)

summary(samples2)
plot_caterpillar(gather_posterior_summary(samples2))  + ylim(c(0.5,1.0))

# Look at Traceplots and Densities
library(ggmcmc)
samples %>% ggs() %>% filter(Parameter %in% c("beta0")) %>% ggs_traceplot() #+ ylim(0,1)
samples %>% ggs() %>% filter(Parameter %in% c("PhiF","PhiM","PF","PM")) %>% ggs_traceplot()#+ ylim(0.01,0.99)
samples %>% ggs() %>% filter(Parameter %in% c("delta")) %>% ggs_traceplot() #+ ylim(0,1)
samples %>% ggs() %>% filter(Parameter %in% c("eps[" %+% 1:ps_data$k %+% "]")) %>% ggs_traceplot() + ylim(0,1)
samples %>% ggs() %>% filter(Parameter %in% c("gamma","rho")) %>% ggs_traceplot() #+ ylim(-1,1)
samples %>% ggs() %>% filter(Parameter %in% c("gamma_raw","rho_raw")) %>% ggs_traceplot() #+ ylim(0,1)
samples %>% ggs() %>% filter(Parameter %in% c("xi")) %>% ggs_traceplot() #+ ylim(0,1)




















# HDUCK EXPERIMENTS-----------------------------------------------------------------------------------------------------------

n <- 314
k <- 28
param_list <- list(
  n = n,
  k = k,
  lf = 20,
  lm = 20,
  prop.female = 0.4652568,
  delta = rep(0.55, k),
  phi.f = rep(0.8, k),
  phi.m = rep(0.75, k),
  gam = rep(0.5, k),
  p.f = rep(0.5, k),
  p.m = rep(0.4, k),
  rho = rep(0.75, k),
  betas = list(beta0 = 1, beta1 = 1.5),
  rand_init = F,
  init = sample(k-1, n, TRUE),
  show_unmated = T,
  data_aug = T
)



DIM <- function(x){return(c(NROW(x),NCOL(x)))}

x <- lapply(1:length(jags_data), function(x) DIM(jags_data[[x]]))
names(x) <- names(jags_data)

#HDUCK Data
cap.data <- gather_hq_data(dat_dir) %>% build_cr_df() %>% populate_missing_mate_data() %>%
  populate_missing_mate_data() %>% 
  add_implied_states() %>%
  add_last_capture() %>% 
  clean_filtered() 

drop_id_yr_filter <- cap.data %>%
  group_by(animal_id) %>%
  summarize(num = sum(recapture_individual)) %>% 
  filter(num < 1) %>%
  pull(animal_id)



cap.data <- cap.data %>%
  # filter(!(animal_id %in% c(24,drop_id_yr_filter)), time > 1, initial_entry < 11) %>% 
  filter(!(animal_id %in% c(drop_id_yr_filter))) %>% 
  assign_ids_bysex()
# 
# cap.data <- cap.data %>% mutate(time = time - 1,
#                     initial_entry = initial_entry-1)

cap.data <- cap.data %>% 
  # filter(initial_entry < 11 & time < 11) %>% 
  filter(initial_entry < 28)


jags_data <- build_jags_data(cap.data, data_aug = T, lf = 20, lm = 20)
cjs_data <- format_to_cjs(jags_data)
js_data <- format_to_js(jags_data)


nimble_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps", "gl", "gu", "ru", "rl", "NF", "NM")
start <- Sys.time()
CpsMCMC_List <- compile_pair_swap_nimble(jags_data, nimble_params)
end <- Sys.time()
print(start-end)

start <- Sys.time()
samples <- run_nimble(CpsMCMC_List$CpsMCMC,niter = 500,nburnin = 1, thin = 1)
end <- Sys.time()
print(start-end)

# jid 138
# aid 558
# pid 453

y <- lapply(1:length(jags_data2), function(y) DIM(jags_data2[[y]]))
names(y) <- names(jags_data2)

# # 
# # animal1cap <- cap.data %>%
# #   group_by(animal_id) %>%
# #   summarize(init = min(initial_entry),
# #             first = min(which(recapture_individual==1)),
# #             last = max(which(recapture_individual==1))) %>%
# #   ungroup() %>%
# #   mutate(known_lifespan = last-first + 1) %>%
# #   filter(known_lifespan <= 1) %>%
# #   pull(animal_id)
# 
# # Drop transients and assume unmated when dropped (NEED TO FIX THIS)
# cap.data <- cap.data %>%
# #   filter(!(animal_id %in% animal1cap)) %>%
# #   mutate(mated = ifelse(partner_id %in% animal1cap,0,mated),
# #          partner_id = ifelse(partner_id %in% animal1cap, 0, partner_id)) %>%
#   populate_missing_mate_data() %>% 
#   #filter(initial_entry <= 10, time <= 10) %>% 
#   add_implied_states() %>%
#   add_last_capture() %>% 
#   clean_filtered() 
# 
# 
# 
# drop_id_yr_filter <- cap.data %>% group_by(animal_id) %>% summarize(num = sum(recapture_individual)) %>% filter(num < 1) %>% pull(animal_id)
# 
# cap.data <- cap.data %>% filter(!(animal_id %in% drop_id_yr_filter)) %>% 
#   assign_ids_bysex()
# 
# 
# # Age diff experiment
# cap.data %>% select(animal_id, time, lower_age, upper_age, initial_entry) %>% rename(partner_id = animal_id,
#                                                                                      lower_age_partner = lower_age,
#                                                                                      upper_age_partner = upper_age,
#                                                                                      initial_entry_partner = initial_entry) -> x
# 
# 
# 
# 
# cap.data <- left_join(cap.data, x, by = c("partner_id", "time"))
# 
# 
# cap.data %>% filter(!is.na(partner_id)) %>% filter(partner_id != 0) %>% select(animal_id, partner_id, time, lower_age, upper_age, lower_age_partner, upper_age_partner)
# 
# y <- cap.data %>% filter(!is.na(partner_id)) %>% filter(partner_id != 0, time >= initial_entry & time>=initial_entry_partner) %>% 
#   mutate(max_age_diff = pmax(upper_age_partner - lower_age, upper_age - lower_age_partner),
#          probable_age_diff = abs((time-initial_entry) - (time-initial_entry_partner)),
#          min_age_diff = pmin(abs(upper_age_partner - lower_age), abs(upper_age_partner - upper_age),
#                              abs(lower_age_partner - upper_age), abs(lower_age_partner - lower_age)))
# 
# barplot(prop.table(table(y$max_age_diff)))
# barplot(prop.table(table(y$probable_age_diff)))
# barplot(prop.table(table(y$min_age_diff)))
# 
# jags_data <- build_jags_data(cap.data)
# js_data <- format_to_js(jags_data)
# cjs_data <- format_to_cjs(jags_data)
# 
# # Investigate max distance between recaptures before attrition
# max(sapply(1:jags_data$nf, function(i) max(diff(which(jags_data$recap_f[i,]==1)))))
# max(sapply(1:jags_data$nm, function(i) max(diff(which(jags_data$recap_m[i,]==1)))))
# 
# # Investigate overall max distance between recaptures 
# max(sapply(1:jags_data$nf, function(i) diff(range(which(jags_data$recap_f[i,]==1)))))
# max(sapply(1:jags_data$nm, function(i) diff(range(which(jags_data$recap_m[i,]==1)))))
# 
# animals_left <- (1:615)[!1:615 %in% animal1cap]
# 
# cap.data %>% filter(animal_id == 28) %>% select(partner_id, time, mated, recapture_individual)
# cap.data %>% filter(animal_id == 116) %>% select(partner_id, time, mated, recapture_individual)
# # 
# # [1]   1   2   3   4   7   8   9  10  11  14  17  18  19  20
# # [15]  21  23  24  28  29  30  34  49  52  53  59  61  65  66
# # [29]  74  79  80  81  84  86  87  89  99 100 101 102 105 110
# # [43] 112 113 115 116 120 121 125 126 127 130 133 149 151 152
# # [57] 155 157 161 163 166 168 172 173 176 180 189 190 194 197
# # [71] 198 199 205 213 220 224 228 231 232 240 242 249 253 255
# # [85] 256 259 268 281 286 292 294 319 320 324 327 331 333 334
# # [99] 335 336 337 338 339 340 341 342 347 356 357 358 359 360
# # [113] 361 382 396 411 421 430 449 451 453 454 455 457 459 462
# # [127] 464 465 466 470 471 480 482 483 485 486 489 490 507 513
# # [141] 514 531 534 537 540 543 553 554 557 558 559 561 563 565
# # [155] 575 576 578 579 580 585 587 615
# 
# # for(i in 1:length(jags_data)){
# #   print(names(jags_data[i]))
# #   if(is.null(dim(jags_data[[i]]))){
# #     print(jags_data[[i]])
# #   } else {
# #     print(dim(jags_data[[i]]))
# #   }
# # }
# 
# for(t in 1:30){
#   index <- which(!is.na(jags_data$af[,t]))
#   print(all(jags_data$sf[index,t] == jags_data$af[index,t]))
# } 
