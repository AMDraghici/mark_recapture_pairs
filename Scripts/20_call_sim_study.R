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
out_dir <- getwd() %+% "/Output/"


k = 10
n = 100

#set.seed(42)
param_list <- list(
  n = n,
  k = k,
  prop.female = 0.5,#0.4652568,
  delta = rep(0.9, k),
  phi.f = rep(0.8, k),
  phi.m = rep(0.8, k),
  gam = rep(0.7, k),
  p.f = rep(0.9, k),
  p.m = rep(0.9, k),
  rho = rep(0.7, k),
  betas = list(beta0 = 1, beta1 = 1.5),
  rand_init = F,
  init = sample(1, n, TRUE),
  show_unmated = T
)


jags_data <- sim_dat(param_list)



CpsMCMC <- compile_pair_swap_nimble(jags_data)
samples <- run_nimble(CpsMCMC,niter = 2e4,nburnin = 1e4, thin = 1)

gather_posterior_summary(samples)
plot_caterpillar(samples)

# 
# saveRDS(samples, "long_run_332_28.rds")
library(ggmcmc)
samples %>% ggs() %>% filter(Parameter %in% c("beta0","beta1")) %>% ggs_traceplot() + ylim(-5,5)
samples %>% ggs() %>% filter(Parameter %in% c("PhiF","PhiM","PF","PM")) %>% ggs_traceplot() #+ ylim(0.70,0.95)
samples %>% ggs() %>% filter(Parameter %in% c("delta")) %>% ggs_traceplot() + ylim(0,1)
samples %>% ggs() %>% filter(Parameter %in% c("eps[" %+% 1:jags_data$k %+% "]")) %>% ggs_traceplot() + ylim(0,1)
samples %>% ggs() %>% filter(Parameter %in% c("gamma","rho")) %>% ggs_traceplot() + ylim(-1,1)
samples %>% ggs() %>% filter(Parameter %in% c("gamma_kappa_raw","rho_kappa_raw")) %>% ggs_traceplot() + ylim(0,10)
samples %>% ggs() %>% filter(Parameter %in% c("gamma_phi_raw","rho_phi_raw")) %>% ggs_traceplot() + ylim(0,1)
samples %>% ggs() %>% filter(Parameter %in% c("gamma_raw","rho_raw")) %>% ggs_traceplot() + ylim(0,1)

n <- 332
k <- 28
param_list <- list(
  n = n,
  k = k,
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
  show_unmated = T
)


#HDUCK Data
cap.data <- gather_hq_data(dat_dir) %>% 
  build_cr_df() %>% 
  populate_missing_mate_data() %>% 
  populate_missing_mate_data() %>%
  add_implied_states() %>%
  add_last_capture() %>% 
  clean_filtered() %>% 
  assign_ids_bysex()


cap.data <- cap.data %>% filter(initial_entry < 28)
jags_data <- build_jags_data(cap.data)
cjs_data <- format_to_cjs(jags_data)


cap.data <- gather_hq_data(dat_dir) %>% build_cr_df() %>% populate_missing_mate_data() 
# 
# animal1cap <- cap.data %>%
#   group_by(animal_id) %>%
#   summarize(init = min(initial_entry),
#             first = min(which(recapture_individual==1)),
#             last = max(which(recapture_individual==1))) %>%
#   ungroup() %>%
#   mutate(known_lifespan = last-first + 1) %>%
#   filter(known_lifespan <= 1) %>%
#   pull(animal_id)

# Drop transients and assume unmated when dropped (NEED TO FIX THIS)
cap.data <- cap.data %>%
#   filter(!(animal_id %in% animal1cap)) %>%
#   mutate(mated = ifelse(partner_id %in% animal1cap,0,mated),
#          partner_id = ifelse(partner_id %in% animal1cap, 0, partner_id)) %>%
  populate_missing_mate_data() %>% 
  #filter(initial_entry <= 10, time <= 10) %>% 
  add_implied_states() %>%
  add_last_capture() %>% 
  clean_filtered() 



drop_id_yr_filter <- cap.data %>% group_by(animal_id) %>% summarize(num = sum(recapture_individual)) %>% filter(num < 1) %>% pull(animal_id)

cap.data <- cap.data %>% filter(!(animal_id %in% drop_id_yr_filter)) %>% 
  assign_ids_bysex()


# Age diff experiment
cap.data %>% select(animal_id, time, lower_age, upper_age, initial_entry) %>% rename(partner_id = animal_id,
                                                                                     lower_age_partner = lower_age,
                                                                                     upper_age_partner = upper_age,
                                                                                     initial_entry_partner = initial_entry) -> x




cap.data <- left_join(cap.data, x, by = c("partner_id", "time"))


cap.data %>% filter(!is.na(partner_id)) %>% filter(partner_id != 0) %>% select(animal_id, partner_id, time, lower_age, upper_age, lower_age_partner, upper_age_partner)

y <- cap.data %>% filter(!is.na(partner_id)) %>% filter(partner_id != 0, time >= initial_entry & time>=initial_entry_partner) %>% 
  mutate(max_age_diff = pmax(upper_age_partner - lower_age, upper_age - lower_age_partner),
         probable_age_diff = abs((time-initial_entry) - (time-initial_entry_partner)),
         min_age_diff = pmin(abs(upper_age_partner - lower_age), abs(upper_age_partner - upper_age),
                             abs(lower_age_partner - upper_age), abs(lower_age_partner - lower_age)))

barplot(prop.table(table(y$max_age_diff)))
barplot(prop.table(table(y$probable_age_diff)))
barplot(prop.table(table(y$min_age_diff)))

jags_data <- build_jags_data(cap.data)
js_data <- format_to_js(jags_data)
cjs_data <- format_to_cjs(jags_data)

# Investigate max distance between recaptures before attrition
max(sapply(1:jags_data$nf, function(i) max(diff(which(jags_data$recap_f[i,]==1)))))
max(sapply(1:jags_data$nm, function(i) max(diff(which(jags_data$recap_m[i,]==1)))))

# Investigate overall max distance between recaptures 
max(sapply(1:jags_data$nf, function(i) diff(range(which(jags_data$recap_f[i,]==1)))))
max(sapply(1:jags_data$nm, function(i) diff(range(which(jags_data$recap_m[i,]==1)))))

animals_left <- (1:615)[!1:615 %in% animal1cap]

cap.data %>% filter(animal_id == 28) %>% select(partner_id, time, mated, recapture_individual)
cap.data %>% filter(animal_id == 116) %>% select(partner_id, time, mated, recapture_individual)
# 
# [1]   1   2   3   4   7   8   9  10  11  14  17  18  19  20
# [15]  21  23  24  28  29  30  34  49  52  53  59  61  65  66
# [29]  74  79  80  81  84  86  87  89  99 100 101 102 105 110
# [43] 112 113 115 116 120 121 125 126 127 130 133 149 151 152
# [57] 155 157 161 163 166 168 172 173 176 180 189 190 194 197
# [71] 198 199 205 213 220 224 228 231 232 240 242 249 253 255
# [85] 256 259 268 281 286 292 294 319 320 324 327 331 333 334
# [99] 335 336 337 338 339 340 341 342 347 356 357 358 359 360
# [113] 361 382 396 411 421 430 449 451 453 454 455 457 459 462
# [127] 464 465 466 470 471 480 482 483 485 486 489 490 507 513
# [141] 514 531 534 537 540 543 553 554 557 558 559 561 563 565
# [155] 575 576 578 579 580 585 587 615

# for(i in 1:length(jags_data)){
#   print(names(jags_data[i]))
#   if(is.null(dim(jags_data[[i]]))){
#     print(jags_data[[i]])
#   } else {
#     print(dim(jags_data[[i]]))
#   }
# }

for(t in 1:30){
  index <- which(!is.na(jags_data$af[,t]))
  print(all(jags_data$sf[index,t] == jags_data$af[index,t]))
} 


#SIM DATA

k = 10
n = 100

#set.seed(42)
param_list <- list(
  n = n,
  k = k,
  prop.female = 0.5,#0.4652568,
  delta = rep(0.9, k),
  phi.f = rep(0.9, k),
  phi.m = rep(0.9, k),
  gam = rep(0.7, k),
  p.f = rep(0.8, k),
  p.m = rep(0.8, k),
  rho = rep(0.7, k),
  betas = list(beta0 = 1.5, beta1 = 2),
  rand_init = F,
  init = sample(1, n, TRUE)
)

# attach(param_list)

# # # # # Pull individual dataset

jags_data <- sim_dat(param_list)

cjs_data <- format_to_cjs(jags_data)
js_data <- format_to_js(jags_data)
# # Multiple Datasets using parallel
data_list <- sim_cr_dat(parameter_list = param_list, iterations =  100)
# shuffled_list <- replicate_shuffled_data(jags_data, 4)
# 

# # Run JAGS
#jags_data <- sim_cr_dat(parameter_list = param_list, iterations =  100)
# 
# ## MCMC parameters  
par_settings <- list('n.iter' = 1e2,
                     'n.thin' = 1,
                     'n.burn' = 1e2,
                     'n.chains' = 1,
                    'n.adapt' = 1e2)

# Vaillancourt 
# ## Jags parameters and model script
# 
# # Run standard Model
# # 
jags_params <- c("pF", "pM", "phiF", "phiM", "eps")
jags_model <- script_dir %+% "/10_js_mod_standard.R"

z <- list()

for(i in 1:100){
  z[[i]] <- run_jags(jags_data = js_data,
                jags_model  = jags_model,
                jags_params = jags_params,
                par_settings = par_settings,
                debug = F)
}


results <- process_simulation_data(z, param_list)

results %>% group_by(Parameter) %>% summarize(coverage_50 = mean(In_50),
                                              coverage_95 = mean(In_95),
                                              avg_range_50 = mean(Range_50),
                                              avg_range_95 = mean(Range_95),
                                              avg_bias = mean(Bias),
                                              avg_cv = mean(coef_var))


jags_params <- c("pF", "pM", "phiF", "phiM")
jags_model <- script_dir %+% "/10_cjs_mod_standard.R"


y <- run_jags(jags_data = cjs_data,
              jags_model  = jags_model,
              jags_params = jags_params,
              par_settings = par_settings,
              debug = F)


jags_data$n <- NULL
jags_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps")
jags_model <- script_dir %+% "/11_mod_pair_swap_notime.R"

x <- run_jags(jags_data = jags_data,
              jags_model  = jags_model,
              jags_params = jags_params,
              par_settings = par_settings,
              debug = F)




#STARTED AT 8:00pm monday

jags_samples <- run_jags_parallel(cjs_data,
                                  jags_model,
                                  jags_params,
                                  par_settings,
                                  out_dir,
                                  save = F,
                                  outname = "T1_CJS_STD")


# TEST DATA WITH PROGRAM MARK TO SEE RESULTS 
# ADD FLEXIBLE INIT FOR CJS RUNS

# Run Full Model + No Groups
## MCMC parameters  
par_settings <- list('n.iter' = 10,
                     'n.thin' = 1,
                     'n.burn' = 10,
                     'n.chains' = 1,
                     'n.adapt' = 10)


jags_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps")
jags_model <- script_dir %+% "/11_mod_pair_swap_notime.R"

x <- run_jags(jags_data = jags_data,
              jags_model  = jags_model,
              jags_params = jags_params,
              par_settings = par_settings,
              debug = F)


# jags_samples2 <- run_jags_parallel(jags_data,
#                                    jags_model,
#                                    jags_params,
#                                    par_settings,
#                                    out_dir,
#                                    outname = "TESTING_MODEL2")

# 
# # x <- run_jags(jags_data,
# #          jags_model,
# #          jags_params, 
# #          par_settings)
# 
gather_posterior_summary(x$jags_samples) %>%
#add_true_values(param_list) %>% 
plot_caterpillar(params = jags_params[c(1:7)])# +
# geom_point(aes(x = Parameter, y = true), size = 3, alpha = 0.75, color = "darkblue")

# To do

## ADD RECRUITMENT LOGIC FOR SIMULATION STUDY

# Update Model Code
# - Check that the pair-swap sampling makes sense (DONE)
# - Once pair-swap is reasonable, add histories and repartner logic (DONE) 
# - Then make sure Hduck data is generated correctly (DONE)
# - If simulation is reasonable look at 
# - Data Augementation -> need enough spots for all birds (double check when doing model code) +
#  - -> sometimes birds observed with unseen mate (MEET w/ Simon to Discuss)
# - Speedup by optimizing code
# - Speedup by NIMBLIZING Code
# - Write base R code to generate simulation study datasets
# - Write code to sample iteratively for simulation study (MEET w/ Simon First to Discuss)
# - Code up many time step option (maybe not needed for this study)?

# To do in a while
# Choose good parameters for simulation study + objective 
# Deploy simulation study to sharcnet 
# Write code to build statistics highlighting key aspects of the study
# Once code is fast enough and tested run hduck data through sharcnet 
# Maybe do the chat (or  log-likelihood) thing from chapter 1 to suggest correlation exists (do it after running hduck)
# Pull together results
# Write 



# REPARTNER PROCESS and AREPARTNER PROCESS NEED TO BE FIXED (DONE!)
# NIMBLE WONT COMPILE...EIGENVALUE ERRORS 
# JAGS RUNS GREAT
# LARGE DATASETS MAY BE A PROBLEM

# TO DO
# 1. NEED CODE REVIEW FROM SIMON AFTER I DOCUMENT ALL THE STEPS
# 2. DO THE HDUCK CONVERSION 

#k=10,n=100, beta = 3, 12k iter + 2 cores = 70 hours


#SIM DATA

k = 5
n = 50

param_list <- list(
  n = n,
  k = k,
  prop.female = 0.5,
  delta = rep(0.9, k),
  phi.f = rep(0.8, k),
  phi.m = rep(0.8, k),
  gam = rep(0.6, k),
  p.f = rep(0.75, k),
  p.m = rep(0.75, k),
  rho = rep(0.6, k),
  betas = list(beta0 = 1.0, beta1 = 1.5),
  rand_sex = F,
  rand_init = F,
  init  = rep(1,n)# sample(k-1, n, TRUE)
)

par_settings <- list('n.iter' = 1e2, 
                     'n.thin' = 10,
                     'n.burn' = 1e2,
                     'n.chains' = 1,
                     'n.adapt' = 1e2)

jags_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps")# "psi_raw", "psi_cond", "psi_cond2", "male_taken_jt")
jags_model <- script_dir %+% "/11_mod_pair_swap_notime.R"
jags_data <- sim_dat(param_list)


x <- run_jags(jags_data = jags_data,
              jags_model  = jags_model,
              jags_params = jags_params,
              par_settings = par_settings,
              debug = F)


jags_data_list <- replicate_shuffled_data(jags_data, 100)
jags_data_list <- sim_cr_dat(param_list, iterations = 100, ncores = 5)

saveRDS(jags_data_list, out_dir %+% "jags_data_list_rep_study.rds")

x <- run_jags_simulation_parallel(jags_data_list = jags_data_list,
                                  jags_model  = jags_model,
                                  jags_params = jags_params,
                                  par_settings = par_settings,
                                  out_dir = out_dir,
                                  outname = "test",
                                  save = F,
                                  ncores = 6)



post_summary <- extract_sim_posterior(samples)

plot_sim_caterpillar(posterior_summary = post_summary,
                     parameter_name = "gamma",
                     slope = 0,
                     intercept = param_list$gam[1])


plot_sim_caterpillar(posterior_summary = post_summary,
                     parameter_name = "rho",
                     slope = 0,
                     intercept = param_list$rho[1])

plot_sim_caterpillar(posterior_summary = post_summary,
                     parameter_name = "PhiF",
                     slope = 0,
                     intercept = param_list$phi.f[1])

