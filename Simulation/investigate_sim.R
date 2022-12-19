## Load Custom Scripts ---------------------------------------------------------------------------------------------
`%+%`      <- function(a, b) paste0(a, b)
src_dir    <- getwd() #"/home/mdraghic/projects/def-sbonner/mdraghic/mark_recapture_pair_swap/"
out_dir    <- src_dir %+% "/Simulation/Study1/Run5/Output/"
source(file.path(src_dir, "Scripts", "fn_generic.R"))
source(file.path(src_dir, "Scripts", "fn_sim_pair_data.R"))
source(file.path(src_dir, "Scripts", "fn_correlation_estimators.R"))

## Options (ECHO FOR LOGS)
options(echo = TRUE)

# Load packages ---------------------------------------------------------------------------------------------------
libs <- c("tidyverse","RMark")
load_packages(libs, FALSE)
# -----------------------------------------------------------------------------------------------------------------

files <- list.files(out_dir)
scenario_grid <- get_scenarios()

corr_list <- list()
cjs_list <- list()
cchat_list <- list()
lrt_list <- list()
aic_list <- list()
n_list <- list()

for(i in 1:length(files)){
  print(i)
  temp <-  readRDS(out_dir %+% files[[i]])
  corr_list[[i]] <-temp$summary_corr
  cjs_list[[i]] <- temp$summary_cjs
  cchat_list[[i]] <- temp$summ_chat
  lrt_list[[i]] <- temp$summ_lrt
  n_list[[i]] <- temp$summ_n
  aic_list[[i]] <- temp$summ_aic
  gc()
}

summ_corr <- do.call(rbind, corr_list)
summ_cjs <- do.call(rbind, cjs_list)
summ_cchat <- do.call(rbind, cchat_list)
summ_lrt <- do.call(rbind, lrt_list)
summ_n <- do.call(rbind, n_list)
summ_aic <- do.call(rbind, aic_list)

summ_corr <- summ_corr %>% left_join(scenario_grid, by = c("Scenario" = "scenario"))
summ_cjs <- summ_cjs %>% left_join(scenario_grid, by = c("scenario" = "scenario"))
summ_lrt <- summ_lrt %>% left_join(scenario_grid, by = c("scenario" = "scenario"))


summ_corr %>% group_by(Parameter, Parametric, Pearson, Scenario, Truth) %>%
  summarize(MedBias = median(Bias),
            MeanBias = mean(Bias),
            mean_est      = mean(Est),
            median_est = median(Est),
            sd_est    = sd(Est),
            med_se   = median(SE),
            mean_se = mean(SE),
            Cover95 = mean(In95),
            Cover50 = mean(In50),
            Range95 = mean(Range95),
            Range50 = mean(Range50),
            Cover095 = mean(Cover095),
            Cover050 = mean(Cover050),
            rho = first(rho_true),
            gamma = first(gam_true),
            PF     = first(PF),
            PM     = first(PM),
            PhiF   = first(PhiF),
            PhiM   = first(PhiM),
            n_obs  = first(n_obs),
            k = first(k)) %>% 
  filter(Parameter == "rho" & n_obs == 250 & k == 25 & PF == 0.75) %>% 
  ggplot(aes(y = MeanBias, x = rho, col = as.factor(gamma))) +
  geom_line() + 
  geom_point() + 
  facet_grid(Parametric ~ Pearson)# + ylim(c(-0.05,0.05))



summ_lrt %>% group_by(test, scenario) %>%
  summarize(cover5 = mean(1 * (F_pval_partial_pearson   <= 0.05)),
            rho = first(rho_true),
            gamma = first(gam_true),
            PF     = first(PF),
            PM     = first(PM),
            PhiF   = first(PhiF),
            PhiM   = first(PhiM),
            n_obs  = first(n_obs),
            k = first(k)) %>% 
  filter(n_obs == 250 & k == 25 & PF == 0.75) %>% 
  ggplot() +
  geom_boxplot(aes(x = as.factor(rho), y = cover5, col = as.factor(test)))

summ_cjs %>% group_by(Parameter, scenario, Truth, Version) %>% 
  summarize(bias = mean(Bias),
            In95 = mean(In95),
            In95_Partial_Pearson = mean(In95_Partial_Pearson),
            rho = first(rho_true),
            gamma = first(gam_true),
            PF     = first(PF),
            PM     = first(PM),
            PhiF   = first(PhiF),
            PhiM   = first(PhiM),
            n_obs  = first(n_obs),
            k = first(k)) %>% 
  filter(n_obs == 250 & Parameter == "P" & k == 25 & PF == 0.75 & Version=="N" &  gamma == 0)  %>% 
  ggplot() +
  geom_line(aes(x = rho, y = In95)) +
  geom_point(aes(x = rho, y = In95)) +
  geom_line(aes(x = rho, y = In95_Partial_Pearson), col = "blue", lty = 2) +
  geom_point(aes(x = rho, y = In95_Partial_Pearson), col = "blue") +
  facet_wrap(Parameter ~ .) + ylim(c(0.85,1))

summ_cjs %>% group_by(Parameter, scenario, Truth, Version) %>% 
  summarize(bias = mean(Bias),
            In95 = mean(In95),
            In95_Partial_Pearson = mean(In95_Partial_Pearson),
            rho = first(rho_true),
            gamma = first(gam_true),
            PF     = first(PF),
            PM     = first(PM),
            PhiF   = first(PhiF),
            PhiM   = first(PhiM),
            n_obs  = first(n_obs),
            k = first(k)) %>% 
  filter(Parameter == "Phi" & k == 25 & Version=="R")  %>% 
  ggplot() +
  geom_line(aes(x = gamma, y = In95, col = as.factor(rho))) +
  geom_point(aes(x = gamma, y = In95, col = as.factor(rho))) +
  geom_line(aes(x = gamma, y = In95_Partial_Pearson, col = as.factor(rho)), lty = 2) +
  geom_point(aes(x = gamma, y = In95_Partial_Pearson, col = as.factor(rho))) +
  facet_wrap(PF ~ n_obs) + ylim(c(0.75,1)) + geom_hline(yintercept = 0.95, col = "red")

summ_cjs %>% group_by(Parameter, scenario, Truth, Version) %>% 
  summarize(bias = mean(Bias),
            In95 = mean(In95),
            In95_Partial_Pearson = mean(In95_Pearson),
            rho = first(rho_true),
            gamma = first(gam_true),
            PF     = first(PF),
            PM     = first(PM),
            PhiF   = first(PhiF),
            PhiM   = first(PhiM),
            n_obs  = first(n_obs),
            k = first(k)) %>% 
  filter(Parameter == "P" & k == 25 & Version=="N")  %>% 
  ggplot() +
  geom_line(aes(x = rho, y = In95, col = as.factor(gamma))) +
  geom_point(aes(x =rho, y = In95, col = as.factor(gamma))) +
  geom_line(aes(x = rho, y = In95_Partial_Pearson, col = as.factor(gamma)), lty = 2) +
  geom_point(aes(x = rho, y = In95_Partial_Pearson, col = as.factor(gamma))) +
  facet_wrap(PF ~ n_obs) + ylim(c(0.75,1)) + geom_hline(yintercept = 0.95, col = "red")
