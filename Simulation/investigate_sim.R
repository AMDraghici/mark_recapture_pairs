## Load Custom Scripts ---------------------------------------------------------------------------------------------
`%+%`      <- function(a, b) paste0(a, b)
src_dir    <- getwd() #"/home/mdraghic/projects/def-sbonner/mdraghic/mark_recapture_pair_swap/"
out_dir    <- src_dir %+% "/Simulation/Run6/Output/"
source(file.path(src_dir, "Scripts", "fn_generic.R"))
source(file.path(src_dir, "Scripts", "fn_sim_pair_data.R"))
source(file.path(src_dir, "Scripts", "fn_correlation_estimators.R"))

## Options (ECHO FOR LOGS)
options(echo = TRUE)

# Load packages ---------------------------------------------------------------------------------------------------
libs <- c("tidyverse","RMark")
load_packages(libs, FALSE)
# -----------------------------------------------------------------------------------------------------------------

# Grab Results
files <- list.files(out_dir)
nruns <- length(files)
scenario_grid <- get_scenarios()

# Unpack Summaries
out_list     <- lapply(1:nruns, function(i) readRDS(out_dir %+% files[i]))
summ_corr    <- do.call(rbind, lapply(1:nruns, function(i) out_list[[i]]$summary_corr))
summ_cjs     <- do.call(rbind, lapply(1:nruns, function(i) out_list[[i]]$summary_cjs))
summ_chat    <- do.call(rbind, lapply(1:nruns, function(i) out_list[[i]]$summ_chat))
summ_lrt     <- do.call(rbind, lapply(1:nruns, function(i) out_list[[i]]$summ_lrt))
summ_aic     <- do.call(rbind, lapply(1:nruns, function(i) out_list[[i]]$summ_aic))
summ_n       <- do.call(rbind, lapply(1:nruns, function(i) out_list[[i]]$summ_n))
summ_gam_delta <- do.call(rbind, lapply(1:nruns, function(i) out_list[[i]]$summ_gam_delta))

# Add Scenario Settings
summ_corr <-  summ_corr %>% left_join(scenario_grid, by = c("scenario"))
summ_cjs  <-  summ_cjs %>% left_join(scenario_grid, by = c("scenario"))
summ_lrt  <-  summ_lrt %>% left_join(scenario_grid, by = c("scenario"))
summ_chat  <- summ_chat %>% left_join(scenario_grid, by = c("scenario"))
summ_aic  <-  summ_aic %>% left_join(scenario_grid, by = c("scenario"))
summ_n  <-    summ_n %>% left_join(scenario_grid, by = c("scenario"))
summ_gam_delta <- summ_gam_delta %>% left_join(scenario_grid, by = c("scenario"))


saveRDS(summ_corr, src_dir %+% "/Output/summ_corr.rds")
saveRDS(summ_cjs,  src_dir %+% "/Output/summ_cjs.rds")
saveRDS(summ_lrt,  src_dir %+% "/Output/summ_lrt.rds")
saveRDS(summ_chat, src_dir %+% "/Output/summ_chat.rds")
saveRDS(summ_aic,  src_dir %+% "/Output/summ_aic.rds")
saveRDS(summ_n,    src_dir %+% "/Output/summ_n.rds")
saveRDS(summ_gam_delta, src_dir %+% "/Output/summ_gam_delta.rds")

summ_corr <- readRDS(src_dir %+% "/Output/summ_corr.rds")
summ_cjs <- readRDS(src_dir %+% "/Output/summ_cjs.rds")
summ_lrt <- readRDS(src_dir %+% "/Output/summ_lrt.rds")
summ_chat <- readRDS(src_dir %+% "/Output/summ_chat.rds")
summ_aic <- readRDS(src_dir %+% "/Output/summ_aic.rds")
summ_n <- readRDS(src_dir %+% "/Output/summ_n.rds")

mc_corr <- summ_corr %>% group_by(Parameter, 
                                  Parametric, 
                                  Pearson,
                                  scenario, 
                                  Truth) %>%
  mutate(Pearson = ifelse(Pearson == 0, "MLE", ifelse(Pearson == 2, "Partial-Pearson", "Pearson")),
         Parametric = ifelse(Parametric == 0, "Non-Parametric", "Semi-Parametric")) %>% 
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
            k      = first(k),
            Beta0 = first(Beta0),
            Accept95_1 = mean(1 * (Pval0_1 > 0.05)),
            Accept95_2 = mean(1 * (Pval0_2 > 0.05)),
            imputed_pairs = first(imputed_pairs))



# Investigating Performance of hat Rho
mc_corr %>% 
  filter(Parameter == "rho" & Beta0 > 1 & PF > 0.4 & PF != 0.48 & PF != 0.7 & imputed_pairs==T) %>% 
  ggplot(aes(y = Cover95, x = rho, col = as.factor(gamma))) +
  geom_line() + 
  geom_point() + 
  facet_grid(Pearson + as.factor(n_obs) + as.factor(k) ~Parametric + as.factor(PF)) + 
  ylim(c(0.6,1)) +
  geom_hline(yintercept = 0.95, lty = 2) + 
  labs(x = expression(rho), y = "Mean Bias", col = expression(gamma)) +
  theme(legend.position = "bottom")

mc_corr %>% 
  filter(Parameter == "rho" & Beta0 > 1 & n_obs == 250 & k == 25 & PF == 0.75 & PM == 0.75 & imputed_pairs==T) %>% 
  ggplot(aes(y = Cover95, x = rho, col = as.factor(gamma))) +
  geom_line() + 
  geom_point() + 
  facet_grid(Pearson ~Parametric) + 
  ylim(c(0.75,1)) +
  geom_hline(yintercept = 0.95, lty = 2) + 
  labs(x = expression(rho), y = "Mean 95% Confidence Interval Coverage", col = expression(gamma)) +
  theme(legend.position = "bottom")



mc_corr %>% 
  filter(Parameter == "rho" & Beta0 > 1 & n_obs == 250 & k == 25 & PF == 0.75 & PM == 0.75 & imputed_pairs==T) %>% 
  ggplot(aes(y = Cover50, x = rho, col = as.factor(gamma))) +
  geom_line() + 
  geom_point() + 
  facet_grid(Pearson ~Parametric) + 
  ylim(c(0,1)) +
  geom_hline(yintercept = 0.5, lty = 2) + 
  labs(x = expression(rho), y = "Mean 50% Confidence Interval Coverage", col = expression(gamma)) +
  theme(legend.position = "bottom")


mc_corr %>% 
  filter(Parameter == "rho" & Beta0 > 1 & k == 25 & n_obs == 250 & Pearson == "Partial-Pearson" & Parametric == "Semi-Parametric" & PhiF == 0.8 & PF == 0.75 & PM == 0.75 & imputed_pairs==T) %>% 
  ggplot(aes(y = 1-Cover095, x = rho, col = as.factor(gamma))) +
  geom_line() + 
  geom_point() + 
  # facet_grid(n_obs ~k) + 
  ylim(c(0,1)) +
  geom_hline(yintercept = c(0.05,0.8), lty = 2) + 
  labs(x = expression(rho), y = "Average Coverage of 0 by 95% Confidence Interval", col = expression(gamma)) +
  theme(legend.position = "bottom")


# Investigating Performance of hat gamma
mc_corr %>% 
  filter(Parameter == "gamma" & Beta0 > 1 & n_obs == 250 & k == 25 &  PhiF == 0.8 &  PF == 0.75 & PM == 0.75 & imputed_pairs==T) %>% 
  ggplot(aes(y = MedBias, x = gamma, col = as.factor(rho))) +
  geom_line() + 
  geom_point() + 
  facet_grid(Pearson ~Parametric) + 
  ylim(c(-0.03,0.03)) +
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = expression(gamma), y = "Mean Bias", col = expression(rho)) +
  theme(legend.position = "bottom")

mc_corr %>% 
  filter(Parameter == "gamma" & Beta0 > 1 & n_obs == 250 & k == 25 & PF == 0.75 & PM == 0.75 & imputed_pairs==T) %>% 
  ggplot(aes(y = Cover95, x = gamma, col = as.factor(rho))) +
  geom_line() + 
  geom_point() + 
  facet_grid(Pearson ~Parametric) + 
  ylim(c(0.75,1)) +
  geom_hline(yintercept = 0.95, lty = 2) + 
  labs(x = expression(gamma), y = "Mean 95% Confidence Interval Coverage", col = expression(rho)) +
  theme(legend.position = "bottom")


summ_gam_delta %>% 
  group_by(Parameter, 
           Pearson,
           scenario, 
           Truth)  %>% 
  summarize(mean_est      = mean(Est),
            median_est = median(Est),
            sd_est    = sd(Est),
            Cover95 = mean(In95),
            Cover50 = mean(In50),
            rho = first(rho_true),
            gamma = first(gam_true),
            PF     = first(PF),
            PM     = first(PM),
            PhiF   = first(PhiF.x),
            PhiM   = first(PhiM.x),
            n_obs  = first(n_obs),
            k      = first(k),
            Beta0 = first(Beta0),
            imputed_pairs = first(imputed_pairs)) %>% 
  filter(Parameter == "gamma" & Beta0 > 1 & n_obs == 250 & k == 25 & PF == 0.75 & PM == 0.75 & imputed_pairs==T) %>% 
  ggplot(aes(y = Cover95, x = gamma, col = as.factor(rho))) +
  geom_line() + 
  geom_point() + 
  facet_grid(Pearson ~.) + 
  ylim(c(0.75,1)) +
  geom_hline(yintercept = 0.95, lty = 2) + 
  labs(x = expression(gamma), y = "Mean 95% Confidence Interval Coverage", col = expression(rho)) +
  theme(legend.position = "bottom")


mc_corr %>% 
  filter(Parameter == "gamma"& Beta0 > 1 & n_obs == 250 & k == 25 & PF == 0.75 & PM == 0.75 & imputed_pairs==T) %>% 
  ggplot(aes(y = Cover50, x = gamma, col = as.factor(rho))) +
  geom_line() + 
  geom_point() + 
  facet_grid(Pearson ~Parametric) + 
  ylim(c(0,1)) +
  geom_hline(yintercept = 0.5, lty = 2) + 
  labs(x = expression(gamma), y = "Mean 50% Confidence Interval Coverage", col = expression(rho)) +
  theme(legend.position = "bottom")


mc_corr %>% 
  filter(Parameter == "gamma" & Beta0 > 1 & n_obs == 250 & k == 25 & PhiF == 0.8 & PF == 0.75 & PM == 0.75 & imputed_pairs==T) %>% 
  ggplot(aes(y = 1 - Accept95_1, x = gamma, col = as.factor(rho))) +
  geom_line() + 
  geom_point() + 
  facet_grid(Pearson ~ Parametric) + 
  ylim(c(0,1)) +
  geom_hline(yintercept = c(0.8), lty = 2) + 
  labs(x = expression(rho), y = "Mean Coverage of 0 by 95% Confidence Interval", col = expression(rho)) +
  theme(legend.position = "bottom")


mc_corr %>% 
  filter(Parameter == "rho" & Beta0 < 1 & n_obs == 250 & k == 25 & PhiF == 0.8 & PF == 0.75 & PM == 0.75 & imputed_pairs==T) %>% 
  ggplot(aes(y = 1 - Accept95_1, x = rho, col = as.factor(gamma))) +
  geom_line() + 
  geom_point() + 
  facet_grid(Pearson ~ Parametric) + 
  ylim(c(0,1)) +
  geom_hline(yintercept = c(0.8), lty = 2) + 
  labs(x = expression(rho), y = "Mean Coverage of 0 by 95% Confidence Interval", col = expression(rho)) +
  theme(legend.position = "bottom")
# 
# summ_cjs2 <- summ_cjs %>% mutate(CChatAdj_Partial_Pearson2 = 2^(log(CChatAdj_Partial_Pearson,base = 2)/2),
#                                  UBAdj_Partial_Pearson2    = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Partial_Pearson2), alpha = 0.05)[["ub"]],
#                                  LBAdj_Partial_Pearson2    = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Partial_Pearson2), alpha = 0.05)[["lb"]],
#                                  In95_Partial_Pearson2     = 1*(Truth <= UBAdj_Partial_Pearson2 & Truth >= LBAdj_Partial_Pearson2))

# Comparing Coverage of N model
mc_cjs <- summ_cjs %>% group_by(Parameter, scenario, Truth, Version) %>% 
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
            Beta0         = first(Beta0),
            k = first(k),
            imputed_pairs = first(imputed_pairs))

mc_cjs2 <- summ_cjs %>% group_by(scenario, Version) %>% 
  summarize(chat          = mean(Deviance_chat),
            chat_pearson  = mean(Pearson_chat),
            chat_fletcher = mean(Fletcher_chat),
            CChat         = mean(CChatAdj_Partial_Pearson),
            rho           = first(rho_true),
            gamma         = first(gam_true),
            PF            = first(PF),
            PM            = first(PM),
            PhiF          = first(PhiF),
            PhiM          = first(PhiM),
            n_obs         = first(n_obs),
            k             = first(k),
            imputed_pairs = first(imputed_pairs),
            Beta0         = first(Beta0)) %>% 
  pivot_longer(cols = c("chat", "chat_pearson", "chat_fletcher"), 
               names_to = "Estimator", values_to = "chat") %>% 
  mutate(Estimator = ifelse(Estimator == "chat", "Deviance",
                            ifelse(Estimator == "chat_pearson", "Pearson","Fletcher")))

mc_cjs %>% 
  filter(Parameter == "P" & Beta0 > 1 & n_obs == 250 & k == 25 & PF == 0.75 & PhiF == 0.8 & imputed_pairs==T)  %>% 
  ggplot() +
  geom_line(aes(x = rho, y = In95, col = as.factor(gamma))) +
  geom_point(aes(x = rho, y = In95, col = as.factor(gamma))) +
  geom_line(aes(x = rho, y = In95_Partial_Pearson, col = as.factor(gamma)), lty = 2) +
  geom_point(aes(x = rho, y = In95_Partial_Pearson, col = as.factor(gamma))) +
  ylim(c(0.8,1)) +
  geom_hline(yintercept = 0.95) +
  facet_grid(. ~ Version) + 
  labs(x = expression(rho), y = "Mean 95% Confidence Interval Coverage", col = expression(rho)) +
  theme(legend.position = "bottom")

mc_cjs %>% 
  filter(Parameter == "Phi" & Beta0 > 1 & k == 25 & n_obs == 250 & PF == 0.75 & PhiF == 0.8 & imputed_pairs==T)  %>% 
  ggplot() +
  geom_line(aes(x = gamma, y = In95, col = as.factor(rho))) +
  geom_point(aes(x = gamma, y = In95, col = as.factor(rho))) +
  geom_line(aes(x = gamma, y = In95_Partial_Pearson, col = as.factor(rho)), lty = 2) +
  geom_point(aes(x = gamma, y = In95_Partial_Pearson, col = as.factor(rho))) +
  ylim(c(0.8,1)) +
  geom_hline(yintercept = 0.95) +
  facet_grid(. ~ Version) + 
  labs(x = expression(gamma), y = "Mean 95% Confidence Interval Coverage", col = expression(rho)) +
  theme(legend.position = "bottom")

mc_cjs2 %>% 
  filter(n_obs == 250 & k == 25 & PF == 0.75 & PhiF == 0.8 & imputed_pairs==T)  %>% 
  ggplot() +
  geom_line(aes(x = rho, y = chat, col = as.factor(gamma))) +
  geom_point(aes(x = rho, y = chat, col = as.factor(gamma))) +
  geom_line(aes(x = rho, y = CChat, col = as.factor(gamma)), lty = 2) +
  # geom_point(aes(x = rho, y = CChat, col = as.factor(gamma)), lty = 2) +
  geom_hline(yintercept = 0.95) +
  facet_grid(Estimator ~ Version) + 
  labs(x = expression(gamma), y = "Mean 95% Confidence Interval Coverage", col = expression(rho)) +
  theme(legend.position = "bottom")



success <- summ_cjs$scenario %>% unique()
scenario_grid %>% filter(!scenario %in% success) %>% pull(scenario)
scenario_grid %>% filter(!scenario %in% success) %>% pull(scenario) -> x

paste0(x, ",")

