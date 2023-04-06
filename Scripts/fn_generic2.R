## Generic Functions -------------------------------------------------------------------------------------

## Concatenate Strings inline 
`%+%` <- function(a,b) paste0(a,b)

## Read in vector of libraries 
## If package is not installed then do so
load_packages <- function(libs, install=TRUE){
  
  if(install==TRUE){
    #Read through vector libraries called by user
    for (pkg in libs) {
      #If package is missing from installed.packages then install it
      if (!pkg %in% rownames(installed.packages())) {
        install.packages(pkg)
      }
    }
  }
  
  ## Attach our libraries
  lapply(libs, require, character.only = TRUE)
}

#Inline Paste0 function
`%+%` <- function(a, b) paste0(a, b)

#Odds function
odds <- function(p){
  out <- p/(1-p)
  return(out)
}

# odds ratio
or <- function(p,q){
  odds_p <- odds(p)
  odds_q <- odds(q)
  out <- odds_p/odds_q
  return(out)
}

# odds product
op <- function(p,q){
  odds_p <- odds(p)
  odds_q <- odds(q)
  out <- odds_p*odds_q
  return(out)
}

#Logistic Function for probability
logit <- function(p){
  out <- log(odds(p))
  return(out)
}

# #Inverse Logistic Function
inv.logit <- function(x){
  out <- (1+exp(-x))^(-1)
  return(out)
}

## S/R Distribution Functions -------------------------------------------------------------------------------------

#Extract joint binomial parameters
compute_jbin_param_cjs <- function(prob.f,prob.m){
  #Compute Components
  prob.prod  <- prob.m * prob.f
  sig.prob.f <- sqrt(prob.f*(1-prob.f))
  sig.prob.m <- sqrt(prob.m*(1-prob.m))
  sig.prod   <- sig.prob.f*sig.prob.m
  ##Upper and Lower bounds for cov of joint binomial
  lub        <-  sqrt(pmin(or(prob.f,prob.m),1/or(prob.f,prob.m)))
  glb        <- -sqrt(pmin(op(prob.f,prob.m),1/op(prob.f,prob.m)))
  #Return values
  out <- list(prob.prod       = prob.prod,
              sig.prob.f      = sig.prob.f,
              sig.prob.m      = sig.prob.m,
              sig.prod        = sig.prod,
              cor_lower_bound = glb,
              cor_upper_bound = lub)
  return(out)
}

#Generate Derived Probabilities
compute_jbin_cjs <- function(prob.f,prob.m,corr){
  
  parameters      <- compute_jbin_param_cjs(prob.f,prob.m)
  cor_upper_bound <- parameters[["cor_upper_bound"]]
  cor_lower_bound <- parameters[["cor_lower_bound"]]
  sig.prob.f      <- parameters[["sig.prob.f"]]
  sig.prob.m      <- parameters[["sig.prob.m"]]
  
  # If we are on a boundary condition update correlation structure 
  boundary_f                <-  which(prob.f == 1| prob.f == 0)
  boundary_m                <-  which(prob.m == 1| prob.m == 0)
  boundary                  <-  sort(unique(c(boundary_f, boundary_m)))
  corr[boundary]            <- 0 #correlation is zero if 1 or 0 for both probs
  cor_upper_bound[boundary] <- 0
  cor_lower_bound[boundary] <- 0
  
  if(any(round(corr,3) > round(cor_upper_bound,3)) || any(round(corr,3) < round(cor_lower_bound,3))){
    stop("Correlation " %+% corr %+% " is outside of the bounds [" %+% 
           as.character(cor_lower_bound) %+% " , " %+%
           as.character(cor_upper_bound) %+% "]")
  }
  
  ###Joint Probability Distn 
  prob.mf <- corr * sig.prob.m * sig.prob.f + (prob.f*prob.m)
  prob.f0 <- prob.f - prob.mf
  prob.m0 <- prob.m - prob.mf
  prob.00 <- 1 - prob.f0 - prob.m0 - prob.mf
  
  #List of parameters
  out <- list(prob.mf = prob.mf,
              prob.f0 = prob.f0,
              prob.m0 = prob.m0,
              prob.00 = prob.00)
  
  #Return values
  return(out)
}



# Run CJS Model using Mark
run_cjs_model_marked <- function(cjs_data,
                                 version = "B"){
  
  # Define survival/recapture parameter space
  phi.grp     <- list(formula = ~sex)
  p.grp       <- list(formula = ~sex)
  
  #Process Mark Data (Extract Recapture/sex)
  # Add sex grouping for B,R,S only
  if(version == "N"){
    dat.process <- cjs_data$x %>% 
      as.data.frame() %>% 
      unite("ch",sep="") 
    
    mark.process <- marked::process.data(data   = dat.process,
                                         model  = "CJS")
  } else {
    dat.process <- cjs_data$x %>% 
      as.data.frame() %>% 
      unite("ch",sep="") %>% 
      mutate(sex = as.factor(ifelse(cjs_data$female == 1, "F", "M"))) %>%
      arrange(sex)
    
    mark.process <-  marked::process.data(data   = dat.process,
                                          model  = "CJS",
                                          groups = "sex")
  }
  
  #Design Data List
  mark.ddl <-  marked::make.design.data(mark.process) 
  
  # Choose Appropriate CJS Model Settings
  if(version == "B"){
    model.parameters <- list(Phi = phi.grp, p   = p.grp)
    phi.name         <- c("PhiF","PhiM")
    p.name           <- c("PF","PM")
    param.names      <- c(phi.name,p.name)
  } else if(version == "S"){
    model.parameters <- list(Phi = phi.grp)
    phi.name         <- c("PhiF","PhiM")
    p.name           <- c("P")
    param.names      <- c(phi.name,p.name)
  } else if(version == "R"){
    model.parameters <- list(p   = p.grp)
    phi.name         <- c("Phi")
    p.name           <- c("PF","PM")
  } else if(version == "N"){
    model.parameters <- list()
    phi.name         <- c("Phi")
    p.name           <- c("P")
  }
  
  # Add names
  param.names <- c(phi.name, p.name)
  
  #Generate Estimates
  mark_out <-  marked::crm(data             = mark.process,
                           ddl              = mark.ddl,
                           model.parameters = model.parameters,
                           hessian          = T)
  
  mark_out <- mark_out$results
  
  #Extract values of use
  gof.stats <- data.frame("lnl"       = mark_out$neg2lnl,
                          "npar"        = (mark_out$AIC - mark_out$neg2lnl)/2)  %>% 
    rename(`-2lnl` = "lnl")
  
  mark.stats <- rbind(mark_out$reals$Phi %>%  select(estimate, se, lcl, ucl),
                      mark_out$reals$p  %>% select(estimate, se, lcl, ucl)) %>% 
    filter(se > 0) %>% 
    rename("Est" = estimate, "SE" = se, "LB" = lcl, "UB" = ucl) %>%
    mutate(Parameter = param.names,
           Version   = version)
  
  # Add Pearson and Flecthers Chat
  if(version == "B"){
    PhiF      <- mark.stats %>% filter(Parameter == "PhiF") %>% pull(Est)
    PhiM      <- mark.stats %>% filter(Parameter == "PhiM") %>% pull(Est)
    PM        <- mark.stats %>% filter(Parameter == "PM") %>% pull(Est)
    PF        <- mark.stats %>% filter(Parameter == "PF") %>% pull(Est)
    chat_list <- compute_chat_mark(cjs_data = cjs_data,
                                   PhiF     = PhiF,
                                   PhiM     = PhiM,
                                   PF       = PF,
                                   PM       = PM,
                                   ungroup = F)
    
  } else if(version == "S"){
    PhiF      <- mark.stats %>% filter(Parameter == "PhiF") %>% pull(Est)
    PhiM      <- mark.stats %>% filter(Parameter == "PhiM") %>% pull(Est)
    PM        <- PF   <- mark.stats %>% filter(Parameter == "P") %>% pull(Est)
    chat_list <- compute_chat_mark(cjs_data = cjs_data,
                                   PhiF     = PhiF,
                                   PhiM     = PhiM,
                                   PF       = PF,
                                   PM       = PM,
                                   ungroup  = F)
  } else if(version == "R"){
    PhiF      <- PhiM <- mark.stats %>% filter(Parameter == "Phi") %>% pull(Est)
    PM        <- mark.stats %>% filter(Parameter == "PM") %>% pull(Est)
    PF        <- mark.stats %>% filter(Parameter == "PF") %>% pull(Est)
    chat_list <- compute_chat_mark(cjs_data = cjs_data,
                                   PhiF     = PhiF,
                                   PhiM     = PhiM,
                                   PF       = PF,
                                   PM       = PM,
                                   ungroup  = F)
  } else if(version == "N"){
    PhiF <- PhiM <- mark.stats %>% filter(Parameter == "Phi") %>% pull(Est)
    PM   <- PF   <- mark.stats %>% filter(Parameter == "P") %>% pull(Est)
    chat_list    <- compute_chat_mark(cjs_data = cjs_data,
                                      PhiF     = PhiF,
                                      PhiM     = PhiM,
                                      PF       = PF,
                                      PM       = PM,
                                      ungroup = T)
  }
  
  mark_results <- cbind(mark.stats,gof.stats) %>% 
    mutate(Pearson_chat  = chat_list[["Pearson_chat"]],
           Fletcher_chat = chat_list[["Fletcher_chat"]])
  
  #Return Output
  return(mark_results)
  
}


# Proposed Chat Estimator 
compute_proposed_chat <- function(corr1,
                                  corr2, 
                                  linear = F){
  
  if(linear){
    cchat1 <- ifelse(corr1 >= 0, 1+corr1, 1 + 0.5  * corr1)
    cchat2 <- ifelse(corr2 >= 0, 1+corr2, 1 + 0.5  * corr2)
    cchat <- cchat1 * cchat2
  } else {
    cchat <- 2^(corr1 + corr2)
  }
  
  
  return(cchat)
}


#Generate MARK CI for Probabilities
# Arg: Prob is value btwn [0, 1]
# Arg: Se is a standard error value
# Arg: Alpha is a confidence level
compute_mark_ci <- function(prob,se,alpha=0.05){
  var.logit <- (se^2)/((prob-1)^2*prob^2)                    # Var Delta Method
  se.logit <- sqrt(var.logit)                                # SD transform
  est.logit <- logit(prob)                                   # Logit
  ub.logit <- est.logit + round(qnorm(1-alpha/2),2)*se.logit # CLT UB
  lb.logit <- est.logit - round(qnorm(1-alpha/2),2)*se.logit # CLT LB
  lb.real <- exp(lb.logit)/(1+exp(lb.logit))                 # Backtransform LB
  ub.real <- exp(ub.logit)/(1+exp(ub.logit))                 # Backtransform UB
  return(list(lb = lb.real,ub = ub.real))                    # Return interval
} 


# Compute AIC using MARKED program logic
compute_aic_mark <- function(ll, k, n, chat, cc = 0){
  aic  <- ll/chat + 2*k + cc * (2*k * (k+1))/(n-k-1)
  return(aic)
}

# Compute AIC Summary 
compute_aic_summ <- function(summ_cjs,
                             n_eff,
                             iter,
                             scenario){
  
  summ_aic <- summ_cjs %>%
    group_by(Version) %>%
    summarize(ll                     = first(`-2lnl`),
              pars                   = first(npar),
              CChat_Likelihood       = prod(unique(CChatAdj_Likelihood))#,
              # CChat_Pearson          = prod(unique(CChatAdj_Pearson)),
              # CChat_Partial_Pearson  = prod(unique(CChatAdj_Partial_Pearson))
    ) %>% 
    ungroup() %>% 
    mutate(AIC                       = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = rep(1,4),              cc = 0),
           AICC                      = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = rep(1,4),              cc = 1),
           AIC_chat_likelihood       = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = CChat_Likelihood,      cc = 0),
           AICC_chat_likelihood      = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = CChat_Likelihood,      cc = 1),
           # AIC_chat_pearson          = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = CChat_Pearson,         cc = 0),
           # AICC_chat_pearson         = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = CChat_Pearson,         cc = 1),
           # AIC_chat_partial_pearson  = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = CChat_Partial_Pearson, cc = 0),
           # AICC_chat_partial_pearson = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = CChat_Partial_Pearson, cc = 1),
           iter                      = iter,
           scenario                  = scenario)
  
  return(summ_aic)
}

# Compute chi_t value for last observation prob in encounter history of CJS
compute_chi_t <- function(phi,p,k,K){
  ifelse(k==K, 1, (1-phi) + phi * (1-p) * compute_chi_t(phi,p, k+1, K))
}

# Helper Function Converts string of numbers to vector numeric format
convert_history_numeric <- function(hist){
  as.integer(unlist(strsplit(hist, "")))
}

# Compute encounter history probability 
encounter_history_prob <- function(hist, phi, p){
  alive <- which(hist == 1)
  et    <- dplyr::first(alive)
  ek    <- dplyr::last(alive)
  K     <- length(hist)
  
  occasions <- (ek-et)
  
  seen <- 1 *(ek > et) * sum(hist[(et+1):ek])
  prob <- (phi ^ occasions) * (p ^ seen) * ((1-p) ^ (occasions - seen)) *  compute_chi_t(phi, p, ek, K)
  
  return(prob)
}

# Convert CJS data from application/sim into History table
prep_encounter_history_df <- function(cjs_data){
  hist_data <- cjs_data$x %>% 
    as.data.frame() %>% 
    unite("ch",sep="") %>% 
    mutate(sex = as.factor(ifelse(cjs_data$female == 1, "F", "M"))) %>%
    arrange(sex) 
  
  temp_tbl <- hist_data[,1:2] %>% table
  histories <- rownames(temp_tbl)
  Oi_F <- unname(temp_tbl[,1])
  Oi_M <- unname(temp_tbl[,2])
  hist_data <- data.frame(History = histories,
                          Oi_F = Oi_F,
                          Oi_M = Oi_M,
                          Oi   = Oi_F + Oi_M) 
  
  
  cohort <- sapply(1:nrow(hist_data), function(x) first(which(convert_history_numeric(hist_data$History[x]) == 1)))
  hist_data$cohort <- cohort
  return(hist_data)
}

# Get C-hat estimator by using Pearson's or Flecther's methods
compute_chat_mark <- function(cjs_data, 
                              PhiF, 
                              PhiM, 
                              PF, 
                              PM,
                              ungroup = F){
  
  # Get encounter history df 
  hist_data <- prep_encounter_history_df(cjs_data)
  
  # Add history probs
  hist_data$prob_f <- sapply(1:nrow(hist_data), 
                             function(x) 
                               encounter_history_prob(convert_history_numeric(hist_data[x,1]), PhiF, PF))
  hist_data$prob_m <- sapply(1:nrow(hist_data), 
                             function(x) 
                               encounter_history_prob(convert_history_numeric(hist_data[x,1]), PhiM, PM))
  
  
  # Encounter Histories
  K <- cjs_data$k
  
  # Count of Observations in each Cohort
  N_cohort_df <- hist_data %>%
    group_by(cohort) %>% 
    summarize(Nf               = sum(Oi_F,na.rm = T),
              Nm               = sum(Oi_M,na.rm = T),
              N                = sum(Oi,na.rm=T),
              f_correction     = Nf - Nf * sum(prob_f),
              m_correction     = Nm - Nm * sum(prob_m),
              ungrp_correction = N - N *sum(prob_f))
  
  
  # Get Expected Cell Counts and Terms for Sums 
  hist_data <- hist_data %>%
    inner_join(N_cohort_df, by = "cohort") %>% 
    mutate(Ei_F          = Nf * prob_f,
           Ei_M          = Nm * prob_m,
           Ei            = N * prob_f,
           pearson_i_f   = ifelse(Nf == 0, 0, (Oi_F - Ei_F)^2/(Ei_F)),
           pearson_i_m   = ifelse(Nm == 0, 0, (Oi_M - Ei_M)^2/(Ei_M)),
           pearson_i     = (Oi - Ei)^2/Ei,
           Obs_Ratio_i_f = Oi_F/Ei_F,
           Obs_Ratio_i_m = Oi_M/Ei_M,
           Obs_Ratio_i   = Oi/Ei) 
  
  # Gather Pearsons Chi-Sq and sum(Oi/Ei)
  pearson_statistics <- hist_data %>% 
    group_by(cohort) %>% 
    summarize(pearson         = sum(pearson_i_m) + sum(pearson_i_f) + first(f_correction) + first(m_correction),
              pearson_ungrp   = sum(pearson_i) + first(ungrp_correction),
              Obs_ratio       = sum(Obs_Ratio_i_f,na.rm = T) + sum(Obs_Ratio_i_m,na.rm =T),
              Obs_ratio_ungrp = sum(Obs_Ratio_i,na.rm = T)) 
  
  # Possible Outcomes
  possible_histories_df <- hist_data %>% 
    mutate(possible_n = 2^(K-cohort)) %>% 
    group_by(cohort) %>% 
    summarize(incl_f = 1 * (max(Oi_F) > 0),
              incl_m = 1 * (max(Oi_M) > 0),
              N_raw_i = first(possible_n)) %>% 
    mutate(N_cohort       = N_raw_i * (incl_f + incl_m),
           N_cohort_ungrp = N_raw_i) %>% 
    ungroup() 
  
  # If Ungroup true then assume no diff between sexes (then PF=PM and PhiF = PhiM)
  if(ungroup){
    cat("Ungroup selected: Assuming PhiF = Phi and PF = P...", "\n")
    pearson_chisq      <- pearson_statistics %>% pull(pearson_ungrp) %>% sum
    obs_ratio          <- pearson_statistics %>% pull(Obs_ratio_ungrp) %>% sum
    possible_histories <- possible_histories_df %>% pull(N_cohort_ungrp) %>% sum
    pearson_df         <- possible_histories - K - 1
  } else {
    pearson_chisq      <- pearson_statistics %>% pull(pearson) %>% sum
    obs_ratio          <- pearson_statistics %>% pull(Obs_ratio) %>% sum
    possible_histories <- possible_histories_df %>% pull(N_cohort) %>% sum
    pearson_df         <- possible_histories - 2*K - 1
  }
  
  # Pearson's Method
  Pearson_chat <- pearson_chisq/pearson_df
  
  # Fletcher's Method
  N              <- possible_histories                  # All possible outcomes
  n              <- length(unique(hist_data$cohort))    # Number of unique cohorts
  s_bar          <- (obs_ratio - N)/(N-n)           # Fletcher's Term
  Fletcher_chat  <- Pearson_chat/(1+s_bar) # Fletcher's Chat
  
  # Return choosen c-hat
  return(list(Pearson_chat = Pearson_chat,
              Fletcher_chat = Fletcher_chat))
}


compute_bca_gamma <- function(btrsp, 
                              est,
                              ybar,
                              n_eff_gamma,
                              pf_mark,
                              pm_mark,
                              rho, 
                              phif_mark,
                              phim_mark,
                              alpha){
  
  # Number of successes 
  n_success <- ybar * n_eff_gamma
  
  # Bias correction 
  z0 <- qnorm(mean(1*(btrsp < est)))
  
  # Prob FM
  PFM    = compute_jbin_cjs(prob.f = pf_mark,
                            prob.m = pm_mark,
                            corr   = rho)$prob.mf
  
  # JackKnife  the Ybar estimate for one iteration 
  ybar_i <- (c(rep(n_success/(PFM * (n_eff_gamma-1)),n_eff_gamma-n_success), # remove only denominator failure number of times
               rep((n_success-1)/(PFM * (n_eff_gamma-1)),n_success))) # remove only numerator sucess number of times
  
  # Convert from estimate of PhiFM * PFM to estimate of gamma 
  gamma_i <- sapply(1:length(ybar_i),
                    function(a) compute_surv_cor(ybar_i[a]*PFM,
                                                 PMF     = PFM,
                                                 PhiF    = phif_mark,
                                                 PhiM    = phim_mark))
  
  
  # Get mean of jackknife estimates
  gamma_i_mean <- mean(gamma_i)
  
  # Compute sum of squares and cubed of jacknife estimate of gamma 
  ss_jackknife <- sum((gamma_i_mean - gamma_i)^2)
  sc_jackknife <- sum((gamma_i_mean - gamma_i)^3)
  
  # Compute accelerated term for interval adjustment 
  a_hat <- sc_jackknife/(6 * (ss_jackknife)^(3/2))
  a_hat <- ifelse(is.nan(a_hat)|is.infinite(a_hat), 0, a_hat)
  
  # Get new percentile values 
  l <- pnorm(z0 + (z0 + qnorm(alpha/2))/(1-a_hat * (z0 + qnorm(alpha/2))))
  u <- pnorm(z0 + (z0 + qnorm((1-alpha/2)))/(1-a_hat * (z0 + qnorm((1-alpha/2)))))
  
  # 1-alpha interval
  out <- quantile(btrsp, c(l, u))
  names(out) <- c(paste0(as.character(100*alpha/2),"%"), paste0(as.character(100* (1-alpha/2)),"%"))
  return(out)
}

# Get Scenario Grid for correlation estimator simulation study
get_scenarios <- function(){
  
  # General Scenarios
  n     <- c(150, 250)
  k     <- c(15,  25)
  PF    <- c(0.45, 0.75)
  PM    <- c(0.45, 0.75)
  PhiF  <- c(0.8)
  PhiM  <- c(0.8)
  rho   <- sort(unique(c(0, 0.9, -0.1, seq(-0.1, 0.9, by = 0.15))))
  gamma <- sort(unique(c(0, 0.9, -0.1, seq(-0.1, 0.9, by = 0.15))))
  
  scenario_grid_base  <- expand.grid(n_obs    = n,
                                     k        = k,
                                     rho_true = rho,
                                     gam_true = gamma, 
                                     PhiF     = PhiF,
                                     PhiM     = PhiM,
                                     PF       = PF)
  
  scenario_grid_base$PM            <- scenario_grid_base$PF
  scenario_grid_base$Beta0         <- 1e3
  scenario_grid_base$Beta1         <- 1e3
  scenario_grid_base$Delta         <- 1
  scenario_grid_base$imputed_pairs <- TRUE
  
  # Testing Hypothesis with Alternative
  
  phi_param <- compute_jbin_param_cjs(0.7,0.8)
  p_param   <-compute_jbin_param_cjs(0.7,0.75)
  
  gl <- phi_param$cor_lower_bound
  gu <- phi_param$cor_upper_bound
  rl <- p_param$cor_lower_bound
  ru <- p_param$cor_upper_bound
  
  scenario_grid_alternative <- expand.grid(n_obs         = 250,
                                           k             = 25,
                                           rho_true      = rho[rho <= ru & rho >= rl],
                                           gam_true      = gamma[gamma <= gu & gamma >= gl],
                                           PhiF          = 0.7,
                                           PhiM          = 0.8,
                                           PF            = 0.70,
                                           PM            = 0.75,
                                           Beta0         = 1e3,
                                           Beta1         = 1e3,
                                           Delta         = 1,
                                           imputed_pairs = TRUE)
  
  # Hduck: Test1
  phi_param <- compute_jbin_param_cjs(0.67,0.74)
  p_param   <-compute_jbin_param_cjs(0.48,0.21)
  
  gl <- phi_param$cor_lower_bound
  gu <- phi_param$cor_upper_bound
  rl <- p_param$cor_lower_bound
  ru <- p_param$cor_upper_bound
  
  scenario_grid_hduck1 <- expand.grid(n_obs         = 250,
                                      k             = 25,
                                      rho_true      = rho[rho <= ru & rho >= rl],
                                      gam_true      = gamma[gamma <= gu & gamma >= gl],
                                      PhiF          = 0.67,
                                      PhiM          = 0.74,
                                      PF            = 0.48,
                                      PM            = 0.21,
                                      Beta0         = 1e3,
                                      Beta1         = 1e3,
                                      Delta         = 1,
                                      imputed_pairs = TRUE)
  
  # Hduck: Test2
  phi_param <- compute_jbin_param_cjs(0.67,0.74)
  p_param   <-compute_jbin_param_cjs(0.33,0.33)
  
  gl <- phi_param$cor_lower_bound
  gu <- phi_param$cor_upper_bound
  rl <- p_param$cor_lower_bound
  ru <- p_param$cor_upper_bound
  
  scenario_grid_hduck2 <- expand.grid(n_obs         = 250,
                                      k             = 25,
                                      rho_true      = rho[rho <= ru & rho >= rl],
                                      gam_true      = gamma[gamma <= gu & gamma >= gl],
                                      PhiF          = 0.67,
                                      PhiM          = 0.74,
                                      PF            = 0.33,
                                      PM            = 0.33,
                                      Beta0         = 1e3,
                                      Beta1         = 1e3,
                                      Delta         = 1,
                                      imputed_pairs = TRUE)
  
  # Hduck: Test3
  phi_param <- compute_jbin_param_cjs(0.7,0.7)
  p_param   <-compute_jbin_param_cjs(0.48, 0.21)
  
  gl <- phi_param$cor_lower_bound
  gu <- phi_param$cor_upper_bound
  rl <- p_param$cor_lower_bound
  ru <- p_param$cor_upper_bound
  
  scenario_grid_hduck3 <- expand.grid(n_obs         = 250,
                                      k             = 25,
                                      rho_true      = rho[rho <= ru & rho >= rl],
                                      gam_true      = gamma[gamma <= gu & gamma >= gl],
                                      PhiF          = 0.7,
                                      PhiM          = 0.7,
                                      PF            = 0.48,
                                      PM            = 0.21,
                                      Beta0         = 1e3,
                                      Beta1         = 1e3,
                                      Delta         = 1,
                                      imputed_pairs = TRUE)
  
  # Hduck: Test4
  phi_param <- compute_jbin_param_cjs(0.7,0.7)
  p_param   <-compute_jbin_param_cjs(0.33,0.33)
  
  gl <- phi_param$cor_lower_bound
  gu <- phi_param$cor_upper_bound
  rl <- p_param$cor_lower_bound
  ru <- p_param$cor_upper_bound
  
  scenario_grid_hduck4 <- expand.grid(n_obs         = 250,
                                      k             = 25,
                                      rho_true      = rho[rho <= ru & rho >= rl],
                                      gam_true      = gamma[gamma <= gu & gamma >= gl],
                                      PhiF          = 0.7,
                                      PhiM          = 0.7,
                                      PF            = 0.33,
                                      PM            = 0.33,
                                      Beta0         = 1e3,
                                      Beta1         = 1e3,
                                      Delta         = 1,
                                      imputed_pairs = TRUE)
  
  # Testing Imputed Mates Off
  scenario_imputed_off <- expand.grid(n_obs         = 250,
                                      k             = 25,
                                      rho_true      = c(-0.10,  0.00,  0.20,  0.50,  0.80,  0.90),
                                      gam_true      = c(-0.10,  0.00,  0.20,  0.50,  0.80,  0.90), 
                                      PhiF          = 0.8,
                                      PhiM          = 0.8,
                                      PF            = 0.75,
                                      PM            = 0.75,
                                      Beta0         = 0.25,
                                      Beta1         = 0,
                                      Delta         = 1,
                                      imputed_pairs = TRUE)
  
  
  
  scenario_grid          <- rbind(scenario_grid_base, 
                                  scenario_grid_hduck1,
                                  scenario_grid_hduck2,
                                  scenario_grid_hduck3,
                                  scenario_grid_hduck4,
                                  scenario_grid_alternative,
                                  scenario_imputed_off)
  scenario_grid$PropF    <- 0.5
  scenario_grid$scenario <- 1:nrow(scenario_grid)
  return(scenario_grid)
}


# Compact ifelse for NA
convert_na <- function(x,y=0){
  ifelse(is.na(x), y, x)
}

binom_test_gamma <- function(test_prob, success, trials){
  p_cutoff <- dbinom(success, trials, test_prob)
  p_state_space <- dbinom(0:trials, trials, test_prob)
  p_set <- p_state_space[p_state_space <= p_cutoff]
  return(sum(p_set))
}


# Execute one replicate of Simulation Study 
execute_iteration  <- function(iter,
                               bstrp_iter,
                               scenario,
                               PM,
                               PF,
                               PhiF,
                               PhiM,
                               gam_true,
                               rho_true,
                               n_pop,
                               k,
                               tau,
                               Betas,
                               Delta,
                               PropF,
                               imputed_pairs,
                               small_out){
  
  # Simulate Data ------------------------------------------------------------------------------------------------
  cat(paste0("Iteration#:", iter , " - Generating pair-swap mark-recapture data..."),"\n")
  
  # Store Seed information to be able to reproduce datasets if needed
  random_seed <- .Random.seed
  
  # Generate One set of pair-swap data
  ps_data <- simulate_ps_data(n             = n_pop,
                              k             = k,
                              propfemale    = PropF,
                              delta         = rep(Delta,k),
                              phif          = rep(PhiF,k),
                              phim          = rep(PhiM,k),
                              gam           = rep(gam_true,k),
                              pf            = rep(PF,k),
                              pm            = rep(PM,k),
                              rho           = rep(rho_true,k),
                              beta          = Betas,
                              tau           = tau,
                              imputed_pairs = imputed_pairs)
  cjs_data <- format_to_cjs(ps_data)
  
  # Effective Sample Size CJS Model
  n_eff    <- sum(colSums(cjs_data$x[,1:(cjs_data$k-1)]))
  
  # Compute ~true P and PHi
  PropF <- mean(cjs_data$female)
  PropM <- 1-PropF
  Phi <- PhiF * PropF + PhiM * PropM
  P   <- PF * PropF + PM * PropM
  
  # Gather True Param DF 
  ru_true            <- compute_jbin_param_cjs(PF,PM)$cor_lower_bound
  rl_true            <- compute_jbin_param_cjs(PF,PM)$cor_upper_bound
  gu_true            <- compute_jbin_param_cjs(PhiF,PhiM)$cor_upper_bound
  gl_true            <- compute_jbin_param_cjs(PhiF,PhiM)$cor_upper_bound
  true_params        <- c(PF,PM,PhiF, PhiM, P, Phi, gam_true, gl_true, gu_true, rho_true, rl_true, ru_true)
  param_names        <- c("PF","PM", "PhiF", "PhiM","P", "Phi", "gamma", "gl", "gu", "rho", "rl", "ru")
  true_param_df      <- data.frame(Truth = true_params, Parameter = param_names)
  
  # Run CJS Model MARK -----------------------------------------------------------------------------------------
  cat(paste0("Iteration#:", iter , " - Computing standard CJS estimates with program MARK..."),"\n")
  
  cjs_list <- list()
  versions <- c("B","S","R","N")
  
  for(i in 1:length(versions)){
    cat(paste0("Iteration#:", iter , " - Computing CJS estimates for Version:", versions[i], " ..."),"\n")
    cjs_list[[i]] <- run_cjs_model_marked(cjs_data = cjs_data,
                                          version  = versions[i])  
  }
  
  cjs_out <- do.call(rbind, cjs_list) %>% 
    left_join(true_param_df,
              by = "Parameter") %>% 
    mutate(Bias     = Truth - Est,
           In95     = 1*(Truth <= UB & Truth >= LB),
           Range95  = UB - LB,
           iter     = iter,
           scenario = scenario) 
  
  # Get Predicted Probs for Correlation Estimators (using full model version)
  phim_mark <- cjs_out %>% filter(Version == "B" & Parameter == "PhiM") %>% pull(Est)
  phif_mark <- cjs_out %>% filter(Version == "B" & Parameter == "PhiF") %>% pull(Est)
  pm_mark   <- cjs_out %>% filter(Version == "B" & Parameter == "PM")   %>% pull(Est)
  pf_mark   <- cjs_out %>% filter(Version == "B" & Parameter == "PF")   %>% pull(Est)
  #-------------------------------------------------------------------------------------------------------------
  
  # Compute Recapture Correlation Estimate----------------------------------------------------------------------
  
  # 1. Likelihood Approach -------------------------------------------------------------------------------------
  cat(paste0("Iteration#:", iter ," - Estimating recapture correlation, rho, using likelihood approach..."),"\n")
  rho_list <- compute_recapture_correlation(ps_data = ps_data, 
                                            PF      = pf_mark,
                                            PM      = pm_mark,
                                            model   = "likelihood")
  
  rho        <- rho_list$rho
  n_eff_rho  <- rho_list$n_eff_rho
  names(rho) <- "Est"
  
  obs_rho <- c(sum(rho_list$joint_recap == 1), 
               sum(rho_list$joint_recap == 2),
               sum(rho_list$joint_recap == 3),
               sum(rho_list$joint_recap == 4))  
  # prob_rho <- as.vector(rev(unlist(compute_jbin_cjs(pf_mark, pm_mark, 0))))
  # size <- sum(obs_rho)
  # groups <- length(obs_rho)
  # numEvents <- choose(size + groups - 1, groups - 1)
  # pval0_rho2 <- 0 #@ EMT::ExactMultinomialTestChisquare(obs_rho,prob_rho, size, groups, numEvents)$p.value
  
  #-------------------------------------------------------------------------------------------------------------
  
  # 2. Pearson Full Approach -----------------------------------------------------------------------------------
  # cat(paste0("Iteration#:", iter ," - Estimating recapture correlation rho using full pearson ..."),"\n")
  # pearson_rho <- compute_recapture_correlation(ps_data = ps_data, 
  #                                              PF      = pf_mark,
  #                                              PM      = pm_mark,
  #                                              model   = "full_pearson")$rho
  # names(pearson_rho) <- "Est"
  # 
  # # 2. Pearson Conditional Approach ----------------------------------------------------------------------------
  # cat(paste0("Iteration#:", iter ," - Estimating recapture correlation rho using pseudo-pearson..."),"\n")
  # pearson_partial_rho <- compute_recapture_correlation(ps_data = ps_data, 
  #                                                      PF      = pf_mark,
  #                                                      PM      = pm_mark,
  #                                                      model   = "partial_pearson")$rho
  # names(pearson_partial_rho) <- "Est"
  # 
  # Compute Survival Correlation Estimate-----------------------------------------------------------------------
  
  # Conditional on Rho from Likelihood Approach ----------------------------------------------------------------
  cat(paste0("Iteration#:", iter ," - Estimating survival correlation gamma|rho-likelihood..."),"\n")
  # If estimate of rho fails (no valid observations) pass dummy values of 10
  if(is.na(rho)|rho == 10){
    gamma <- 0
    n_eff_gamma <- 0
    ybar <- 0
    pval0_gamma <- 0
  } else {
    # Estimate Gamma from observed data
    gamma_list <- compute_survival_correlation(ps_data = ps_data,
                                               PFM     = compute_jbin_cjs(prob.f = pf_mark,
                                                                          prob.m = pm_mark,
                                                                          corr   = rho)$prob.mf,
                                               PhiF    = phif_mark,
                                               PhiM    = phim_mark)
    
    # Get Gamma and Effective Sample Size
    gamma             <- gamma_list$gamma
    n_eff_gamma       <- gamma_list$n_eff_gamma 
    n_success         <- gamma_list$n_success
    ybar              <- gamma_list$ybar
    # pval0_gamma2      <- 0 #binom_test_gamma(test_prob = compute_jbin_cjs(prob.f = pf_mark,
                                          #                              prob.m = pm_mark,
                                          #                              corr   = rho)$prob.mf * phif_mark * phim_mark, 
                                          # success   = n_success, 
                                          # trials    = n_eff_gamma)
  }
  #-------------------------------------------------------------------------------------------------------------
  
  # # Conditional on Rho from Likelihood Approach ----------------------------------------------------------------
  # cat(paste0("Iteration#:", iter ," - Estimating survival correlation gamma|rho-pearson..."),"\n")
  # 
  # # If estimate of rho fails (no valid observations) pass dummy values of 10
  # if(is.na(pearson_rho)|pearson_rho == 10){
  #   pearson_gamma <- 10
  # } else {
  #   
  #   # Estimate Gamma from observed data
  #   pearson_gamma <- compute_survival_correlation(ps_data = ps_data,
  #                                                 PFM     = compute_jbin_cjs(prob.f = pf_mark,
  #                                                                            prob.m = pm_mark,
  #                                                                            corr   = pearson_rho)$prob.mf,
  #                                                 PhiF    = phif_mark,
  #                                                 PhiM    = phim_mark)$gamma
  #   
  #   pearson_pval0_gamma <- binom_test_gamma(test_prob = compute_jbin_cjs(prob.f = pf_mark,
  #                                                                        prob.m = pm_mark,
  #                                                                        corr   = pearson_rho)$prob.mf * phif_mark * phim_mark, 
  #                                           success   = round(ybar*n_eff_gamma,0), 
  #                                           trials    = n_eff_gamma)
  #   
  # }
  # 
  # # Collect Results
  # names(pearson_gamma) <- "Est"
  # #-------------------------------------------------------------------------------------------------------------
  # 
  # # Conditional on Rho from Likelihood Approach ----------------------------------------------------------------
  # cat(paste0("Iteration#:", iter , " - Estimating survival correlation gamma|rho-psuedo-pearson..."),"\n")
  # 
  # # If estimate of rho fails (no valid observations) pass dummy values of 10
  # if(is.na(pearson_partial_rho)|pearson_partial_rho == 10){
  #   pearson_partial_gamma <- 10
  # } else {
  #   
  #   # Estimate Gamma from observed data
  #   pearson_partial_gamma <- compute_survival_correlation(ps_data = ps_data,
  #                                                         PFM     = compute_jbin_cjs(prob.f = pf_mark,
  #                                                                                    prob.m = pm_mark,
  #                                                                                    corr   = pearson_partial_rho)$prob.mf,
  #                                                         PhiF    = phif_mark,
  #                                                         PhiM    = phim_mark)$gamma
  #   
  #   pearson_partial_pval0_gamma <- binom_test_gamma(test_prob = compute_jbin_cjs(prob.f = pf_mark,
  #                                                                                prob.m = pm_mark,
  #                                                                                corr   = pearson_partial_rho)$prob.mf * phif_mark * phim_mark, 
  #                                                   success   = round(ybar*n_eff_gamma,0), 
  #                                                   trials    = n_eff_gamma)
  #   
  # }
  # 
  # # Collect Results
  # names(pearson_partial_gamma) <- "Est"
  #-------------------------------------------------------------------------------------------------------------
  
  # Run Bootstraps
  
  # Run Bootstrap Estimates for Likelihood Approach
  summ_corr_lik_list <- compute_full_bootstrap(iterations    = bstrp_iter,
                                               full_result   = !small_out,
                                               PM            = pm_mark,
                                               PF            = pf_mark,
                                               PhiF          = phif_mark,
                                               PhiM          = phim_mark, 
                                               gamma         = gamma,
                                               rho           = rho,
                                               n_pop         = n_pop,
                                               k             = k,
                                               tau           = tau,
                                               Betas         = Betas,
                                               Delta         = Delta,
                                               PropF         = PropF,
                                               imputed_pairs = imputed_pairs,
                                               method_rho    = "likelihood",
                                               test0rho = FALSE,
                                               test0gam = FALSE)
  
  summ_corr_lik <- summ_corr_lik_list$results
  
  
  # Get Pval for H0: Rho = 0 vs Ha rho != 0 (conditional on all other parameters)
  pval0_rho2 <- compute_full_bootstrap(iterations    = bstrp_iter,
                                       full_result   = FALSE,
                                       PM            = pm_mark,
                                       PF            = pf_mark,
                                       PhiF          = phif_mark,
                                       PhiM          = phim_mark, 
                                       gamma         = gamma,
                                       rho           = rho,
                                       n_pop         = n_pop,
                                       k             = k,
                                       tau           = tau,
                                       Betas         = Betas,
                                       Delta         = Delta,
                                       PropF         = PropF,
                                       imputed_pairs = imputed_pairs,
                                       method_rho    = "likelihood",
                                       test0rho      = TRUE,
                                       test0gam      = FALSE)$results$pval
  
  
  # Get Pval for H0: Rho = 0 vs Ha rho != 0 (conditional on all other parameters)
  pval0_gamma2 <- compute_full_bootstrap(iterations      = bstrp_iter,
                                         full_result     = FALSE,
                                         PM              = pm_mark,
                                         PF              = pf_mark,
                                         PhiF            = phif_mark,
                                         PhiM            = phim_mark, 
                                         gamma           = gamma,
                                         rho             = rho,
                                         n_pop           = n_pop,
                                         k               = k,
                                         tau             = tau,
                                         Betas           = Betas,
                                         Delta           = Delta,
                                         PropF           = PropF,
                                         imputed_pairs   = imputed_pairs,
                                         method_rho      = "likelihood",
                                         test0rho        = FALSE,
                                         test0gam        = TRUE)$results$pval
  
  # # Run Bootstrap Estimates for Pearson Approach
  # summ_corr_fp_list <- compute_full_bootstrap(iterations    = bstrp_iter,
  #                                             full_result   = !small_out,
  #                                             PM            = pm_mark,
  #                                             PF            = pf_mark,
  #                                             PhiF          = phif_mark,
  #                                             PhiM          = phim_mark, 
  #                                             gamma         = pearson_gamma,
  #                                             rho           = pearson_rho,
  #                                             n_pop         = n_pop,
  #                                             k             = k,
  #                                             tau           = tau,
  #                                             Betas         = Betas,
  #                                             Delta         = Delta,
  #                                             PropF         = PropF,
  #                                             imputed_pairs = imputed_pairs,
  #                                             method_rho    = "full_pearson")
  # 
  # 
  # summ_corr_fp <- summ_corr_fp_list$results
  # 
  # # Run Bootstrap Estimates for Partial-Pearson Approach
  # summ_corr_pp_list <- compute_full_bootstrap(iterations    = bstrp_iter,
  #                                             full_result    = !small_out,
  #                                             PM            = pm_mark,
  #                                             PF            = pf_mark,
  #                                             PhiF          = phif_mark,
  #                                             PhiM          = phim_mark, 
  #                                             gamma         = pearson_partial_gamma,
  #                                             rho           = pearson_partial_rho,
  #                                             n_pop         = n_pop,
  #                                             k             = k,
  #                                             tau           = tau,
  #                                             Betas         = Betas,
  #                                             Delta         = Delta,
  #                                             PropF         = PropF,
  #                                             imputed_pairs = imputed_pairs,
  #                                             method_rho    = "partial_pearson")
  # 
  # summ_corr_pp <- summ_corr_pp_list$results
  
  # Return Results----------------------------------------------------------------------------------------------
  cat("Formatting output ...","\n")
  # Gather Correlation Results
  # summ_corr           <- as.data.frame(rbind(summ_corr_lik,
  #                                            summ_corr_fp, 
  #                                            summ_corr_pp))
  # 
  
  summ_corr           <- as.data.frame(rbind(summ_corr_lik))
  # summ_corr$Pval02    <- c(pval0_rho2, pval0_gamma2)#,0,pearson_pval0_gamma,0,pearson_partial_pval0_gamma)
  summ_corr$Pval02    <- c(pval0_rho2, pval0_gamma2)
  summ_corr           <- summ_corr[,c("Parameter","Est","Est_Btstrp","SE", 
                                      "2.5%","25%","50%","75%","97.5%", 
                                      "Rho_Estimator", "Pval0","Pval02","num_failures")]
  summ_corr           <- summ_corr %>% 
    left_join(true_param_df, by = "Parameter") %>% 
    mutate(Est = as.numeric(Est),
           Est_Btstrp = as.numeric(Est_Btstrp),
           SE = as.numeric(SE),
           `2.5%` = as.numeric(`2.5%`),
           `25%` = as.numeric(`25%`),
           `50%` = as.numeric(`50%`),
           `75%` = as.numeric(`75%`),
           `97.5%` = as.numeric(`97.5%`),
           Pval0  = as.numeric(Pval0),
           Pval02  = as.numeric(Pval02),
           num_failures = as.numeric(num_failures),
           Bias                = Truth - Est,
           Bias_Btstrp1_Mean   = Truth - Est_Btstrp,
           Bias_Btstrp2_Mean   = Est   - Est_Btstrp,
           Bias_Btstrp1_Median = Truth - `50%`,
           Bias_Btstrp2_Median = Est   - `50%`,
           In95                = 1*(Truth <= `97.5%` & Truth >= `2.5%`),
           In50                = 1*(Truth <= `75%`   & Truth  >= `25%`),
           Range95             = `97.5%` - `2.5%`,
           Range50             = `75%` - `25%`,
           Cover095            = 1*(0 <= `97.5%` & 0 >= `2.5%`),    
           Cover050            = 1*(0 <= `75%`   & 0 >= `25%`),    
           iter                = iter,
           scenario            = scenario) 
  
  # Compute CChat Adjustment
  summ_chat <- data.frame(Method     = c("Likelihood", 
                                         # "Pearson", "Psuedo-Pearson",
                                         "Truth"),
                          CChatRho   = c(compute_proposed_chat(rho, 0), 
                                         #compute_proposed_chat(pearson_rho, 0), 
                                         #compute_proposed_chat(pearson_partial_rho, 0),
                                         compute_proposed_chat(rho_true, 0)),
                          CChatGamma = c(compute_proposed_chat(0, gamma), 
                                         #compute_proposed_chat(0, pearson_gamma), 
                                         #compute_proposed_chat(0, pearson_partial_gamma),
                                         compute_proposed_chat(0, gam_true)),
                          CChat      = c(compute_proposed_chat(rho, gamma), 
                                         #compute_proposed_chat(pearson_rho, pearson_gamma), 
                                         #compute_proposed_chat(pearson_partial_rho, pearson_partial_gamma),
                                         compute_proposed_chat(rho_true, gam_true)),
                          iter       = rep(iter, 2),
                          scenario   = rep(scenario, 2))
  
  # Add Correlation Chat Correction to Mark-Recapture Results
  summ_cjs <- cjs_out %>% 
    mutate(CChatAdj_Likelihood      = ifelse(Parameter == "P", 
                                             compute_proposed_chat(rho, 0), 
                                             ifelse(Parameter == "Phi",
                                                    compute_proposed_chat(0, gamma), 
                                                    1)),
           # CChatAdj_Pearson         = ifelse(Parameter == "P", 
           #                                   compute_proposed_chat(pearson_rho, 0), 
           #                                   ifelse(Parameter == "Phi",
           #                                          compute_proposed_chat(0, pearson_gamma), 
           #                                          1)),
           # CChatAdj_Partial_Pearson = ifelse(Parameter == "P", 
           #                                   compute_proposed_chat(pearson_partial_rho, 0), 
           #                                   ifelse(Parameter == "Phi",
           #                                          compute_proposed_chat(0, pearson_partial_gamma), 
           #                                          1)),
           CChatAdj_Truth           = ifelse(Parameter == "P", 
                                             compute_proposed_chat(rho_true, 0), 
                                             ifelse(Parameter == "Phi",
                                                    compute_proposed_chat(0, gam_true), 
                                                    1)),
           UBAdj_Likelihood         = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Likelihood), alpha = 0.05)[["ub"]],
           LBAdj_Likelihood         = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Likelihood), alpha = 0.05)[["lb"]],
           # UBAdj_Pearson            = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Pearson), alpha = 0.05)[["ub"]],
           # LBAdj_Pearson            = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Pearson), alpha = 0.05)[["lb"]],
           # UBAdj_Partial_Pearson    = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Partial_Pearson), alpha = 0.05)[["ub"]],
           # LBAdj_Partial_Pearson    = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Partial_Pearson), alpha = 0.05)[["lb"]],
           In95_Likelihood          = 1*(Truth <= UBAdj_Likelihood & Truth >= LBAdj_Likelihood),
           # In95_Pearson             = 1*(Truth <= UBAdj_Pearson & Truth >= LBAdj_Pearson),
           # In95_Partial_Pearson     = 1*(Truth <= UBAdj_Partial_Pearson & Truth >= LBAdj_Partial_Pearson),
           Range95_Likelihood       = UBAdj_Likelihood - LBAdj_Likelihood#,
           # Range95_Pearson          = UBAdj_Pearson - LBAdj_Pearson,
           # Range95_Partial_Pearson  = UBAdj_Partial_Pearson - LBAdj_Partial_Pearson
    ) 
  
  # Produce AICC comparisons
  summ_aic <- compute_aic_summ(summ_cjs = summ_cjs,
                               n_eff    = n_eff,
                               iter     = iter,
                               scenario = scenario)
  
  # Summarize Sample Sizes
  summ_n <- data.frame(n_eff                       = n_eff,
                       n_eff_rho                   = n_eff_rho,
                       n_eff_gamma                 = n_eff_gamma,
                       iter                        = iter,
                       scenario                    = scenario)
  
  cat("Success, returning results ...","\n")
  
  if(small_out){
    results <- list(random_seed                 = random_seed,
                    summ_n                      = summ_n,
                    summ_cjs                    = summ_cjs,
                    summ_corr                   = summ_corr,
                    summ_chat                   = summ_chat,
                    summ_aic                    = summ_aic)
  } else {
    results <- list(random_seed                 = random_seed,
                    summ_n                      = summ_n,
                    summ_cjs                    = summ_cjs,
                    summ_corr                   = summ_corr,
                    summ_chat                   = summ_chat,
                    summ_aic                    = summ_aic,
                    bstrp_lik_rep               = summ_corr_lik_list$replicates)#,
    # bstrp_fp_rep                = summ_corr_fp_list$replicates,
    #  bstrp_pp_rep                = summ_corr_pp_list$replicates)
  }
  
  
  return(results)
}


run_bootstrap_replicate <- function(PM,
                                    PF,
                                    PhiF,
                                    PhiM,
                                    gamma,
                                    rho,
                                    n_pop,
                                    k,
                                    tau,
                                    Betas,
                                    Delta,
                                    PropF,
                                    imputed_pairs,
                                    method_rho){
  
  
  
  # Generate One set of pair-swap data
  ps_data <- simulate_ps_data(n             = n_pop,
                              k             = k,
                              propfemale    = PropF,
                              delta         = rep(Delta,k),
                              phif          = rep(PhiF,k),
                              phim          = rep(PhiM,k),
                              gam           = rep(gamma,k),
                              pf            = rep(PF,k),
                              pm            = rep(PM,k),
                              rho           = rep(rho,k),
                              beta          = Betas,
                              tau           = tau,
                              imputed_pairs = imputed_pairs)
  cjs_data <- format_to_cjs(ps_data)
  # Run CJS Model MARK -----------------------------------------------------------------------------------------
  
  # Process CJS Data for Mark Version B
  dat.process <- cjs_data$x %>% 
    as.data.frame() %>% 
    unite("ch",sep="") %>% 
    mutate(sex = as.factor(ifelse(cjs_data$female == 1, "F", "M"))) %>%
    arrange(sex)
  
  mark.process <- suppressWarnings(suppressMessages(marked::process.data(data   = dat.process,
                                                                         model  = "CJS",
                                                                         groups = "sex")))
  
  model.parameters <- list(Phi = list(formula = ~sex), p = list(formula = ~sex))
  mark.ddl <-  suppressWarnings(suppressMessages(marked::make.design.data(mark.process)))
  #Generate Estimates
  cjs_out <-  suppressWarnings(suppressMessages(crm2(data             = mark.process,
                                                     ddl              = mark.ddl,
                                                     model.parameters = model.parameters)))
  
  # Get Predicted Probs for Correlation Estimators (using full model version)
  phim_boot <- cjs_out[["PhiM"]]
  phif_boot <- cjs_out[["PhiF"]]
  pm_boot   <- cjs_out[["PM"]]
  pf_boot   <- cjs_out[["PF"]]
  #-------------------------------------------------------------------------------------------------------------
  
  # Compute Recapture Correlation Estimate----------------------------------------------------------------------
  rho_boot <- compute_recapture_correlation(ps_data = ps_data, 
                                            PF      = pf_boot,
                                            PM      = pm_boot,
                                            model   = method_rho)$rho
  names(rho_boot) <- "Est"
  
  #-------------------------------------------------------------------------------------------------------------
  
  # Compute Survival Correlation Estimate-----------------------------------------------------------------------
  
  # Conditional on Rho from Likelihood Approach ----------------------------------------------------------------
  # If estimate of rho fails (no valid observations) pass dummy values of 10
  if(is.na(rho_boot)|rho_boot == 10){
    gamma_boot <- 10
  } else {
    
    # Estimate Gamma from observed data
    gamma_boot <- compute_survival_correlation(ps_data = ps_data,
                                               PFM     = compute_jbin_cjs(prob.f = pf_boot,
                                                                          prob.m = pm_boot,
                                                                          corr   = rho_boot)$prob.mf,
                                               PhiF    = phif_boot,
                                               PhiM    = phim_boot)$gamma
    names(gamma_boot) <- "Est"
  }
  
  #-------------------------------------------------------------------------------------------------------------
  return(c(rho_boot,gamma_boot))
}

compute_full_bootstrap <- function(iterations,
                                   full_result,
                                   PM,
                                   PF,
                                   PhiF,
                                   PhiM,
                                   gamma,
                                   rho,
                                   n_pop,
                                   k,
                                   tau,
                                   Betas,
                                   Delta,
                                   PropF,
                                   imputed_pairs,
                                   method_rho,
                                   test0rho = FALSE,
                                   test0gam = FALSE){
  
  
  # Test if rho = 0
  if(test0rho){
    
    replicates <- lapply(1:iterations,
                         function(x) run_bootstrap_replicate(PM = PM,
                                                             PF = PF,
                                                             PhiF = PhiF,
                                                             PhiM = PhiM,
                                                             gamma = gamma,
                                                             rho = 0,
                                                             n_pop = n_pop,
                                                             k = k,
                                                             tau = tau,
                                                             Betas = Betas,
                                                             Delta = Delta,
                                                             PropF = PropF,
                                                             imputed_pairs = imputed_pairs,
                                                             method_rho = method_rho))
    
    
    replicates <- do.call(rbind, replicates)
    pval       <- mean(1 * (abs(replicates[,1]-0) > abs(rho - 0)))
    se         <- sd(replicates[,1])
    results    <- data.frame(Parameter = "rho",pval = pval, se = se)
    
    # Test if gamma = 0
  } else if(test0gam){
    
    replicates <- lapply(1:iterations,
                         function(x) run_bootstrap_replicate(PM = PM,
                                                             PF = PF,
                                                             PhiF = PhiF,
                                                             PhiM = PhiM,
                                                             gamma = 0,
                                                             rho = rho,
                                                             n_pop = n_pop,
                                                             k = k,
                                                             tau = tau,
                                                             Betas = Betas,
                                                             Delta = Delta,
                                                             PropF = PropF,
                                                             imputed_pairs = imputed_pairs,
                                                             method_rho = method_rho))
    
    replicates <- do.call(rbind, replicates)
    pval       <- mean(1 * (abs(replicates[,2]-0) > abs(gamma - 0)))
    se         <- sd(replicates[,2])
    results    <- data.frame(Parameter = "gamma",pval = pval, se = se)
    
    # Get Confidence Intervals around estimates
  } else {
    replicates <- lapply(1:iterations,
                         function(x) run_bootstrap_replicate(PM = PM,
                                                             PF = PF,
                                                             PhiF = PhiF,
                                                             PhiM = PhiM,
                                                             gamma = gamma,
                                                             rho = rho,
                                                             n_pop = n_pop,
                                                             k = k,
                                                             tau = tau,
                                                             Betas = Betas,
                                                             Delta = Delta,
                                                             PropF = PropF,
                                                             imputed_pairs = imputed_pairs,
                                                             method_rho = method_rho))
    
    replicates <- do.call(rbind, replicates)
    
    results <- as.data.frame(rbind(compute_btsrp_summary(rho,   replicates[,1], method_rho), 
                                   compute_btsrp_summary(gamma, replicates[,2], method_rho))) 
    
    results$Parameter <- c("rho","gamma")
  }
  
  if(full_result){
    return(list(results = results, replicates = replicates))
  } else {
    return(list(results = results))
  }
  
}

# Summarize Bootstrap output
compute_btsrp_summary <- function(estimate, 
                                  bstrp_replicates, 
                                  method_rho){
  # Non-Parametric conditional estimator results
  n                     <- length(bstrp_replicates)
  bstrp_replicates      <- bstrp_replicates[bstrp_replicates != 10]
  n_failures            <- n - length(bstrp_replicates)
  names(n_failures)     <- "num_failures"
  mean_bstrp            <- mean(bstrp_replicates)
  names(mean_bstrp)     <- "Est_Btstrp"
  se_bstrp              <- sd(bstrp_replicates)
  names(se_bstrp)       <- "SE"
  quantiles_btstrp      <- quantile(bstrp_replicates, c(0.025, 0.25, 0.5, 0.75, 0.975))
  method_rho            <- method_rho
  names(method_rho)     <- "Rho_Estimator"
  Pval0                 <- mean(1 * (abs(bstrp_replicates - estimate) > abs(estimate -0)))
  names(Pval0)          <- "Pval0"
  names(estimate)       <- "Est"
  summ_btstrp           <- c(estimate, 
                             mean_bstrp, 
                             se_bstrp,
                             quantiles_btstrp,
                             method_rho,
                             Pval0,
                             n_failures)
  
  return(summ_btstrp)
} 

# Strip down CRM function to produce only CJS estimates to speed up bootstrap
crm2 <- function (data, ddl = NULL, begin.time = 1, model = "CJS", title = "", 
                  model.parameters = list(), design.parameters = list(), initial = NULL, 
                  groups = NULL, time.intervals = NULL, debug = FALSE, method = NULL, 
                  hessian = FALSE, accumulate = TRUE, chunk_size = 1e+07, control = list(), 
                  refit = 1, itnmax = 5000, scale = NULL, run = TRUE, burnin = 100, 
                  iter = 1000, use.admb = FALSE, use.tmb = FALSE, crossed = NULL, 
                  reml = FALSE, compile = FALSE, extra.args = NULL, strata.labels = NULL, 
                  clean = NULL, save.matrices = FALSE, simplify = FALSE, getreals = FALSE, 
                  real.ids = NULL, check = FALSE, prior = FALSE, prior.list = NULL, 
                  useHess = FALSE, optimize = TRUE, ...){
  
  data.proc = data
  model = data$model
  number.of.groups = 2
  par.list =  marked::setup.parameters(data.proc$model, check = TRUE)
  parameters =  marked::setup.parameters(data.proc$model, model.parameters, 
                                         data$nocc, number.of.groups = number.of.groups)
  parameters = parameters[par.list]
  re = FALSE
  
  for (i in 1:length(parameters)) {
    if (is.null(parameters[[i]]$formula)) 
      parameters[[i]]$formula = ~1
  }
  
  design.parameters = ddl$design.parameters
  ddl = marked:::set.fixed(ddl, parameters)
  fullddl = ddl
  dml = marked:::create.dml(ddl, model.parameters = parameters, design.parameters = design.parameters, 
                            chunk_size = chunk_size, simplify = simplify, use.admb = use.admb)
  initial.list = NULL
  fulldml = dml
  
  nocc = data.proc$nocc
  if (!is.null(ddl$Phi$time.interval)){
    time.intervals = matrix(ddl$Phi$time.interval, nrow(x$data), 
                            ncol = nocc - 1, byrow = TRUE)
  } else if(is.vector(data.proc$time.intervals)){
    time.intervals = matrix(data.proc$time.intervals, nrow = nrow(data.proc$data), 
                            ncol = nocc - 1, byrow = TRUE)
  } else {
    time.intervals = data.proc$time.intervals
  }
  parameters = marked:::create.fixed.matrix(ddl, parameters)
  x = data.proc$data
  ch = x$ch
  freq = NULL
  if (!is.null(x$freq)) 
    freq = x$freq
  imat = process.ch(ch, freq, all = FALSE)
  if (is.null(initial)){
    par = marked:::cjs.initial(dml, imat)
  } else {
    par = marked:::set.initial(names(dml), dml, initial)$par
  } 
  
  initial = par
  model_data = list(Phi.dm = dml$Phi$fe, p.dm = dml$p$fe, imat = imat, 
                    Phi.fixed = parameters$Phi$fixed, p.fixed = parameters$p$fixed, 
                    time.intervals = time.intervals)
  
  model_data.save = model_data
  model_data = marked:::cjs.accumulate(x, model_data, nocc, freq, 
                                       chunk_size = chunk_size)
  scale = marked:::set.scale(names(dml), model_data, scale)
  model_data = marked:::scale.dm(model_data, scale)
  par = marked:::scale.par(par, scale)
  
  markedfunc_eval = 0
  cjsenv = environment()
  par = optim(par, marked:::cjs.lnl, model_data = model_data, 
              method = method, hessian = FALSE, debug = debug, 
              control = control, cjsenv = cjsenv)$par
  
  names(par) <- c("PhiF","PhiM","PF","PM")
  par[2] <- par[1] + par[2]
  par[4] <- par[3] + par[4]
  par <- boot::inv.logit(par)
  
  return(par)
  
}

# Execute entire simulation study for correlation estimators
execute_simulation <- function(niter,
                               bstrp_iter,
                               scenario,
                               PM,
                               PF,
                               PhiF,
                               PhiM,
                               gamma,
                               rho,
                               n_pop,
                               k,
                               tau,
                               Betas,
                               Delta,
                               PropF,
                               imputed_pairs,
                               method_rho,
                               small_out     = TRUE){
  # Run niter replicates 
  results_list <- lapply(1:niter, function(iter) execute_iteration(iter          = iter,
                                                                   bstrp_iter    = bstrp_iter,
                                                                   scenario      = scenario,
                                                                   PM            = PM,
                                                                   PF            = PF,
                                                                   PhiF          = PhiF,
                                                                   PhiM          = PhiM,
                                                                   gam_true      = gamma,
                                                                   rho_true      = rho,
                                                                   n_pop         = n_pop,
                                                                   k             = k,
                                                                   Betas         = Betas,
                                                                   Delta         = Delta,
                                                                   PropF         = PropF,
                                                                   tau           = tau,
                                                                   imputed_pairs = imputed_pairs,
                                                                   small_out     = small_out))
  
  # Construct Summaries of Probs and Correlations
  summary_corr    <- do.call(rbind, lapply(1:niter, function(iter) results_list[[iter]]$summ_corr))
  summary_cjs     <- do.call(rbind, lapply(1:niter, function(iter) results_list[[iter]]$summ_cjs))
  summ_chat       <- do.call(rbind, lapply(1:niter, function(iter) results_list[[iter]]$summ_chat))
  summ_n          <- do.call(rbind, lapply(1:niter, function(iter) results_list[[iter]]$summ_n))
  summ_aic        <- do.call(rbind, lapply(1:niter, function(iter) results_list[[iter]]$summ_aic))
  random_seeds    <- lapply(1:niter, function(iter) results_list[[iter]]$random_seed)
  
  # Return Results
  if(small_out){
    # Only summaries and random.seeds
    out <- list(summary_corr   = summary_corr,
                summary_cjs    = summary_cjs,
                summ_chat      = summ_chat,
                summ_n         = summ_n,
                summ_aic       = summ_aic,
                random_seeds   = random_seeds)
  } else {
    # Summaries, data, bootstrap estimates
    out <- list(results_list   = results_list,
                summary_corr   = summary_corr,
                summary_cjs    = summary_cjs,
                summ_chat      = summ_chat,
                summ_n         = summ_n,
                summ_aic       = summ_aic,
                random_seeds   = random_seeds)
  }
  
  return(out)
  
}


# Execute one replicate of Simulation Study 
execute_application  <- function(ps_data, 
                                 bstrp_iter,
                                 small_out = F){
  
  # Simulate Data ------------------------------------------------------------------------------------------------
  
  
  # Store Seed information to be able to reproduce datasets if needed
  random_seed <- .Random.seed
  
  cjs_data <- format_to_cjs(ps_data)
  
  # Effective Sample Size CJS Model
  n_eff    <- sum(colSums(cjs_data$x[,1:(cjs_data$k-1)]))
  
  # Run CJS Model MARK -----------------------------------------------------------------------------------------
  cat(paste0("Computing standard CJS estimates with marked..."),"\n")
  
  cjs_list <- list()
  versions <- c("B","S","R","N")
  
  for(i in 1:length(versions)){
    cat(paste0("Computing CJS estimates for Version:", versions[i], " ..."),"\n")
    cjs_list[[i]] <- run_cjs_model_marked(cjs_data = cjs_data,
                                          version  = versions[i])  
  }
  
  cjs_out <- do.call(rbind, cjs_list) %>% 
    mutate(Range95  = UB - LB) 
  
  # Get Predicted Probs for Correlation Estimators (using full model version)
  phim_mark <- cjs_out %>% filter(Version == "B" & Parameter == "PhiM") %>% pull(Est)
  phif_mark <- cjs_out %>% filter(Version == "B" & Parameter == "PhiF") %>% pull(Est)
  pm_mark   <- cjs_out %>% filter(Version == "B" & Parameter == "PM")   %>% pull(Est)
  pf_mark   <- cjs_out %>% filter(Version == "B" & Parameter == "PF")   %>% pull(Est)
  #-------------------------------------------------------------------------------------------------------------
  
  # Compute Recapture Correlation Estimate----------------------------------------------------------------------
  
  # 1. Likelihood Approach -------------------------------------------------------------------------------------
  cat(paste0("Estimating recapture correlation, rho, using likelihood approach..."),"\n")
  rho_list <- compute_recapture_correlation(ps_data = ps_data, 
                                            PF      = pf_mark,
                                            PM      = pm_mark,
                                            model   = "likelihood")
  
  rho        <- rho_list$rho
  n_eff_rho  <- rho_list$n_eff_rho
  names(rho) <- "Est"
  
  #-------------------------------------------------------------------------------------------------------------
  
  # # 2. Pearson Full Approach -----------------------------------------------------------------------------------
  # cat(paste0("Estimating recapture correlation rho using full pearson ..."),"\n")
  # pearson_rho <- compute_recapture_correlation(ps_data = ps_data, 
  #                                              PF      = pf_mark,
  #                                              PM      = pm_mark,
  #                                              model   = "full_pearson")$rho
  # names(pearson_rho) <- "Est"
  # 
  # # 2. Pearson Conditional Approach ----------------------------------------------------------------------------
  # cat(paste0("Estimating recapture correlation rho using pseudo-pearson..."),"\n")
  # pearson_partial_rho <- compute_recapture_correlation(ps_data = ps_data, 
  #                                                      PF      = pf_mark,
  #                                                      PM      = pm_mark,
  #                                                      model   = "partial_pearson")$rho
  # names(pearson_partial_rho) <- "Est"
  
  # Compute Survival Correlation Estimate-----------------------------------------------------------------------
  
  # Conditional on Rho from Likelihood Approach ----------------------------------------------------------------
  cat(paste0("Estimating survival correlation gamma|rho-likelihood..."),"\n")
  # If estimate of rho fails (no valid observations) pass dummy values of 10
  if(is.na(rho)|rho == 10){
    gamma <- 10
    gamma_bs_np <- rep(10, 1000)
    gamma_bs_sp <- rep(10, 1000)
    gamma_bs_null <- rep(10, 1000)
    n_eff_gamma <- 0
    ybar <- NA
  } else {
    # Estimate Gamma from observed data
    gamma_list <- compute_survival_correlation(ps_data = ps_data,
                                               PFM     = compute_jbin_cjs(prob.f = pf_mark,
                                                                          prob.m = pm_mark,
                                                                          corr   = rho)$prob.mf,
                                               PhiF    = phif_mark,
                                               PhiM    = phim_mark)
    
    # Get Gamma and Effective Sample Size
    gamma       <- gamma_list$gamma
    n_eff_gamma <- gamma_list$n_eff_gamma 
    ybar        <- gamma_list$ybar
  }
  #-------------------------------------------------------------------------------------------------------------
  
  # # Conditional on Rho from Likelihood Approach ----------------------------------------------------------------
  # cat(paste0("Estimating survival correlation gamma|rho-pearson..."),"\n")
  # 
  # # If estimate of rho fails (no valid observations) pass dummy values of 10
  # if(is.na(pearson_rho)|pearson_rho == 10){
  #   pearson_gamma <- 10
  #   pearson_gamma_bs_np <- rep(10, 1000)
  #   pearson_gamma_bs_sp <- rep(10, 1000)
  #   gamma_bs_pearson_null <- rep(10, 1000)
  # } else {
  #   
  #   # Estimate Gamma from observed data
  #   pearson_gamma <- compute_survival_correlation(ps_data = ps_data,
  #                                                 PFM     = compute_jbin_cjs(prob.f = pf_mark,
  #                                                                            prob.m = pm_mark,
  #                                                                            corr   = pearson_rho)$prob.mf,
  #                                                 PhiF    = phif_mark,
  #                                                 PhiM    = phim_mark)$gamma
  #   
  # }
  # 
  # # Collect Results
  # names(pearson_gamma) <- "Est"
  # #-------------------------------------------------------------------------------------------------------------
  # 
  # # Conditional on Rho from Likelihood Approach ----------------------------------------------------------------
  # cat(paste0("Estimating survival correlation gamma|rho-psuedo-pearson..."),"\n")
  # 
  # # If estimate of rho fails (no valid observations) pass dummy values of 10
  # if(is.na(pearson_partial_rho)|pearson_partial_rho == 10){
  #   pearson_partial_gamma <- 10
  #   pearson_partial_gamma_bs_np <- rep(10, 1000)
  #   pearson_partial_gamma_bs_sp <- rep(10, 1000)
  #   gamma_bs_pearson_partial_null <- rep(10, 1000)
  # } else {
  #   
  #   # Estimate Gamma from observed data
  #   pearson_partial_gamma <- compute_survival_correlation(ps_data = ps_data,
  #                                                         PFM     = compute_jbin_cjs(prob.f = pf_mark,
  #                                                                                    prob.m = pm_mark,
  #                                                                                    corr   = pearson_partial_rho)$prob.mf,
  #                                                         PhiF    = phif_mark,
  #                                                         PhiM    = phim_mark)$gamma
  #   
  # }
  # 
  # # Collect Results
  # names(pearson_partial_gamma) <- "Est"
  #-------------------------------------------------------------------------------------------------------------
  
  # Run Bootstraps
  
  tau <- prop.table(table(cjs_data$initial_entry))
  tk <- names(tau)
  tau <- unname(tau)
  tau <- as.numeric(tau)
  
  df_tau <- data.frame(t=as.numeric(tk), tau1 = tau)
  df_zero <- data.frame(t = 1:(cjs_data$k), tau2 = 0)
  tau <- full_join(df_tau, df_zero, by = "t") %>% 
    replace_na(list(tau1 = 0, tau2 = 0)) %>% 
    mutate(tau = tau1+tau2) %>% pull(tau)
  
  # Run Bootstrap Estimates for Likelihood Approach
  summ_corr_lik_list <- compute_full_bootstrap(iterations    = bstrp_iter,
                                               full_result   = !small_out,
                                               PM            = pm_mark,
                                               PF            = pf_mark,
                                               PhiF          = phif_mark,
                                               PhiM          = phim_mark, 
                                               gamma         = gamma,
                                               rho           = rho,
                                               n_pop         = cjs_data$n,
                                               k             = cjs_data$k,
                                               tau           = tau,
                                               Betas         = c(1000,1000),
                                               Delta         = 1.0,
                                               PropF         = mean(cjs_data$female),
                                               imputed_pairs = TRUE,
                                               method_rho    = "likelihood",
                                               test0rho      = FALSE,
                                               test0gam      = FALSE)
  
  summ_corr_lik <- summ_corr_lik_list$results
  
  # Get Pval for H0: Rho = 0 vs Ha rho != 0 (conditional on all other parameters)
  pval0_rho2 <- compute_full_bootstrap(iterations    = bstrp_iter,
                                       full_result   = FALSE,
                                       PM            = pm_mark,
                                       PF            = pf_mark,
                                       PhiF          = phif_mark,
                                       PhiM          = phim_mark, 
                                       gamma         = gamma,
                                       rho           = rho,
                                       n_pop         = cjs_data$n,
                                       k             = cjs_data$k,
                                       tau           = tau,
                                       Betas         = c(1000,1000),
                                       Delta         = 1.0,
                                       PropF         = mean(cjs_data$female),
                                       imputed_pairs = TRUE,
                                       method_rho    = "likelihood",
                                       test0rho      = TRUE,
                                       test0gam      = FALSE)$results$pval
  
  
  # Get Pval for H0: Rho = 0 vs Ha rho != 0 (conditional on all other parameters)
  pval0_gamma2 <- compute_full_bootstrap(iterations      = bstrp_iter,
                                         full_result     = FALSE,
                                         PM              = pm_mark,
                                         PF              = pf_mark,
                                         PhiF            = phif_mark,
                                         PhiM            = phim_mark, 
                                         gamma           = gamma,
                                         rho             = rho,
                                         n_pop           = cjs_data$n,
                                         k               = cjs_data$k,
                                         tau             = tau,
                                         Betas           = c(1000,1000),
                                         Delta           = 1.0,
                                         PropF           = mean(cjs_data$female),
                                         imputed_pairs   = TRUE,
                                         method_rho      = "likelihood",
                                         test0rho        = FALSE,
                                         test0gam        = TRUE)$results$pval
  
  # # Run Bootstrap Estimates for Pearson Approach
  # summ_corr_fp_list <- compute_full_bootstrap(iterations    = bstrp_iter,
  #                                             full_result   = !small_out,
  #                                             PM            = pm_mark,
  #                                             PF            = pf_mark,
  #                                             PhiF          = phif_mark,
  #                                             PhiM          = phim_mark, 
  #                                             gamma         = pearson_gamma,
  #                                             rho           = pearson_rho,
  #                                             n_pop         = cjs_data$n,
  #                                             k             = cjs_data$k,
  #                                             tau           = tau,
  #                                             Betas         = c(1000,1000),
  #                                             Delta         = 1.0,
  #                                             PropF         = mean(cjs_data$female),
  #                                             imputed_pairs = TRUE,
  #                                             method_rho    = "full_pearson")
  # 
  # summ_corr_fp <- summ_corr_fp_list$results
  # 
  # # Run Bootstrap Estimates for Partial-Pearson Approach
  # summ_corr_pp_list <- compute_full_bootstrap(iterations    = bstrp_iter,
  #                                             full_result   = !small_out,
  #                                             PM            = pm_mark,
  #                                             PF            = pf_mark,
  #                                             PhiF          = phif_mark,
  #                                             PhiM          = phim_mark, 
  #                                             gamma         = pearson_partial_gamma,
  #                                             rho           = pearson_partial_rho,
  #                                             n_pop         = cjs_data$n,
  #                                             k             = cjs_data$k,
  #                                             tau           = tau,
  #                                             Betas         = c(1000,1000),
  #                                             Delta         = 1.0,
  #                                             PropF         = mean(cjs_data$female),
  #                                             imputed_pairs = TRUE,
  #                                             method_rho    = "partial_pearson")
  # 
  # summ_corr_pp <- summ_corr_pp_list$results
  
  # Return Results----------------------------------------------------------------------------------------------
  cat("Formatting output ...","\n")
  # Gather Correlation Results
  summ_corr           <- as.data.frame(rbind(summ_corr_lik))
  # summ_corr$Pval02    <- c(pval0_rho2, pval0_gamma2)#,0,pearson_pval0_gamma,0,pearson_partial_pval0_gamma)
  summ_corr$Pval02    <- c(pval0_rho2, pval0_gamma2)
  summ_corr           <- summ_corr[,c("Parameter","Est","Est_Btstrp","SE", 
                                      "2.5%","25%","50%","75%","97.5%", 
                                      "Rho_Estimator", "Pval0","Pval02","num_failures")]
  
  summ_corr           <- summ_corr %>% 
    mutate(Est = as.numeric(Est),
           Est_Btstrp = as.numeric(Est_Btstrp),
           SE = as.numeric(SE),
           `2.5%` = as.numeric(`2.5%`),
           `25%` = as.numeric(`25%`),
           `50%` = as.numeric(`50%`),
           `75%` = as.numeric(`75%`),
           `97.5%` = as.numeric(`97.5%`),
           Pval0  = as.numeric(Pval0),
           Pval02  = as.numeric(Pval02),
           num_failures = as.numeric(num_failures)) 
  
  # Compute CChat Adjustment
  summ_chat <- data.frame(Method     = c("Likelihood"),
                          CChatRho   = c(compute_proposed_chat(rho, 0)),
                          CChatGamma = c(compute_proposed_chat(0, gamma)),
                          CChat      = c(compute_proposed_chat(rho, gamma)))
  
  # Add Correlation Chat Correction to Mark-Recapture Results
  summ_cjs <- cjs_out %>% 
    mutate(CChatAdj_Likelihood      = ifelse(Parameter == "P", 
                                             compute_proposed_chat(rho, 0), 
                                             ifelse(Parameter == "Phi",
                                                    compute_proposed_chat(0, gamma), 
                                                    1)),
           UBAdj_Likelihood         = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Likelihood), alpha = 0.05)[["ub"]],
           LBAdj_Likelihood         = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Likelihood), alpha = 0.05)[["lb"]],
           Range95_Likelihood       = UBAdj_Likelihood - LBAdj_Likelihood) 
  
  # Summarize Sample Sizes
  summ_n <- data.frame(n_eff                       = n_eff,
                       n_eff_rho                   = n_eff_rho,
                       n_eff_gamma                 = n_eff_gamma)
  
  # Produce AICC comparisons
  summ_aic <- compute_aic_summ(summ_cjs = summ_cjs,
                               n_eff    = n_eff,
                               iter     = 0,
                               scenario = 0)
  
  cat("Success, returning results ...","\n")
  
  # Return Results
  if(small_out){
    # Only summaries and random.seeds
    out <- list(random_seed                 = random_seed,
                summ_n                      = summ_n,
                summ_cjs                    = summ_cjs,
                summ_corr                   = summ_corr,
                summ_chat                   = summ_chat,
                summ_aic                    = summ_aic)
  } else {
    # Summaries, data, bootstrap estimates
    out <- list(random_seed                 = random_seed,
                summ_n                      = summ_n,
                summ_cjs                    = summ_cjs,
                summ_corr                   = summ_corr,
                summ_chat                   = summ_chat,
                summ_aic                    = summ_aic,
                bstrp_lik_rep               = summ_corr_lik_list$replicates)
  }
  
  return(out)
  
}


format_to_cjs <- function(model_data){
  
  # Unpack data
  recap_f         <- model_data$recap_f
  recap_m         <- model_data$recap_m
  k               <- model_data$k
  af              <- model_data$af
  am              <- model_data$am
  first_capture_f <- model_data$first_capture_f
  first_capture_m <- model_data$first_capture_m
  nf              <- model_data$nf
  nm              <- model_data$nm
  
  # Build CJS data 
  x <- rbind(recap_f[1:nf,1:k], recap_m[1:nm,1:k])
  a <- rbind(af[1:nf,1:k],      am[1:nm,1:k])
  
  initial_entry <- c(first_capture_f[1:nf],first_capture_m[1:nm])
  
  # Store results in list
  results <- list(n             = nf + nm,
                  k             = k,
                  female        = c(rep(1, nf),rep(0, nm)), 
                  initial_entry = initial_entry,
                  x             = x,
                  a             = a)
  
  # Return Standard CJS Data
  return(results)
}