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

# Functions to RUn Experiments --------------------------------------------------------------------------------------------

# Run CJS Model using Mark
run_cjs_model_mark <- function(cjs_data,
                               title = NULL,
                               version = "B"){
  
  if(is.null(title)) title <- "mark_out_" %+% version
  
  # Define survival/recapture parameter space
  phi.grp     <- list(formula = ~sex)
  p.grp       <- list(formula = ~sex)
  
  #Process Mark Data (Extract Recapture/sex)
  # Add sex grouping for B,R,S only
  if(version == "N"){
    dat.process <- cjs_data$x %>% 
      as.data.frame() %>% 
      unite("ch",sep="") 
    
    mark.process <- process.data(data   = dat.process,
                                 model  = "CJS")
  } else {
    dat.process <- cjs_data$x %>% 
      as.data.frame() %>% 
      unite("ch",sep="") %>% 
      mutate(sex = as.factor(ifelse(cjs_data$female == 1, "F", "M"))) %>%
      arrange(sex)
    
    mark.process <- process.data(data   = dat.process,
                                 model  = "CJS",
                                 groups = "sex")
  }
  
  #Design Data List
  mark.ddl <- make.design.data(mark.process) 
  
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
  mark_out <- mark(data             = mark.process,
                   ddl              = mark.ddl,
                   model.parameters = model.parameters,
                   profile.int      = FALSE,
                   invisible        = TRUE,
                   brief            = TRUE,
                   delete           = TRUE,
                   prefix           = title,
                   output           = FALSE)
  
  mark_out <- mark_out$results
  
  #Extract values of use
  gof.stats <- data.frame("lnl"         = mark_out$lnl,
                          "npar"        = mark_out$npar,
                          "deviance"    = mark_out$deviance,
                          "deviance.df" = mark_out$deviance.df)  %>%
    rename(`-2lnl` = "lnl")
  
  mark.stats <- mark_out$real[,1:4] %>% t() %>% t() %>%
    unname() %>%
    data.frame() %>%
    rename("Est"= X1,
           "SE" = X2,
           "LB" = X3,
           "UB" = X4) %>%
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
    mutate(Deviance_chat = deviance/deviance.df, 
           Pearson_chat  = chat_list[["Pearson_chat"]],
           Fletcher_chat = chat_list[["Fletcher_chat"]])
  #Return Output
  return(mark_results)
  
}

# Compute Confidence Intervals of correlation using Fishers Method
compute_fisher_intervals <- function(r,N, alpha){
  fishers_z <- log((1+r)/(1-r))/2
  interval_width <- qnorm(1 - alpha/2, 0, 1)/sqrt(N-3)
  U <- fishers_z + interval_width
  L <- fishers_z - interval_width
  LB <- (exp(2 * L)- 1)/(exp(2*L) + 1)
  UB <- (exp(2 * U)- 1)/(exp(2*U) + 1)
  return(c(LB,UB))
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

# Summarize Bootstrap output
compute_btsrp_summary <- function(estimate, bstrp_replicates, parameteric, pearson){
  # Non-Parametric conditional estimator results
  mean_bstrp            <- mean(bstrp_replicates)
  names(mean_bstrp)     <- "Est_Btstrp"
  se_bstrp              <- sd(bstrp_replicates)
  names(se_bstrp)       <- "SE"
  quantiles_btstrp      <- quantile(bstrp_replicates, c(0.025, 0.25, 0.5, 0.75, 0.975))
  status_param          <- parameteric
  status_pearson        <- pearson
  names(status_param)   <- "Parametric"
  names(status_pearson) <- "Pearson"
  summ_btstrp           <- c(estimate, 
                             mean_bstrp, 
                             se_bstrp,
                             quantiles_btstrp,
                             status_param, 
                             status_pearson)
  
  return(summ_btstrp)
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

# Compute Likelihood Ratio Test and Quasi-Likelihood Ratio Tests (F-test)
compute_lrt_summ <- function(summ_cjs,
                             n_eff,
                             iter,
                             scenario){
  
  # Get nll values by model parameterization and cchat corrections
  ll_data <- summ_cjs %>% 
    group_by(Version) %>% 
    summarize(ll                    = first(`-2lnl`),
              df                    = first(npar),
              dev_df                = first(deviance.df),
              CChat_Likelihood      = prod(unique(CChatAdj_Likelihood)),
              CChat_Pearson         = prod(unique(CChatAdj_Pearson)),
              CChat_Partial_Pearson = prod(unique(CChatAdj_Partial_Pearson)))
  
  # Compute lrt for TEST1 
  test                       <- c("SB","RB","NB","NR","NS")
  null                       <- c("S","R","N","N","N")
  alt                        <- c("B","B","B","R","S")
  n_tests                    <- length(null) 
  lrt_stat                   <- rep(0,n_tests)
  lrt_df                     <- rep(0,n_tests)
  lrt_pval                   <- rep(0,n_tests)
  F_stat_likelihood          <- rep(0,n_tests)
  F_pval_likelihood          <- rep(0,n_tests)
  F_stat_pearson             <- rep(0,n_tests)
  F_pval_pearson             <- rep(0,n_tests)
  F_stat_partial_pearson     <- rep(0,n_tests)
  F_pval_partial_pearson     <- rep(0,n_tests)
  aic_delta                  <- rep(0,n_tests)
  aicc_delta                 <- rep(0,n_tests)
  aic_delta_likelihood       <- rep(0,n_tests)
  aicc_delta_likelihood      <- rep(0,n_tests)
  aic_delta_pearson          <- rep(0,n_tests)
  aicc_delta_pearson         <- rep(0,n_tests)
  aic_delta_partial_pearson  <- rep(0,n_tests)
  aicc_delta_partial_pearson <- rep(0,n_tests)
  
  
  for(i in 1:n_tests){
    # Grab relevant statistics for test i
    null_temp <- ll_data[ll_data$Version == null[i],]
    alt_temp  <- ll_data[ll_data$Version == alt[i],]
    lrt_delta <- (null_temp$ll - alt_temp$ll)
    # COmpute standard LRT 
    lrt_df[i]   <- ((alt_temp$df - null_temp$df))
    lrt_stat[i] <- lrt_delta/lrt_df[i]
    lrt_pval[i] <- pchisq(lrt_stat[i],df=lrt_df[i],lower.tail=FALSE)
    # Compute F-Test using CChat_Likelihood
    F_stat_likelihood[i] <- lrt_delta/(null_temp$CChat_Likelihood * lrt_df[i])  
    F_pval_likelihood[i] <- pf(F_stat_likelihood[i],df1 = lrt_df[i] ,df2 = alt_temp$df, lower.tail=FALSE)
    # Compute F-Test using CChat_Pearson
    F_stat_pearson[i] <- lrt_delta/(null_temp$CChat_Pearson * lrt_df[i])  
    F_pval_pearson[i] <- pf(F_stat_pearson[i],df1 = lrt_df[i] ,df2 = alt_temp$df,lower.tail=FALSE)
    # Compute F-Test using CChat_Partial_Pearson
    F_stat_partial_pearson[i] <- lrt_delta/(null_temp$CChat_Partial_Pearson * lrt_df[i])  
    F_pval_partial_pearson[i] <- pf(F_stat_partial_pearson[i],df1 = lrt_df[i],df2 = alt_temp$df,lower.tail=FALSE)
    # AIC Delta
    aic_delta[i] <- compute_aic_mark(ll = null_temp$ll, k = null_temp$df, n = n_eff, chat = 1, cc = 0) -
      compute_aic_mark(ll = alt_temp$ll,  k = alt_temp$df,  n = n_eff, chat = 1, cc = 0)
    
    aicc_delta[i] <- compute_aic_mark(ll = null_temp$ll, k = null_temp$df, n = n_eff, chat = 1, cc = 1) -
      compute_aic_mark(ll = alt_temp$ll,  k = alt_temp$df,  n = n_eff, chat = 1, cc = 1)
    
    # AIC Delta Likelihood
    aic_delta_likelihood[i] <- compute_aic_mark(ll = null_temp$ll, k = null_temp$df, n = n_eff, chat = null_temp$CChat_Likelihood,  cc = 0) -
      compute_aic_mark(ll = alt_temp$ll,  k = alt_temp$df,  n = n_eff, chat = null_temp$CChat_Likelihood,  cc = 0)
    
    aicc_delta_likelihood[i] <- compute_aic_mark(ll = null_temp$ll, k = null_temp$df, n = n_eff, chat = null_temp$CChat_Likelihood, cc = 1) -
      compute_aic_mark(ll = alt_temp$ll,  k =alt_temp$df,   n = n_eff, chat = null_temp$CChat_Likelihood, cc = 1)
    
    # AIC Delta Pearson
    aic_delta_pearson[i] <- compute_aic_mark(ll = null_temp$ll, k = null_temp$df, n = n_eff, chat = null_temp$CChat_Pearson,cc = 0) -
      compute_aic_mark(ll = alt_temp$ll,  k =alt_temp$df,  n = n_eff, chat = null_temp$CChat_Pearson, cc = 0)
    
    aicc_delta_pearson[i] <- compute_aic_mark(ll = null_temp$ll, k = null_temp$df, n = n_eff, chat = null_temp$CChat_Pearson, cc = 1) -
      compute_aic_mark(ll = alt_temp$ll,  k =alt_temp$df,   n = n_eff, chat = null_temp$CChat_Pearson, cc = 1)
    
    # AIC Delta Partial Pearson
    aic_delta_partial_pearson[i] <- compute_aic_mark(ll = null_temp$ll, k = null_temp$df, n = n_eff, chat = null_temp$CChat_Partial_Pearson, cc = 0) -
      compute_aic_mark(ll = alt_temp$ll,  k =alt_temp$df,   n = n_eff, chat = null_temp$CChat_Partial_Pearson, cc = 0)
    
    aicc_delta_partial_pearson[i] <- compute_aic_mark(ll = null_temp$ll, k = null_temp$df, n = n_eff, chat = null_temp$CChat_Partial_Pearson,  cc = 1) -
      compute_aic_mark(ll = alt_temp$ll,  k =alt_temp$df,   n = n_eff, chat = null_temp$CChat_Partial_Pearson,  cc = 1)
  }
  
  summ_lrt <- data.frame(test,
                         null,
                         alt,
                         lrt_stat,
                         lrt_df,
                         lrt_pval,
                         F_stat_likelihood,
                         F_pval_likelihood,
                         F_stat_pearson,
                         F_pval_pearson,
                         F_stat_partial_pearson,
                         F_pval_partial_pearson,
                         aic_delta,
                         aicc_delta,
                         aic_delta_likelihood,
                         aicc_delta_likelihood,
                         aic_delta_pearson,
                         aicc_delta_pearson,
                         aic_delta_partial_pearson,
                         aicc_delta_partial_pearson,
                         iter,
                         scenario)
  
  return(summ_lrt)
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
              dev_df                 = first(deviance.df),
              CChat_Likelihood       = prod(unique(CChatAdj_Likelihood)),
              CChat_Pearson          = prod(unique(CChatAdj_Pearson)),
              CChat_Partial_Pearson  = prod(unique(CChatAdj_Partial_Pearson))) %>% 
    ungroup() %>% 
    mutate(AIC                       = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = rep(1,4),              cc = 0),
           AICC                      = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = rep(1,4),              cc = 1),
           AIC_chat_likelihood       = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = CChat_Likelihood,      cc = 0),
           AICC_chat_likelihood      = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = CChat_Likelihood,      cc = 1),
           AIC_chat_pearson          = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = CChat_Pearson,         cc = 0),
           AICC_chat_pearson         = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = CChat_Pearson,         cc = 1),
           AIC_chat_partial_pearson  = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = CChat_Partial_Pearson, cc = 0),
           AICC_chat_partial_pearson = compute_aic_mark(ll = ll, k = pars, n = n_eff, chat = CChat_Partial_Pearson, cc = 1),
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


# Execute one replicate of Simulation Study 
execute_iteration  <- function(iter,
                               scenario,
                               PM,
                               PF,
                               PhiF,
                               PhiM,
                               gam_true,
                               rho_true,
                               n_pop,
                               k,
                               init          = NULL,
                               Betas         = list(beta0 = 1000, 
                                                    beta1 = 1000),
                               Delta         = 1,
                               PropF         = 0.5,
                               imputed_pairs = TRUE,
                               small_out     = TRUE){
  
  # Simulate Data ------------------------------------------------------------------------------------------------
  cat(paste0("Iteration#:", iter , " - Generating pair-swap mark-recapture data..."),"\n")
  # Parameter Grid 
  param_list <- list(
    n             = n_pop,              # Number of Animals
    k             = k,                  # Occasions
    prop.female   = PropF,              # Proportion of simulated individuals to be female
    delta         = rep(Delta, k),      # Probability that mating is attempted
    phi.f         = rep(PhiF, k),       # Marginal Prob of Female Survival
    phi.m         = rep(PhiM, k),       # Marginal Prob of Male Survival
    gam           = rep(gam_true, k),   # Correlation in Survival Prob of Mates
    p.f           = rep(PF, k),         # Marginal Prob of Female Recapture
    p.m           = rep(PM, k),         # Marginal Prob of Male Recapture
    rho           = rep(rho_true, k),   # Correlation in male survival rates
    betas         = Betas,              # logit pair reform params, beta1 is history coef, assume always repartner
    rand_init     = F,                  # Randomize Initial Entry (just leave as F)
    init          = init,               # Initial Entry into population for individual n
    show_unmated  = T,                  # Include unmated observations in attempt to mate step,
    imputed_pairs = imputed_pairs       # Impute between observed partnerships 
  )
  
  # Store Seed information to be able to reproduce datasets if needed
  random_seed <- .Random.seed
  
  # Generate One set of pair-swap data
  ps_data  <- sim_dat(param_list) 
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
    cjs_list[[i]] <- run_cjs_model_mark(cjs_data = cjs_data,
                                        title    = "mark_scenario_" %+% scenario %+% "_iter_" %+% iter %+% "_" %+% versions[i] %+% "_",
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
  
  # Non-Parametric Bootstrap To Estimate SE 
  cat(paste0("Iteration#:", iter ," - Non-Parametric Bootstrapping to get standard error estimates of rho (likelihood)..."),"\n")
  rho_bs_np <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                                 iter       = 1000,
                                                                 PF         = pf_mark,
                                                                 PM         = pm_mark,
                                                                 rho        = rho,
                                                                 use_block  = FALSE,
                                                                 parametric = FALSE,
                                                                 model      = "likelihood")
  
  # Non-Parametric conditional estimator results
  summ_rho_np              <- compute_btsrp_summary(rho, rho_bs_np, parameteric = 0, pearson = 0)             
  
  # Parametric Bootstrap To Estimate SE 
  cat(paste0("Iteration#:", iter ," - Semi-Parametric Bootstrapping to get standard error estimates of rho (likelihood)..."),"\n")
  rho_bs_sp <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                                 iter       = 1000,
                                                                 PF         = pf_mark,
                                                                 PM         = pm_mark,
                                                                 rho        = rho,
                                                                 use_block  = FALSE,
                                                                 parametric = TRUE,
                                                                 model      = "likelihood")
  
  # Semi-Parametric conditional estimator
  summ_rho_sp              <- compute_btsrp_summary(rho, rho_bs_sp, parameteric = 1, pearson = 0)       
  
  #-------------------------------------------------------------------------------------------------------------
  
  # 2. Pearson Full Approach -----------------------------------------------------------------------------------
  cat(paste0("Iteration#:", iter ," - Estimating recapture correlation rho using full pearson ..."),"\n")
  pearson_rho <- compute_recapture_correlation(ps_data = ps_data, 
                                               PF      = pf_mark,
                                               PM      = pm_mark,
                                               model   = "full_pearson")$rho
  names(pearson_rho) <- "Est"
  
  # Non-Parametric Bootstrap To Estimate SE 
  cat(paste0("Iteration#:", iter ," - Non-Parametric Bootstrapping to get standard error estimates of rho (pearson) ..."),"\n")
  rho_bs_np_pearson <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                                         iter       = 1000,
                                                                         PF         = pf_mark,
                                                                         PM         = pm_mark,
                                                                         rho        = pearson_rho,
                                                                         use_block  = FALSE,
                                                                         parametric = FALSE,
                                                                         model      = "full_pearson")
  
  # Non-Parametric conditional pearson estimator
  summ_rho_np_pearson   <- compute_btsrp_summary(pearson_rho, rho_bs_np_pearson, parameteric = 0, pearson = 1)    
  
  # Parametric Bootstrap To Estimate SE 
  cat(paste0("Iteration#:", iter ," - Semi-Parametric Bootstrapping to get standard error estimates of rho (pearson)..."),"\n")
  rho_bs_sp_pearson <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                                         iter       = 1000,
                                                                         PF         = pf_mark,
                                                                         PM         = pm_mark,
                                                                         rho        = pearson_rho,
                                                                         use_block  = FALSE,
                                                                         parametric = TRUE,
                                                                         model      = "full_pearson")
  # Semi-Parametric conditional pearson estimator
  summ_rho_sp_pearson   <- compute_btsrp_summary(pearson_rho, rho_bs_sp_pearson, parameteric = 1, pearson = 1)   
  
  #-------------------------------------------------------------------------------------------------------------
  
  # 2. Pearson Conditional Approach ----------------------------------------------------------------------------
  cat(paste0("Iteration#:", iter ," - Estimating recapture correlation rho using pseudo-pearson..."),"\n")
  pearson_partial_rho <- compute_recapture_correlation(ps_data = ps_data, 
                                                       PF      = pf_mark,
                                                       PM      = pm_mark,
                                                       model   = "partial_pearson")$rho
  names(pearson_partial_rho) <- "Est"
  
  # Non-Parametric Bootstrap To Estimate SE 
  cat(paste0("Iteration#:", iter ," - Non-Parametric Bootstrapping to get standard error estimates of rho (Psuedo-Pearson)..."),"\n")
  rho_bs_np_pearson_partial <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                                                 iter       = 1000,
                                                                                 PF         = pf_mark,
                                                                                 PM         = pm_mark,
                                                                                 rho        = pearson_partial_rho,
                                                                                 use_block  = FALSE,
                                                                                 parametric = FALSE,
                                                                                 model      = "partial_pearson")
  
  # Non-Parametric conditional pearson estimator
  summ_rho_np_pearson_partial   <- compute_btsrp_summary(pearson_partial_rho, rho_bs_np_pearson_partial, parameteric = 0, pearson = 2)   
  
  # Parametric Bootstrap To Estimate SE 
  cat(paste0("Iteration#:", iter , " - Semi-Parametric Bootstrapping to get standard error estimates of rho (Psuedo-Pearson)..."),"\n")
  rho_bs_sp_pearson_partial <- compute_bootstrap_estimates_recapture_correlation(ps_data    = ps_data,
                                                                                 iter       = 1000,
                                                                                 PF         = pf_mark,
                                                                                 PM         = pm_mark,
                                                                                 rho        = pearson_partial_rho,
                                                                                 use_block  = FALSE,
                                                                                 parametric = TRUE,
                                                                                 model      = "partial_pearson")
  # Semi-Parametric conditional pearson estimator
  summ_rho_sp_pearson_partial   <- compute_btsrp_summary(pearson_partial_rho, rho_bs_sp_pearson_partial, parameteric = 1, pearson = 2)   
  #-------------------------------------------------------------------------------------------------------------
  
  # Compute Survival Correlation Estimate-----------------------------------------------------------------------
  
  # Conditional on Rho from Likelihood Approach ----------------------------------------------------------------
  cat(paste0("Iteration#:", iter ," - Estimating survival correlation gamma|rho-likelihood..."),"\n")
  
  # If estimate of rho fails (no valid observations) pass dummy values of 10
  if(rho == 10){
    gamma <- 10
    gamma_bs_np <- rep(10, 1000)
    gamma_bs_sp <- rep(10, 1000)
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
    
    # Non-Parametric Bootstrap To Estimate SE 
    cat(paste0("Iteration#:", iter ," - Non-Parametric Bootstrapping to get standard error estimates of gamma|rho-likelihood..."),"\n")
    gamma_bs_np <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                    iter                  = 1000,
                                                                    rho                   = rho,
                                                                    PF                    = pf_mark,
                                                                    PM                    = pm_mark,
                                                                    gamma                 = gamma,
                                                                    PhiF                  = phif_mark,
                                                                    PhiM                  = phim_mark,
                                                                    use_block             = FALSE,
                                                                    parametric            = FALSE)
    
    
    # Semi-Parametric Bootstrap To Estimate SE 
    cat(paste0("Iteration#:", iter ," - Semi-Parametric Bootstrapping to get standard error estimates of gamma|rho-likelihood..."),"\n")
    gamma_bs_sp <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                    iter                  = 1000,
                                                                    rho                   = rho,
                                                                    PF                    = pf_mark,
                                                                    PM                    = pm_mark,
                                                                    gamma                 = gamma,
                                                                    PhiF                  = phif_mark,
                                                                    PhiM                  = phim_mark,
                                                                    use_block             = FALSE,
                                                                    parametric            = TRUE)
  }
  
  # Collect Results
  names(gamma) <- "Est"
  summ_gamma_np   <- compute_btsrp_summary(gamma, gamma_bs_np, parameteric = 0, pearson = 0)   
  summ_gamma_sp   <- compute_btsrp_summary(gamma, gamma_bs_sp, parameteric = 1, pearson = 0)   
  #-------------------------------------------------------------------------------------------------------------
  
  # Conditional on Rho from Likelihood Approach ----------------------------------------------------------------
  cat(paste0("Iteration#:", iter ," - Estimating survival correlation gamma|rho-pearson..."),"\n")
  
  # If estimate of rho fails (no valid observations) pass dummy values of 10
  if(is.na(pearson_rho)|pearson_rho == 10){
    pearson_gamma <- 10
    pearson_gamma_bs_np <- rep(10, 1000)
    pearson_gamma_bs_sp <- rep(10, 1000)
  } else {
    
    # Estimate Gamma from observed data
    pearson_gamma <- compute_survival_correlation(ps_data = ps_data,
                                                  PFM     = compute_jbin_cjs(prob.f = pf_mark,
                                                                             prob.m = pm_mark,
                                                                             corr   = pearson_rho)$prob.mf,
                                                  PhiF    = phif_mark,
                                                  PhiM    = phim_mark)$gamma
    
    # Non-Parametric Bootstrap To Estimate SE 
    cat(paste0("Iteration#:", iter , " - Non-Parametric Bootstrapping to get standard error estimates of gamma|rho-pearson..."),"\n")
    pearson_gamma_bs_np <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                            iter                  = 1000,
                                                                            rho                   = pearson_rho,
                                                                            PF                    = pf_mark,
                                                                            PM                    = pm_mark,
                                                                            gamma                 = pearson_gamma,
                                                                            PhiF                  = phif_mark,
                                                                            PhiM                  = phim_mark,
                                                                            use_block             = FALSE,
                                                                            parametric            = FALSE)
    
    
    # Semi-Parametric Bootstrap To Estimate SE 
    cat(paste0("Iteration#:", iter , " - Semi-Parametric Bootstrapping to get standard error estimates of gamma|rho-pearson..."),"\n")
    pearson_gamma_bs_sp <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                            iter                  = 1000,
                                                                            rho                   = pearson_rho,
                                                                            PF                    = pf_mark,
                                                                            PM                    = pm_mark,
                                                                            gamma                 = pearson_gamma,
                                                                            PhiF                  = phif_mark,
                                                                            PhiM                  = phim_mark,
                                                                            use_block             = FALSE,
                                                                            parametric            = TRUE)
  }
  
  # Collect Results
  names(pearson_gamma) <- "Est"
  summ_pearson_gamma_np   <- compute_btsrp_summary(pearson_gamma, pearson_gamma_bs_np, parameteric = 0, pearson = 1)   
  summ_pearson_gamma_sp   <- compute_btsrp_summary(pearson_gamma, pearson_gamma_bs_sp, parameteric = 1, pearson = 1) 
  #-------------------------------------------------------------------------------------------------------------
  
  # Conditional on Rho from Likelihood Approach ----------------------------------------------------------------
  cat(paste0("Iteration#:", iter , " - Estimating survival correlation gamma|rho-psuedo-pearson..."),"\n")
  
  # If estimate of rho fails (no valid observations) pass dummy values of 10
  if(is.na(pearson_partial_rho)|pearson_partial_rho == 10){
    pearson_partial_gamma <- 10
    pearson_partial_gamma_bs_np <- rep(10, 1000)
    pearson_partial_gamma_bs_sp <- rep(10, 1000)
  } else {
    
    # Estimate Gamma from observed data
    pearson_partial_gamma <- compute_survival_correlation(ps_data = ps_data,
                                                          PFM     = compute_jbin_cjs(prob.f = pf_mark,
                                                                                     prob.m = pm_mark,
                                                                                     corr   = pearson_partial_rho)$prob.mf,
                                                          PhiF    = phif_mark,
                                                          PhiM    = phim_mark)$gamma
    
    # Non-Parametric Bootstrap To Estimate SE 
    cat(paste0("Iteration#:", iter , " - Non-Parametric Bootstrapping to get standard error estimates of gamma|rho-psuedo-pearson..."),"\n")
    pearson_partial_gamma_bs_np <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                                    iter                  = 1000,
                                                                                    rho                   = pearson_partial_rho,
                                                                                    PF                    = pf_mark,
                                                                                    PM                    = pm_mark,
                                                                                    gamma                 = pearson_partial_gamma,
                                                                                    PhiF                  = phif_mark,
                                                                                    PhiM                  = phim_mark,
                                                                                    use_block             = FALSE,
                                                                                    parametric            = FALSE)
    
    
    # Semi-Parametric Bootstrap To Estimate SE 
    cat(paste0("Iteration#:", iter , " - Semi-Parametric Bootstrapping to get standard error estimates of gamma|rho-psuedo-pearson..."),"\n")
    pearson_partial_gamma_bs_sp <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                                                    iter                  = 1000,
                                                                                    rho                   = pearson_partial_rho,
                                                                                    PF                    = pf_mark,
                                                                                    PM                    = pm_mark,
                                                                                    gamma                 = pearson_partial_gamma,
                                                                                    PhiF                  = phif_mark,
                                                                                    PhiM                  = phim_mark,
                                                                                    use_block             = FALSE,
                                                                                    parametric            = TRUE)
  }
  
  # Collect Results
  names(pearson_partial_gamma) <- "Est"
  summ_pearson_partial_gamma_np   <- compute_btsrp_summary(pearson_partial_gamma, pearson_partial_gamma_bs_np, parameteric = 0, pearson = 2)   
  summ_pearson_partial_gamma_sp   <- compute_btsrp_summary(pearson_partial_gamma, pearson_partial_gamma_bs_sp, parameteric = 1, pearson = 2) 
  #-------------------------------------------------------------------------------------------------------------
  
  # Return Results----------------------------------------------------------------------------------------------
  
  cat("Formatting output ...","\n")
  
  # Gather Correlation Results
  summ_corr           <- as.data.frame(rbind(summ_rho_np,
                                             summ_rho_sp,
                                             summ_rho_np_pearson,
                                             summ_rho_sp_pearson, 
                                             summ_rho_np_pearson_partial,
                                             summ_rho_sp_pearson_partial,
                                             summ_gamma_np, 
                                             summ_gamma_sp,
                                             summ_pearson_gamma_np,
                                             summ_pearson_gamma_sp,
                                             summ_pearson_partial_gamma_np,
                                             summ_pearson_partial_gamma_sp))
  
  param_true          <- c(rep(rho_true, 6), rep(gam_true,6))
  summ_corr$Parameter <- c(rep("rho", 6), rep("gamma",6))
  summ_corr           <- summ_corr[,c("Parameter","Est","Est_Btstrp","SE", "2.5%","25%","50%","75%","97.5%", "Parametric", "Pearson")]
  summ_corr           <- summ_corr %>% 
    left_join(true_param_df, by = "Parameter") %>% 
    mutate(Bias                = Truth - Est,
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
  
  rownames(summ_corr) <- NULL
  
  # Compute CChat Adjustment
  summ_chat <- data.frame(Method     = c("Likelihood", "Pearson", "Psuedo-Pearson", "Truth"),
                          CChatRho   = c(compute_proposed_chat(rho, 0), 
                                         compute_proposed_chat(pearson_rho, 0), 
                                         compute_proposed_chat(pearson_partial_rho, 0),
                                         compute_proposed_chat(rho_true, 0)),
                          CChatGamma = c(compute_proposed_chat(0, gamma), 
                                         compute_proposed_chat(0, pearson_gamma), 
                                         compute_proposed_chat(0, pearson_partial_gamma),
                                         compute_proposed_chat(0, gam_true)),
                          CChat      = c(compute_proposed_chat(rho, gamma), 
                                         compute_proposed_chat(pearson_rho, pearson_gamma), 
                                         compute_proposed_chat(pearson_partial_rho, pearson_partial_gamma),
                                         compute_proposed_chat(rho_true, gam_true)),
                          iter       = rep(iter, 4),
                          scenario   = rep(scenario, 4))
  
  # Add Correlation Chat Correction to Mark-Recapture Results
  summ_cjs <- cjs_out %>% 
    mutate(CChatAdj_Likelihood      = ifelse(Parameter == "P", 
                                             compute_proposed_chat(rho, 0), 
                                             ifelse(Parameter == "Phi",
                                                    compute_proposed_chat(0, gamma), 
                                                    1)),
           CChatAdj_Pearson         = ifelse(Parameter == "P", 
                                             compute_proposed_chat(pearson_rho, 0), 
                                             ifelse(Parameter == "Phi",
                                                    compute_proposed_chat(0, pearson_gamma), 
                                                    1)),
           CChatAdj_Partial_Pearson = ifelse(Parameter == "P", 
                                             compute_proposed_chat(pearson_partial_rho, 0), 
                                             ifelse(Parameter == "Phi",
                                                    compute_proposed_chat(0, pearson_partial_gamma), 
                                                    1)),
           CChatAdj_Truth           = ifelse(Parameter == "P", 
                                             compute_proposed_chat(rho_true, 0), 
                                             ifelse(Parameter == "Phi",
                                                    compute_proposed_chat(0, gam_true), 
                                                    1)),
           UBAdj_Likelihood         = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Likelihood), alpha = 0.05)[["ub"]],
           LBAdj_Likelihood         = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Likelihood), alpha = 0.05)[["lb"]],
           UBAdj_Pearson            = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Pearson), alpha = 0.05)[["ub"]],
           LBAdj_Pearson            = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Pearson), alpha = 0.05)[["lb"]],
           UBAdj_Partial_Pearson    = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Partial_Pearson), alpha = 0.05)[["ub"]],
           LBAdj_Partial_Pearson    = compute_mark_ci(prob =  Est, se = SE * sqrt(CChatAdj_Partial_Pearson), alpha = 0.05)[["lb"]],
           In95_Likelihood          = 1*(Truth <= UBAdj_Likelihood & Truth >= LBAdj_Likelihood),
           In95_Pearson             = 1*(Truth <= UBAdj_Pearson & Truth >= LBAdj_Pearson),
           In95_Partial_Pearson     = 1*(Truth <= UBAdj_Partial_Pearson & Truth >= LBAdj_Partial_Pearson),
           Range95_Likelihood       = UBAdj_Likelihood - LBAdj_Likelihood,
           Range95_Pearson          = UBAdj_Pearson - LBAdj_Pearson,
           Range95_Partial_Pearson  = UBAdj_Partial_Pearson - LBAdj_Partial_Pearson
    ) 
  
  # Produce Likelihood Ratio Tests
  summ_lrt <- compute_lrt_summ(summ_cjs = summ_cjs,
                               n_eff    = n_eff,
                               iter     = iter,
                               scenario = scenario)
  
  
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
                    summ_lrt                    = summ_lrt,
                    summ_aic                    = summ_aic)
  } else {
    results <- list(random_seed                 = random_seed,
                    ps_data                     = ps_data,
                    cjs_data                    = cjs_data,
                    summ_n                      = summ_n,
                    summ_cjs                    = summ_cjs,
                    summ_corr                   = summ_corr,
                    summ_chat                   = summ_chat,
                    summ_lrt                    = summ_lrt,
                    summ_aic                    = summ_aic,
                    rho_bs_np                   = unname(rho_bs_np),
                    rho_bs_sp                   = unname(rho_bs_sp),
                    rho_bs_np_pearson           = unname(rho_bs_np_pearson),
                    rho_bs_sp_pearson           = unname(rho_bs_sp_pearson),
                    rho_bs_np_pearson_partial   = unname(rho_bs_np_pearson_partial),
                    rho_bs_sp_pearson_partial   = unname(rho_bs_sp_pearson_partial),
                    gamma_bs_np                 = unname(gamma_bs_np),
                    gamma_bs_sp                 = unname(gamma_bs_sp),
                    pearson_gamma_bs_np         = unname(pearson_gamma_bs_np),
                    pearson_gamma_bs_sp         = unname(pearson_gamma_bs_sp),
                    pearson_partial_gamma_bs_np = unname(pearson_partial_gamma_bs_np),
                    pearson_partial_gamma_bs_sp = unname(pearson_partial_gamma_bs_sp))
  }
  
  
  return(results)
}

# Execute entire simulation study for correlation estimators
execute_simulation <- function(niter,
                               scenario,
                               PM,
                               PF,
                               PhiF,
                               PhiM,
                               gam_true,
                               rho_true,
                               n_pop,
                               k,
                               Betas         = list(beta0 = 1000, 
                                                    beta1 = 1000),
                               Delta         = 1,
                               PropF         = 0.5,
                               init          = NULL,
                               imputed_pairs = TRUE,
                               small_out     = TRUE){
  # Run niter replicates 
  results_list <- lapply(1:niter, function(iter) execute_iteration(iter          = iter,
                                                                   scenario      = scenario,
                                                                   PM            = PM,
                                                                   PF            = PF,
                                                                   PhiF          = PhiF,
                                                                   PhiM          = PhiM,
                                                                   gam_true      = gam_true,
                                                                   rho_true      = rho_true,
                                                                   n_pop         = n_pop,
                                                                   k             = k,
                                                                   Betas         = Betas,
                                                                   Delta         = Delta,
                                                                   PropF         = PropF,
                                                                   init          = init,
                                                                   imputed_pairs = imputed_pairs,
                                                                   small_out     = small_out))
  
  # Construct Summaries of Probs and Correlations
  summary_corr <- do.call(rbind, lapply(1:niter, function(iter) results_list[[iter]]$summ_corr))
  summary_cjs  <- do.call(rbind, lapply(1:niter, function(iter) results_list[[iter]]$summ_cjs))
  summ_chat    <- do.call(rbind, lapply(1:niter, function(iter) results_list[[iter]]$summ_chat))
  summ_lrt     <- do.call(rbind, lapply(1:niter, function(iter) results_list[[iter]]$summ_lrt))
  summ_aic     <- do.call(rbind, lapply(1:niter, function(iter) results_list[[iter]]$summ_aic))
  summ_n       <- do.call(rbind, lapply(1:niter, function(iter) results_list[[iter]]$summ_n))
  random_seeds <- lapply(1:niter, function(iter) results_list[[iter]]$random_seed)
  
  # Return Results
  if(small_out){
    # Only summaries and random.seeds
    out <- list(summary_corr = summary_corr,
                summary_cjs  = summary_cjs,
                summ_chat    = summ_chat,
                summ_lrt     = summ_lrt,
                summ_aic     = summ_aic,
                summ_n       = summ_n,
                random_seeds = random_seeds)
  } else {
    # Summaries, data, bootstrap estimates
    out <- list(results_list = results_list,
                summary_corr = summary_corr,
                summary_cjs  = summary_cjs,
                summ_chat    = summ_chat,
                summ_lrt     = summ_lrt,
                summ_aic     = summ_aic,
                summ_n       = summ_n)
  }
  
  
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
  rho   <- sort(unique(c(0, seq(-0.1, 0.9, by = 0.25))))
  gamma <- sort(unique(c(0, seq(-0.1, 0.9, by = 0.25))))
  
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
  
  rho2   <- sort(unique(c(0, seq(rl + 0.01, ru - 0.01, by = 0.25))))
  gamma2 <- sort(unique(c(0, seq(gl + 0.01, gu - 0.01, by = 0.25))))
  
  scenario_grid_alternative <- expand.grid(n_obs         = 250,
                                           k             = 25,
                                           rho_true      = rho2,
                                           gam_true      = gamma2, 
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
  
  rho3   <- sort(unique(c(0,rl + 0.01, ru - 0.01, seq(rl + 0.01, ru - 0.01, by = 0.25))))
  gamma3 <- sort(unique(c(0, gl + 0.01, gu - 0.01,seq(gl + 0.01, gu - 0.01, by = 0.25))))
  
  scenario_grid_hduck1 <- expand.grid(n_obs         = 250,
                                      k             = 25,
                                      rho_true      = rho3,
                                      gam_true      = gamma3, 
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
  
  rho4   <- sort(unique(c(0,rl + 0.01, ru - 0.01, seq(rl + 0.01, ru - 0.01, by = 0.25))))
  gamma4 <- sort(unique(c(0,gl + 0.01, gu - 0.01,seq(gl + 0.01, gu - 0.01, by = 0.25))))
  
  scenario_grid_hduck2 <- expand.grid(n_obs         = 250,
                                      k             = 25,
                                      rho_true      = rho4,
                                      gam_true      = gamma4, 
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
  
  rho5   <- sort(unique(c(0,rl + 0.01, ru - 0.01, seq(rl + 0.01, ru - 0.01, by = 0.25))))
  gamma5 <- sort(unique(c(0, gl + 0.01, gu - 0.01,seq(gl + 0.01, gu - 0.01, by = 0.25))))
  
  scenario_grid_hduck3 <- expand.grid(n_obs         = 250,
                                      k             = 25,
                                      rho_true      = rho5,
                                      gam_true      = gamma5, 
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
  
  rho6   <- sort(unique(c(0,rl + 0.01, ru - 0.01, seq(rl + 0.01, ru - 0.01, by = 0.25))))
  gamma6 <- sort(unique(c(0, gl + 0.01, gu - 0.01,seq(gl + 0.01, gu - 0.01, by = 0.25))))
  
  scenario_grid_hduck4 <- expand.grid(n_obs         = 250,
                                      k             = 25,
                                      rho_true      = rho6,
                                      gam_true      = gamma6, 
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
                                      rho_true      = rho,
                                      gam_true      = gamma, 
                                      PhiF          = 0.8,
                                      PhiM          = 0.8,
                                      PF            = 0.75,
                                      PM            = 0.75,
                                      Beta0         = 0.25,
                                      Beta1         = 0,
                                      Delta         = 1,
                                      imputed_pairs = FALSE)
  
  
  
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