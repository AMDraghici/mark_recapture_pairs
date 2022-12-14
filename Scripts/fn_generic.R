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


# Numerically Stable Version
softmax <- function(par){
  n.par <- length(par)
  par1  <- sort(par, decreasing = TRUE)
  Lk    <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
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
                               title = NULL){
  
  
  if(is.null(title)) title <- "mark_out"
  
  #Choose Appropriate CJS Model Settings
  phi.grp     <- list(formula = ~sex)
  p.grp       <- list(formula = ~sex)
  phi.name    <- c("PhiM","PhiF")
  p.name      <- c("PM","PF")
  param.names <- c(phi.name,p.name)
  
  #Process Mark Data (Extract Recapture/sex)
  dat.process <- cjs_data$x %>% 
    as.data.frame() %>% 
    unite("ch",sep="") %>% 
    mutate(sex = as.factor(cjs_data$female)) %>%
    arrange(sex) 
  
  mark.process <- process.data(data   = dat.process,
                               model  = "CJS",
                               groups = "sex")
  
  #Design Data List
  mark.ddl <- make.design.data(mark.process) 
  
  #Generate Estimates
  mark_out <- mark(mark.process,
                   mark.ddl,
                   model.parameters =list(Phi = phi.grp,
                                          p   = p.grp),
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
    mutate(Parameter = param.names)
  
  mark_results <- cbind(mark.stats,gof.stats)
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
                               init = NULL){
  
  # Simulate Data ------------------------------------------------------------------------------------------------
  cat(paste0("Iteration#:", iter , " - Generating pair-swap mark-recapture data..."),"\n")
  # Parameter Grid 
  param_list <- list(
    n            = n_pop,              # Number of Animals
    k            = k,                  # Occasions
    prop.female  = 0.5,                # Proportion of simulated individuals to be female
    delta        = rep(1, k),          # Probability that mating is attempted
    phi.f        = rep(PhiF, k),       # Marginal Prob of Female Survival
    phi.m        = rep(PhiM, k),       # Marginal Prob of Male Survival
    gam          = rep(gam_true, k),   # Correlation in Survival Prob of Mates
    p.f          = rep(PF, k),         # Marginal Prob of Female Recapture
    p.m          = rep(PM, k),         # Marginal Prob of Male Recapture
    rho          = rep(rho_true, k),   # Correlation in male survival rates
    betas        = list(beta0 = 1000, 
                        beta1 = 1000), # logit pair reform params, beta1 is history coef, assume always repartner
    rand_init    = F,                  # Randomize Initial Entry (just leave as F)
    init         = init,               # Initial Entry into population for individual n
    show_unmated = T                   # Include unmated observations in attempt to mate step 
  )
  
  # Store Seed information to be able to reproduce datasets if needed
  random_seed <- .Random.seed
  
  # Generate One set of pair-swap data
  ps_data <- sim_dat(param_list) 
  cjs_data <- format_to_cjs(ps_data)
  
  # Gather True Param DF 
  ru_true            <- compute_jbin_param_cjs(PF,PM)$cor_lower_bound
  rl_true            <- compute_jbin_param_cjs(PF,PM)$cor_upper_bound
  gu_true            <- compute_jbin_param_cjs(PhiF,PhiM)$cor_upper_bound
  gl_true            <- compute_jbin_param_cjs(PhiF,PhiM)$cor_upper_bound
  true_params        <- c(PF,PM,PhiF, PhiM, gam_true, gl_true, gu_true, rho_true, rl_true, ru_true)
  param_names        <- c("PF","PM", "PhiF", "PhiM", "gamma", "gl", "gu", "rho", "rl", "ru")
  true_param_df      <- data.frame(Truth = true_params, Parameter = param_names)
  
  # Run CJS Model MARK -----------------------------------------------------------------------------------------
  cat(paste0("Iteration#:", iter , " - Computing standard CJS estimates with program MARK..."),"\n")
  cjs_out <- run_cjs_model_mark(cjs_data = cjs_data,
                                title    = "mark_" %+% iter %+% "_" %+% scenario) %>% 
    left_join(true_param_df,
              by = "Parameter") %>% 
    mutate(Bias     = Truth - Est,
           In95     = 1*(Truth <= UB & Truth >= LB),
           Range95  = UB - LB,
           iter     = iter,
           scenario = scenario) 
  
  
  # Get Predicted Probs for Correlation Estimators 
  phim_mark <- cjs_out %>% filter(Parameter == "PhiM") %>% pull(Est)
  phif_mark <- cjs_out %>% filter(Parameter == "PhiF") %>% pull(Est)
  pm_mark   <- cjs_out %>% filter(Parameter == "PM") %>% pull(Est)
  pf_mark   <- cjs_out %>% filter(Parameter == "PF") %>% pull(Est)
  #-------------------------------------------------------------------------------------------------------------
  
  # Compute Recapture Correlation Estimate----------------------------------------------------------------------
  
  # 1. Likelihood Approach -------------------------------------------------------------------------------------
  cat(paste0("Iteration#:", iter ," - Estimating recapture correlation, rho, using likelihood approach..."),"\n")
  rho <- compute_recapture_correlation(ps_data = ps_data, 
                                       PF      = pf_mark,
                                       PM      = pm_mark,
                                       model   = "likelihood")
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
  mean_bstrp_rho_np        <- mean(rho_bs_np)
  names(mean_bstrp_rho_np) <- "Est_Btstrp"
  se_bstrp_rho_np          <- sd(rho_bs_np)
  names(se_bstrp_rho_np)   <- "SE"
  quantiles_rho_np         <- quantile(rho_bs_np, c(0.025, 0.25, 0.5, 0.75, 0.975))
  status_np                <- 0
  pearson_status_np        <- 0
  names(status_np)         <- "Parametric"
  names(pearson_status_np) <- "Pearson"
  summ_rho_np              <- c(rho, 
                                mean_bstrp_rho_np, 
                                se_bstrp_rho_np,
                                quantiles_rho_np,
                                status_np, 
                                pearson_status_np)
  
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
  mean_bstrp_rho_sp        <- mean(rho_bs_sp)
  names(mean_bstrp_rho_sp) <- "Est_Btstrp"
  se_bstrp_rho_sp          <- sd(rho_bs_sp)
  names(se_bstrp_rho_sp)   <- "SE"
  quantiles_rho_sp         <- quantile(rho_bs_sp, c(0.025, 0.25, 0.5, 0.75, 0.975))
  status_sp                <- 1
  pearson_status_sp        <- 0
  names(status_sp)         <- "Parametric"
  names(pearson_status_sp) <- "Pearson"
  summ_rho_sp              <- c(rho, 
                                mean_bstrp_rho_sp,
                                se_bstrp_rho_sp, 
                                quantiles_rho_sp,
                                status_sp,
                                pearson_status_sp)
  
  #-------------------------------------------------------------------------------------------------------------
  
  # 2. Pearson Full Approach -----------------------------------------------------------------------------------
  cat(paste0("Iteration#:", iter ," - Estimating recapture correlation rho using full pearson ..."),"\n")
  pearson_rho <- compute_recapture_correlation(ps_data = ps_data, 
                                               PF      = pf_mark,
                                               PM      = pm_mark,
                                               model   = "full_pearson")
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
  mean_bstrp_rho_np_pearson        <- mean(rho_bs_np_pearson)
  names(mean_bstrp_rho_np_pearson) <- "Est_Btstrp"
  se_bstrp_rho_np_pearson          <- sd(rho_bs_np_pearson)
  names(se_bstrp_rho_np_pearson)   <- "SE"
  quantiles_rho_np_pearson         <- quantile(rho_bs_np_pearson, c(0.025, 0.25, 0.5, 0.75, 0.975))
  status_np_pearson                <- 0
  pearson_status_pearson_np        <- 1
  names(status_np_pearson)         <- "Parametric"
  names(pearson_status_pearson_np) <- "Pearson"
  summ_rho_np_pearson              <- c(pearson_rho, 
                                        mean_bstrp_rho_np_pearson, 
                                        se_bstrp_rho_np_pearson, 
                                        quantiles_rho_np_pearson,
                                        status_np_pearson, 
                                        pearson_status_pearson_np)
  
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
  mean_bstrp_rho_sp_pearson         <- mean(rho_bs_sp_pearson)
  names(mean_bstrp_rho_sp_pearson)  <- "Est_Btstrp"
  se_bstrp_rho_sp_pearson           <- sd(rho_bs_sp_pearson)
  names(se_bstrp_rho_sp_pearson)    <- "SE"
  quantiles_rho_sp_pearson         <- quantile(rho_bs_sp_pearson, c(0.025, 0.25, 0.5, 0.75, 0.975))
  status_sp_pearson                <- 1
  pearson_status_pearson_sp        <- 1
  names(status_sp_pearson)         <- "Parametric"
  names(pearson_status_pearson_sp) <- "Pearson"
  summ_rho_sp_pearson              <- c(pearson_rho, 
                                        mean_bstrp_rho_sp_pearson,
                                        se_bstrp_rho_sp_pearson, 
                                        quantiles_rho_sp_pearson, 
                                        status_sp_pearson, 
                                        pearson_status_pearson_sp)
  
  #-------------------------------------------------------------------------------------------------------------
  
  # 2. Pearson Conditional Approach ----------------------------------------------------------------------------
  cat(paste0("Iteration#:", iter ," - Estimating recapture correlation rho using pseudo-pearson..."),"\n")
  pearson_partial_rho <- compute_recapture_correlation(ps_data = ps_data, 
                                                       PF      = pf_mark,
                                                       PM      = pm_mark,
                                                       model   = "partial_pearson")
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
  mean_bstrp_rho_np_pearson_partial        <- mean(rho_bs_np_pearson_partial)
  names(mean_bstrp_rho_np_pearson_partial) <- "Est_Btstrp"
  se_bstrp_rho_np_pearson_partial          <- sd(rho_bs_np_pearson_partial)
  names(se_bstrp_rho_np_pearson_partial)   <- "SE"
  quantiles_rho_np_pearson_partial         <- quantile(rho_bs_np_pearson_partial, c(0.025, 0.25, 0.5, 0.75, 0.975))
  status_np_pearson_partial                <- 0
  pearson_status_pearson_partial_np        <- 2
  names(status_np_pearson_partial)         <- "Parametric"
  names(pearson_status_pearson_partial_np) <- "Pearson"
  summ_rho_np_pearson_partial              <- c(pearson_partial_rho, 
                                                mean_bstrp_rho_np_pearson_partial,
                                                se_bstrp_rho_np_pearson_partial,
                                                quantiles_rho_np_pearson_partial, 
                                                status_np_pearson_partial, 
                                                pearson_status_pearson_partial_np)
  
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
  mean_bstrp_rho_sp_pearson_partial         <- mean(rho_bs_sp_pearson_partial)
  names(mean_bstrp_rho_sp_pearson_partial)  <- "Est_Btstrp"
  se_bstrp_rho_sp_pearson_partial           <- sd(rho_bs_sp_pearson_partial)
  names(se_bstrp_rho_sp_pearson_partial)    <- "SE"
  quantiles_rho_sp_pearson_partial         <- quantile(rho_bs_sp_pearson_partial, c(0.025, 0.25, 0.5, 0.75, 0.975))
  status_sp_pearson_partial                <- 1
  pearson_status_pearson_partial_sp        <- 2
  names(status_sp_pearson_partial)         <- "Parametric"
  names(status_sp_pearson_partial)         <- "Pearson"
  summ_rho_sp_pearson_partial              <- c(pearson_partial_rho, 
                                                mean_bstrp_rho_sp_pearson_partial, 
                                                se_bstrp_rho_sp_pearson_partial, 
                                                quantiles_rho_sp_pearson_partial, 
                                                status_sp_pearson_partial, 
                                                pearson_status_pearson_partial_sp)
  
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
    gamma <- compute_survival_correlation(ps_data = ps_data,
                                          PFM     = compute_jbin_cjs(prob.f = pf_mark,
                                                                     prob.m = pm_mark,
                                                                     corr   = rho)$prob.mf,
                                          PhiF    = phif_mark,
                                          PhiM    = phim_mark)
    
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
  names(gamma)                   <- "Est"
  mean_bstrp_gamma_np            <- mean(gamma_bs_np)
  mean_bstrp_gamma_sp            <- mean(gamma_bs_sp)
  names(mean_bstrp_gamma_np)     <- "Est_Btstrp"
  names(mean_bstrp_gamma_sp)     <- "Est_Btstrp"
  se_bstrp_gamma_np              <- sd(gamma_bs_np)
  se_bstrp_gamma_sp              <- sd(gamma_bs_sp)
  names(se_bstrp_gamma_np)       <- "SE"
  names(se_bstrp_gamma_sp)       <- "SE"
  quantiles_gamma_np             <- quantile(gamma_bs_np, c(0.025, 0.25, 0.5, 0.75, 0.975))
  quantiles_gamma_sp             <- quantile(gamma_bs_sp, c(0.025, 0.25, 0.5, 0.75, 0.975))
  status_gamma_np                <- 0
  status_gamma_sp                <- 1
  names(status_gamma_np)         <- "Parametric"
  names(status_gamma_sp)         <- "Parametric"
  pearson_status_gamma_np        <- 0
  pearson_status_gamma_sp        <- 0
  names(pearson_status_gamma_np) <- "Pearson"
  names(pearson_status_gamma_sp) <- "Pearson"
  summ_gamma_np                  <- c(gamma, 
                                      mean_bstrp_gamma_np,
                                      se_bstrp_gamma_np,
                                      quantiles_gamma_np,
                                      status_gamma_np,
                                      pearson_status_gamma_np)
  summ_gamma_sp                  <- c(gamma, 
                                      mean_bstrp_gamma_sp,
                                      se_bstrp_gamma_sp,
                                      quantiles_gamma_sp,
                                      status_gamma_sp,
                                      pearson_status_gamma_sp)
  
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
                                                  PhiM    = phim_mark)
    
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
  names(pearson_gamma)                       <- "Est"
  mean_bstrp_pearson_gamma_np                <- mean(pearson_gamma_bs_np)
  mean_bstrp_pearson_gamma_sp                <- mean(pearson_gamma_bs_sp)
  names(mean_bstrp_pearson_gamma_np)         <- "Est_Btstrp"
  names(mean_bstrp_pearson_gamma_sp)         <- "Est_Btstrp"
  se_bstrp_pearson_gamma_np                  <- sd(pearson_gamma_bs_np)
  se_bstrp_pearson_gamma_sp                  <- sd(pearson_gamma_bs_sp)
  names(se_bstrp_pearson_gamma_np)           <- "SE"
  names(se_bstrp_pearson_gamma_sp)           <- "SE"
  quantiles_pearson_gamma_np                 <- quantile(pearson_gamma_bs_np, c(0.025, 0.25, 0.5, 0.75, 0.975))
  quantiles_pearson_gamma_sp                 <- quantile(pearson_gamma_bs_sp, c(0.025, 0.25, 0.5, 0.75, 0.975))
  status_pearson_gamma_np                    <- 0
  status_pearson_gamma_sp                    <- 1
  names(status_pearson_gamma_np)             <- "Parametric"
  names(status_pearson_gamma_sp)             <- "Parametric"
  pearson_status_pearson_gamma_np            <- 1
  pearson_status_pearson_gamma_sp            <- 1
  names(pearson_status_pearson_gamma_np)     <- "Pearson"
  names(pearson_status_pearson_gamma_sp)     <- "Pearson"
  
  summ_pearson_gamma_np                      <- c(pearson_gamma, 
                                                  mean_bstrp_pearson_gamma_np, 
                                                  se_bstrp_pearson_gamma_np,
                                                  quantiles_pearson_gamma_np,
                                                  status_pearson_gamma_np,
                                                  pearson_status_pearson_gamma_np)
  summ_pearson_gamma_sp                     <- c(pearson_gamma, 
                                                 mean_bstrp_pearson_gamma_sp,
                                                 se_bstrp_pearson_gamma_sp,
                                                 quantiles_pearson_gamma_sp,
                                                 status_pearson_gamma_sp,
                                                 pearson_status_pearson_gamma_sp)
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
                                                          PhiM    = phim_mark)
    
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
  names(pearson_partial_gamma)                       <- "Est"
  mean_bstrp_pearson_partial_gamma_np                <- mean(pearson_partial_gamma_bs_np)
  mean_bstrp_pearson_partial_gamma_sp                <- mean(pearson_partial_gamma_bs_sp)
  names(mean_bstrp_pearson_partial_gamma_np)         <- "Est_Btstrp"
  names(mean_bstrp_pearson_partial_gamma_sp)         <- "Est_Btstrp"
  se_bstrp_pearson_partial_gamma_np                  <- sd(pearson_partial_gamma_bs_np)
  se_bstrp_pearson_partial_gamma_sp                  <- sd(pearson_partial_gamma_bs_sp)
  names(se_bstrp_pearson_partial_gamma_np)           <- "SE"
  names(se_bstrp_pearson_partial_gamma_sp)           <- "SE"
  quantiles_pearson_partial_gamma_np                 <- quantile(pearson_partial_gamma_bs_np, c(0.025, 0.25, 0.5, 0.75, 0.975))
  quantiles_pearson_partial_gamma_sp                 <- quantile(pearson_partial_gamma_bs_sp, c(0.025, 0.25, 0.5, 0.75, 0.975))
  status_pearson_partial_gamma_np                    <- 0
  status_pearson_partial_gamma_sp                    <- 1
  names(status_pearson_partial_gamma_np)             <- "Parametric"
  names(status_pearson_partial_gamma_sp)             <- "Parametric"
  pearson_status_pearson_partial_gamma_np            <- 2
  pearson_status_pearson_partial_gamma_sp            <- 2
  names(pearson_status_pearson_partial_gamma_np)     <- "Pearson"
  names(pearson_status_pearson_partial_gamma_sp)     <- "Pearson"
  summ_pearson_partial_gamma_np                      <- c(pearson_partial_gamma, 
                                                          mean_bstrp_pearson_partial_gamma_np,
                                                          se_bstrp_pearson_partial_gamma_np, 
                                                          quantiles_pearson_partial_gamma_np,
                                                          status_pearson_partial_gamma_np,
                                                          pearson_status_pearson_partial_gamma_np)
  summ_pearson_partial_gamma_sp                      <- c(pearson_partial_gamma,
                                                          mean_bstrp_pearson_partial_gamma_sp, 
                                                          se_bstrp_pearson_partial_gamma_sp, 
                                                          quantiles_pearson_partial_gamma_sp,
                                                          status_pearson_partial_gamma_sp,
                                                          pearson_status_pearson_partial_gamma_sp)
  #-------------------------------------------------------------------------------------------------------------
  
  
  # Return Results----------------------------------------------------------------------------------------------
  
  cat("Success, returning results ...","\n")
  
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
           Scenario            = scenario) 
  
  rownames(summ_corr) <- NULL
  
  results <- list(random_seed                 = random_seed,
                  cjs_out                     = cjs_out,
                  summ_corr                   = summ_corr,
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
                               init = NULL){
  
  
  # Run niter replicates 
  results_list <- lapply(1:niter, function(iter) execute_iteration(iter,
                                                                   scenario,
                                                                   PM,
                                                                   PF,
                                                                   PhiF,
                                                                   PhiM,
                                                                   gam_true,
                                                                   rho_true,
                                                                   n_pop,
                                                                   k,
                                                                   init = NULL))
  
  
  
  # Construct Summaries of Probs and Correlations
  summary_corr <- do.call(rbind, lapply(1:niter, function(iter) results_list[[iter]]$summ_corr))
  summary_mark <- do.call(rbind, lapply(1:niter, function(iter) results_list[[iter]]$cjs_out))
  
  # Return Results
  out <- list(results_list = results_list,
              summary_corr = summary_corr,
              summary_mark = summary_mark)
  
  return(out)
  
}

# Get Scenario Grid for correlation estimator simulation study
get_scenarios <- function(){
  n <- c(150, 350)
  k <- c(15,  30)
  PF <- c(0.45, 0.75)
  PM <- c(0.45, 0.75)
  PhiF <- c(0.8)
  PhiM <- c(0.8)
  rho <- sort(unique(c(0, seq(-0.1, 0.9, by = 0.25))))
  gamma <-  sort(unique(c(0,seq(-0.1, 0.9, by = 0.25))))
  
  scenario_grid          <- expand.grid(n_obs    = n,
                                        k        = k,
                                        rho_true = rho,
                                        gam_true = gamma, 
                                        PhiF     = PhiF,
                                        PhiM     = PhiM,
                                        PF       = PF)
  
  scenario_grid$PM       <- scenario_grid$PF
  scenario_grid$scenario <- 1:nrow(scenario_grid)
  return(scenario_grid)
}


# Compact ifelse for NA
convert_na <- function(x,y=0){
  ifelse(is.na(x), y, x)
}

# Get Summary Statistics from CODA MCMC object
gather_posterior_summary <- function(fit, nchains){
  
  # Summarize Results
  summ <- summary(fit)
  # nchains <- length(fit)
  
  if(nchains > 1){
    niter   <- nrow(fit[[1]])
  } else {
    niter <- nrow(fit)
  }
  
  # Put results into dataframe
  summ_stats <- as.data.frame(summ$statistics) %>%
    tibble::rownames_to_column(var = "Parameter")
  summ_quant <- as.data.frame(summ$quantiles)  %>% 
    tibble::rownames_to_column(var = "Parameter")
  
  if(nchains > 1){
    rhat <- gelman.diag(fit,multivariate = FALSE) 
    rhat_df <- rhat$psrf %>% 
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Parameter") %>% 
      rename("rhat" = "Point est.", 
             rhat_upper = "Upper C.I.")
  }
  
  effective_vector <- effectiveSize(fit) 
  effective_df <- data.frame(Parameter = names(effective_vector),
                             n_eff = unname(effective_vector)) %>% 
    mutate(normalized_n_eff = n_eff/(nchains * niter))
  
  summ_stats <- as.data.frame(summ$statistics) %>%
    tibble::rownames_to_column(var = "Parameter")
  
  summ_quant <- as.data.frame(summ$quantiles)  %>%
    tibble::rownames_to_column(var = "Parameter")
  
  post_stats <- suppressWarnings(inner_join(summ_stats, 
                                            summ_quant, 
                                            by = "Parameter") %>% 
                                   mutate(unique_pars      = as.factor(gsub(Parameter, pattern = "\\[.*]", replacement = "")),
                                          par_index        = convert_na(as.double(sub("^.*\\[(.*?)\\]","\\1",Parameter)))) %>% 
                                   inner_join(effective_df, by = "Parameter") %>%
                                   mutate(mcse_n_eff = SD/sqrt(n_eff)) %>% 
                                   rename(mcse = "Naive SE",
                                          mcse_ts = "Time-series SE"))
  
  if(nchains > 1){
    post_stats <- post_stats %>% 
      inner_join(rhat_df, by = "Parameter")
  }
  
  # Return Results
  return(post_stats)
}
