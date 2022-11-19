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
run_cjs_model_mark <- function(cjs_data       , 
                               PhiF     = NULL,
                               PhiM     = NULL,
                               PF       = NULL,
                               PM       = NULL,
                               Iter     = NULL,
                               scenario = NULL){
  
  #Choose Appropriate CJS Model Settings
  phi.grp     <- list(formula = ~sex)
  p.grp       <- list(formula = ~sex)
  phi.name    <- c("phi.f","phi.m")
  p.name      <- c("p.f","p.m")
  param.names <- c(phi.name,p.name)
  
  # Add True Parameters if they exist
  phi.true   <- c(PhiF,PhiM)
  p.true     <- c(PF,PM)
  param.true <- c(phi.true,p.true)
  
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
    mutate(Truth     = param.true,
           Bias      = Est - param.true,
           Range     = UB - LB,
           Parameter = param.names,
           Iter      = Iter,
           Scenario  = scenario)
  
  mark_results <- cbind(mark.stats,gof.stats)
  #Return Output
  return(mark_results)
  
}


# Execute Simulation Study 
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
  cat("Generating pair-swap mark-recapture data...","\n")
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
  
  # Generate One set of Data
  ps_data <- sim_dat(param_list) # pair-swap data
  cjs_data <- format_to_cjs(ps_data)

  # Run CJS Model MARK -----------------------------------------------------------------------------------------
  cat("Estimate standard CJS estimates with program MARK...","\n")
  cjs_out <- run_cjs_model_mark(cjs_data = cjs_data,
                                PhiF     = PhiF,
                                PhiM     = PhiM,
                                PF       = PF,
                                PM       = PM,
                                Iter     = iter,
                                scenario = scenario)
  
  pred_probs <- cjs_out$Est
  
  #-------------------------------------------------------------------------------------------------------------
  
  # Compute Recapture Correlation Estimate----------------------------------------------------------------------
  cat("Estimating recapture correlation rho...","\n")
  rho <- compute_recapture_correlation(ps_data = ps_data, 
                                       PF      = pred_probs[3],
                                       PM      = pred_probs[4])
  names(rho) <- "Est"
  
  # Bootstrap To Estimate SE 
  cat("Bootstrapping to get standard error estimates of rho...","\n")
  rho_bs <- compute_bootstrap_estimates_recapture_correlation(ps_data = ps_data,
                                                              iter    = 10000,
                                                              PF      = pred_probs[3],
                                                              PM      = pred_probs[4])
  # Collect Results
  mean_bstrp_rho      <- mean(rho_bs)
  names(mean_bstrp_rho) <- "Est_Btstrp"
  se_bstrp_rho        <- sd(rho_bs)
  names(se_bstrp_rho) <- "SE"
  quantiles_rho       <- quantile(rho_bs, c(0.025, 0.5, 0.75, 0.975))
  summ_rho <- c(rho, mean_bstrp_rho, se_bstrp_rho, quantiles_rho)
  #-------------------------------------------------------------------------------------------------------------
  
  # Compute Survival Correlation Estimate-----------------------------------------------------------------------
  cat("Estimating recapture correlation rho...","\n")
  gamma <- compute_survival_correlation(ps_data = ps_data,
                                        PFM     = compute_jbin_cjs(prob.f = pred_probs[3],
                                                                   prob.m = pred_probs[4],
                                                                   corr   = rho)$prob.mf,
                                        PhiF    = pred_probs[1],
                                        PhiM    = pred_probs[2])
  names(gamma) <- "Est"
  
  # Bootstrap To Estimate SE 
  cat("Bootstrapping to get standard error estimates of gamma...","\n")
  gamma_bs <- compute_bootstrap_estimates_survival_correlation(ps_data               = ps_data,
                                                               iter                  = 10000,
                                                               recapture_correlation = NULL,
                                                               PF                    = pred_probs[3],
                                                               PM                    = pred_probs[4],
                                                               PhiF                  = pred_probs[1],
                                                               PhiM                  = pred_probs[2])
  
  # Collect Results
  mean_bstrp_gamma      <- mean(gamma_bs)
  names(mean_bstrp_gamma) <- "Est_Btstrp"
  se_bstrp_gamma        <- sd(gamma_bs)
  names(se_bstrp_gamma) <- "SE"
  quantiles_gamma       <- quantile(gamma_bs, c(0.025, 0.5, 0.75, 0.975))
  summ_gamma <- c(gamma, mean_bstrp_gamma, se_bstrp_gamma, quantiles_gamma)
  #-------------------------------------------------------------------------------------------------------------
  
  # Return Results----------------------------------------------------------------------------------------------
  
  cat("Success, returning results ...","\n")
  
  # Gather Correlation Results
  summ_corr           <- as.data.frame(rbind(summ_rho, summ_gamma))
  param_true          <- c(rho_true, gam_true)
  summ_corr$Parameter <- c("rho","gamma")
  summ_corr           <- summ_corr[,c("Parameter","Est", "Est_Btstrp", "SE", "2.5%","50%","75%","97.5%")]
  summ_corr           <- summ_corr %>% 
    mutate(Truth        = param_true,
           Bias         = Est - param_true,
           Bias_Btstrp1 = Est_Btstrp - param_true,
           Bias_Btstrp2 = Est - Est_Btstrp,
           Range95     = `97.5%` - `2.5%`,
           Range50     = `75%`   - `50%`,
           Iter        = iter,
           Scenario    = scenario)
  
  rownames(summ_corr) <- NULL
  
  results <- list(random_seed = random_seed,
                  cjs_out     = cjs_out,
                  summ_corr   = summ_corr,
                  rho_bs      = unname(rho_bs),
                  gamma_bs    = unname(gamma_bs))
  
  return(results)
}

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

# Get Scenario Grid
# Convenient Function to get scenarios explored in this manuscript
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
