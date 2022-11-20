# Recapture Correlation Functions --------------------------------------------------------------------------------------------------------------------------------

# Partial log-likelihood function for recapture correlation 
partial_likelihood_recapture_correlation <- function(pars,
                                                     PF,
                                                     PM,
                                                     recap_f_filter,
                                                     recap_m_filter){
  
  # Compute joint bernoulli density with proposed correlation
  recap_dist <- compute_jbin_cjs(prob.f = PF, 
                                 prob.m = PM, 
                                 corr   = pars)
  
  # Convert to log-likelihood for a pair 
  lli <- log(c(recap_dist$prob.00,
               recap_dist$prob.m0,
               recap_dist$prob.f0,
               recap_dist$prob.mf))
  
  # Map conditioned recapture states to log-likelihood and compute negative sum
  return(-sum(lli[1 + 2 * recap_f_filter + recap_m_filter]))
}

# Compute recapture Correlation for one data set
compute_recapture_correlation <- function(ps_data,
                                          PF, 
                                          PM){
  
  # Unpack data
  rl               <- compute_jbin_param_cjs(PF,PM)$cor_lower_bound + 1e-7
  ru               <- compute_jbin_param_cjs(PF,PM)$cor_upper_bound - 1e-7 
  nf               <- ps_data$nf
  k                <- ps_data$k
  nm               <- ps_data$nm
  first_capture_f  <- ps_data$first_capture_f
  recap_f          <- ps_data$recap_f
  recap_m          <- ps_data$recap_m
  apairs_f         <- ps_data$apairs_f
  
  # Compute filters for likelihood function 
  pairs_seen <- !is.na(ps_data$apairs_f) # check if NA
  
  # Define filter matrix objects 
  pairs_seen_index_list <- list()
  pairs_mask <- matrix(F, nrow = nf, ncol = k-2)
  
  # Check if partnered at time before and after this one (we condition on first/last paired captures for this estimator)
  # Also condition on singles (dummy status nm+1)
  for(t in 2:(k-1)){
    pairs_seen_index_list[[t-1]] <- which(pairs_seen[1:nf,t] & pairs_seen[1:nf,t-1] & pairs_seen[1:nf,t+1])
    pairs_mask[pairs_seen_index_list[[t-1]] ,t-1] <- (apairs_f[pairs_seen_index_list[[t-1]] ,t] == apairs_f[pairs_seen_index_list[[t-1]] ,t-1]) & 
                                                     (apairs_f[pairs_seen_index_list[[t-1]] ,t] == apairs_f[pairs_seen_index_list[[t-1]] ,t+1]) & 
                                                     (apairs_f[pairs_seen_index_list[[t-1]] ,t] != nm+1)
  }
  
  # Grab male recapture histories that are partnered with valid females on all times 
  male_index <- apairs_f[1:nf,2:(k-1)][pairs_mask]
  times <- t((2:(k-1))* t(pairs_mask))[pairs_mask]
  recap_f_filter <- recap_f[1:nf, 2:(k-1)][pairs_mask]
  recap_m_filter <- sapply(1:length(male_index), function(x) recap_m[male_index[x], times[x]])
  
  # Compute rho using non-linear in-box optimization 
  rho <- nlminb(start            = runif(1, min = rl,max = ru),
                objective        = partial_likelihood_recapture_correlation,
                PF               = PF,
                PM               = PM,
                recap_f_filter   = recap_f_filter,
                recap_m_filter   = recap_m_filter,
                lower            = rl,
                upper            = ru)$par
  
  return(rho)
}

# Compute recapture correlation for n replciates 
compute_recapture_correlation_simulation <- function(ps_data_list, PF, PM){
  
  # Run simulation study on list of data 
  recapture_correlations <- lapply(1:length(ps_data_list), 
                                   function(i) compute_recapture_correlation(ps_data = ps_data_list[[i]],
                                                                             PF      = PF[i],
                                                                             PM      = PM[i]))
  return(unlist(recapture_correlations))
} 

# Bootstrapping for Error Estimates of One Replicate
bootstrap_dataset_recapture <- function(recap_f,
                                        recap_m,
                                        apairs_f,
                                        first_capture_f,
                                        nf,
                                        nm,
                                        k,
                                        size,
                                        index_viable_f){
  
  # sample without replacement and build new recapture dataset
  index_sample   <- sample(index_viable_f, replace = T, size = size)
  
  return(list(recap_f         = recap_f[index_sample, 1:k],
              recap_m         = recap_m,
              apairs_f        = apairs_f[index_sample, 1:k],
              first_capture_f = first_capture_f[index_sample],
              nf              = size,
              nm              = nm, 
              k               = k))
}

# Run bootstrap simulation 
generate_bootstrap_replicates_recapture <- function(ps_data, 
                                                    iter){
  
  # Unpack data
  nf               <- ps_data$nf
  k                <- ps_data$k
  nm               <- ps_data$nm
  first_capture_f  <- ps_data$first_capture_f
  recap_f          <- ps_data$recap_f
  recap_m          <- ps_data$recap_m
  apairs_f         <- ps_data$apairs_f
  
  # Compute filters for likelihood function 
  pairs_seen <- !is.na(ps_data$apairs_f)
  
  # Define filter matrix objects 
  pairs_seen_index_list <- list()
  pairs_mask <- matrix(F, nrow = nf, ncol = k-2)
  
  # Check if partnered at time before and after this one (we condition on first/last paired captures for this estimator)
  # Also condition on singles (dummy status nm+1)
  for(t in 2:(k-1)){
    pairs_seen_index_list[[t-1]] <- which(pairs_seen[1:nf,t] & pairs_seen[1:nf,t-1] & pairs_seen[1:nf,t+1])
    pairs_mask[pairs_seen_index_list[[t-1]] ,t-1] <- (apairs_f[pairs_seen_index_list[[t-1]] ,t] == apairs_f[pairs_seen_index_list[[t-1]] ,t-1]) & 
                                                     (apairs_f[pairs_seen_index_list[[t-1]] ,t] == apairs_f[pairs_seen_index_list[[t-1]] ,t+1]) & 
                                                     (apairs_f[pairs_seen_index_list[[t-1]] ,t] != nm+1)
  }
  
  # Get viable rows (if conditions remove individual altogether no point in sampling)
  index_viable_f <- which(rowSums(pairs_mask) != 0)
  size <- length(index_viable_f)
  
  # Run iter replicates
  replicates <- lapply(1:iter, function(x) bootstrap_dataset_recapture(recap_f         = recap_f,
                                                                       recap_m         = recap_m,
                                                                       apairs_f        = apairs_f,
                                                                       first_capture_f = first_capture_f,
                                                                       nf              = nf,
                                                                       nm              = nm,
                                                                       k               = k,
                                                                       size            = size,
                                                                       index_viable_f  = index_viable_f))
  return(replicates)
}

# Caller for bootstrapping of recapture correaltion
compute_bootstrap_estimates_recapture_correlation <- function(ps_data, 
                                                              iter,
                                                              PF,
                                                              PM){
  
  bootstrap_data_replicates <- generate_bootstrap_replicates_recapture(ps_data, iter)
  rho <- compute_recapture_correlation_simulation(bootstrap_data_replicates, rep(PF, iter), rep(PM, iter))
  return(rho)
}

# Survival Correlation Functions ----------------------------------------------------------------------------------------------------------------------------------

# Estimate Survival Correlation for One Replicate 
compute_surv_cor <- function(x, PMF, PhiF, PhiM){
  
  # Compute all relevant parameters 
  params_cjs <- compute_jbin_param_cjs(PhiF, PhiM)
  
  # Standard Deviation of PhiF, PhiM
  sigF <- params_cjs$sig.prob.f
  sigM <- params_cjs$sig.prob.m
  
  # Upper and Lower bounds of correlation 
  gl <- params_cjs$cor_lower_bound
  gu <- params_cjs$cor_upper_bound
  
  y <- ((x/PMF) - (PhiF * PhiM))/(sigF * sigM)
  
  return(pmax(pmin(gu, y),gl))
}

# Compute Survival Correlation for One Replicate Using Derived Formula
compute_survival_correlation <- function(ps_data, 
                                         PFM, 
                                         PhiF, 
                                         PhiM){
  
  # Initialize Success/Failures
  N <- 0 #num trials 
  n <- 0 #num success
  
  # Unpack data 
  nm              <- ps_data$nm
  nf              <- ps_data$nf
  k               <- ps_data$k
  first_capture_f <- ps_data$first_capture_f
  recap_f         <- ps_data$recap_f
  recap_m         <- ps_data$recap_m
  apairs_f        <- ps_data$apairs_f
  
  # Filter on known pairs that are not single
  pairs_known <- !is.na(apairs_f[1:nf,1:k])
  pairs_taken <- (apairs_f[1:nf,1:k] != (nm+1)) 
  
  # Iterate on females (which map to males who paired through apairs_f)
  for(i in 1:nf){
    if(first_capture_f[i] == k) next # skip if first capture is on final sampling occasion
    # Iterate over time 
    for(j in (first_capture_f[i]+1):k){
      if(pairs_known[i,j-1]){
        if(pairs_taken[i,j-1]){
          if(recap_f[i,j-1] == 1 & recap_m[apairs_f[i,j-1],j-1] == 1){
            if(pairs_known[i,j]){
              if((apairs_f[i,j] == apairs_f[i,j-1])){
                if(recap_f[i,j] == 1 & recap_m[apairs_f[i,j],j] == 1){
                  # Recapture at t-1 and t as a pair therefore success
                  n <- n + 1
                  N <- N +1
                } else {
                  # Recaptured as a pair at t-1 but not t therefore fail 
                  N <- N+1
                } 
              } else {
                # Recaptured as a pair at t-1 but pairs changed at t therefore fail 
                N <- N + 1
              }
            }else {
              # Recaptured as a pair at t-1 but pairs unknown at t therefore fail 
              N <- N + 1
            } 
          }
        }
      }
    }
  }
  
  # Use invariance property of mle to compute gamma
  obs <- compute_surv_cor(n/N,PFM, PhiF, PhiM)
  
  return(obs)
}

# Compute survival correlation for n replciates 
compute_survival_correlation_simulation <- function(ps_data_list, 
                                                    recapture_correlations, 
                                                    PF, 
                                                    PM, 
                                                    PhiF, 
                                                    PhiM){
  
  # Get Joint Probability of Recapture
  PFM_list <- list()
  
  # Pass recapture correlation and marginals to get joint recapture probabilities
  for(i in 1:length(recapture_correlations)){
    PFM_list[[i]] <- compute_jbin_cjs(PF[i],PM[i],recapture_correlations[i])$prob.mf 
  }
  
  PFM <- unlist(PFM_list)
  
  # Compute survival correlation
  survival_correlations <- lapply(1:length(ps_data_list), 
                                  function(i) compute_survival_correlation(ps_data = ps_data_list[[i]],
                                                                           PFM     = PFM[i],
                                                                           PhiF    = PhiF[i],
                                                                           PhiM    = PhiM[i]))
  return(unlist(survival_correlations))
} 

# Bootstrap datasets
bootstrap_dataset_survival <- function(recap_f,
                                       recap_m,
                                       apairs_f,
                                       first_capture_f,
                                       nm,
                                       k,
                                       index_viable_f,
                                       size){
  
  # sample without replacement and build new recapture dataset
  index_sample <- sample(index_viable_f, replace = T, size = size)
  
  return(list(recap_f         = recap_f[index_sample, 1:k],
              recap_m         = recap_m,
              apairs_f        = apairs_f[index_sample, 1:k],
              first_capture_f = first_capture_f[index_sample],
              nf              = size,
              nm              = nm, 
              k               = k))
}

# Run bootstrap to get standard error estimates of survival
generate_bootstrap_replicates_surv <- function(ps_data, iter){
  
  # Grab Relevant Data
  recap_f         <- ps_data$recap_f
  recap_m         <- ps_data$recap_m
  apairs_f        <- ps_data$apairs_f
  first_capture_f <- ps_data$first_capture_f
  nf              <- ps_data$nf
  nm              <- ps_data$nm
  k               <- ps_data$k
  
  # Flag females with pair histories that are viable 
  viable_f <- vector(length = nf)
  
  for(i in 1:nf){
    for(j in (first_capture_f[i]+1):k){
      if(!is.na(apairs_f[i,j-1])){
        if(apairs_f[i,j-1] == (nm+1)) next
        if(recap_f[i,j-1] == 1 & recap_m[apairs_f[i,j-1],j-1] == 1){
          viable_f[i] <- T
          break
        }
      }
    }
  }
  
  size <- sum(1 * viable_f)
  index_viable_f <- which(viable_f)
  
  replicates <- lapply(1:iter, function(x) bootstrap_dataset_survival(recap_f         = recap_f,
                                                                      recap_m         = recap_m,
                                                                      apairs_f        = apairs_f,
                                                                      first_capture_f = first_capture_f,
                                                                      nm              = nm,
                                                                      k               = k,
                                                                      size            = size,
                                                                      index_viable_f  = index_viable_f))
  
  return(replicates)
}

# Call bootstrap
compute_bootstrap_estimates_survival_correlation <- function(ps_data, 
                                                             iter,
                                                             recapture_correlation = NULL,
                                                             PF,
                                                             PM,
                                                             PhiF,
                                                             PhiM){
  
  bootstrap_data_replicates <- generate_bootstrap_replicates_surv(ps_data, iter)
  
  # Can regenerate new values of Rho for Gamma estimation
  if(is.null(recapture_correlation)){
    rho <- compute_recapture_correlation_simulation(bootstrap_data_replicates, rep(PF, iter), rep(PM, iter))
  } else {
    rho <- rep(recapture_correlation, iter)
  }
 
  gamma <- compute_survival_correlation_simulation(ps_data_list           = bootstrap_data_replicates,
                                                   recapture_correlations = rho,
                                                   PF                     = rep(PF, iter),
                                                   PM                     = rep(PM, iter),
                                                   PhiF                   = rep(PhiF, iter),
                                                   PhiM                   = rep(PhiM, iter))
  return(gamma)
}