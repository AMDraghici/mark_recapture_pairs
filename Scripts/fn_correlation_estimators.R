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
                                          PM,
                                          model){
  
  # Make model type lower case
  model <- tolower(trimws(model))
  
  # Check for any bad names
  if(!model %in% (c("likelihood", "full_pearson", "partial_pearson"))) stop("Need to specify model as: likelihood, full_pearson or partial_pearson")
  
  # Unpack data
  rl               <- compute_jbin_param_cjs(PF,PM)$cor_lower_bound + 1e-6
  ru               <- compute_jbin_param_cjs(PF,PM)$cor_upper_bound - 1e-6 
  nf               <- ps_data$nf
  k                <- ps_data$k
  nm               <- ps_data$nm
  first_capture_f  <- ps_data$first_capture_f
  recap_f          <- ps_data$recap_f
  recap_m          <- ps_data$recap_m
  apairs_f         <- ps_data$apairs_f
  
  if(dim(apairs_f)[1] != nf & dim(apairs_f)[2] != k) browser()
  
  # Compute filters for likelihood function 
  pairs_seen <- !is.na(apairs_f) # check if NA
  
  if(dim(pairs_seen)[1] != nf & dim(pairs_seen)[2] != k) browser()
  
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

  # Correlation Calculations break on boundaries
  if(PF == 0|PM == 0|PF == 1|PM == 1){
    return(list(rho       = 0, 
                n_eff_rho = length(recap_f_filter)))
  }
  
  if(sum(pairs_mask * 1) == 0){
    rho <- 0
    # Return correlation and effective sample size
    return(list(rho       = rho, 
                n_eff_rho = length(recap_f_filter)))
  } else {
    
    if(model == "likelihood"){
      # Compute rho using non-linear in-box optimization 
      rho <- nlminb(start            = runif(1, min = rl,max = ru),
                    objective        = partial_likelihood_recapture_correlation,
                    PF               = PF,
                    PM               = PM,
                    recap_f_filter   = recap_f_filter,
                    recap_m_filter   = recap_m_filter,
                    lower            = rl,
                    upper            = ru)$par
    } else if(model == "partial_pearson"){
      # 
      # sdf <- sqrt(PF * (1-PF))
      # sdm <- sqrt(PM * (1-PM))
      # PFM <- sum(recap_f_filter * recap_m_filter)/length(recap_f_filter) 
      # 
      # rho <- (PFM - PF*PM)/(sdf * sdm)
      # 
      dev_f <- (recap_f_filter - PF)
      dev_m <- (recap_m_filter - PM)
      rho <- sum(dev_f * dev_m)/sqrt(sum(dev_f^2) * sum(dev_m^2))
    } else if(model == "full_pearson"){
      rho <- cor(recap_f_filter, recap_m_filter)
      if(is.na(rho)) rho <- 0
    }
  }
  
  rho <- ifelse(rho >= ru, ru, ifelse(rho <= rl, rl, rho))
  
  # Return correlation and effective sample size
  return(list(rho       = rho, 
              n_eff_rho = length(recap_f_filter),
              recap_f   = recap_f_filter,
              recap_m   = recap_m_filter,
              joint_recap = 1 + 2 * recap_f_filter + recap_m_filter))
}

# Bootstrapping for Error Estimates of One Replicate
bootstrap_dataset_recapture <- function(recap_f,
                                        recap_m,
                                        apairs_f,
                                        first_capture_f,
                                        nf,
                                        nm,
                                        k,
                                        blocks,
                                        block_matrix,
                                        PF,
                                        PM,
                                        rho,
                                        parametric){
  
  # sample without replacement and build new recapture dataset
  block_sample   <- sample(blocks, replace = T, size = blocks)
  
  row_list <- list()
  col_list <- list()
  
  for(b in 1:blocks){
    sample_index_b <-   which(block_matrix == block_sample[b])
    row_list[[b]] <-  1 + (sample_index_b-1) %% (nf)
    col_list[[b]] <- ceiling(sample_index_b/(nf)) +1
  }
  
  row_index <- unlist(row_list)
  col_index <- unlist(col_list)
  
  num_samples <- length(row_index)
  
  if(parametric){
    
    rcat <- function(n, prob){
      y <- sapply(1:n, function(i) which(rmultinom(1,1,prob)==1))
      return(y)
    }
    
    # Convert from joint distn to marginal F
    convert_to_recap_f <- function(x){
      ifelse(x == 4|x == 3, 1, 0)
    }
    # Convert from joint distn to marginal M
    convert_to_recap_m <- function(x){
      ifelse(x == 4|x == 2, 1, 0)
    }
    
    replicate_joint <- rcat(n = num_samples, rev(unlist(compute_jbin_cjs(PF,PM,rho))))
    
    return(list(recap_f_filter         = convert_to_recap_f(replicate_joint),
                recap_m_filter         = convert_to_recap_m(replicate_joint)))
    
  } else {
    recap_f_filter <- rep(0,num_samples)
    recap_m_filter <- rep(0,num_samples)
    
    for(i in 1:num_samples){
      recap_f_filter[i] <- recap_f[row_index[i], col_index[i]]
      recap_m_filter[i] <- recap_m[apairs_f[row_index[i], col_index[i]], col_index[i]]
      
    }
    
    return(list(recap_f_filter         = recap_f_filter,
                recap_m_filter         = recap_m_filter))
  }
}

#Run bootstrap simulation
generate_bootstrap_replicates_recapture <- function(ps_data,
                                                    iter,
                                                    PF,
                                                    PM,
                                                    rho,
                                                    parametric,
                                                    use_block){
  
  # Unpack data
  nf               <- ps_data$nf
  k                <- ps_data$k
  nm               <- ps_data$nm
  first_capture_f  <- ps_data$first_capture_f
  recap_f          <- ps_data$recap_f
  recap_m          <- ps_data$recap_m
  apairs_f         <- ps_data$apairs_f
  
  if(dim(apairs_f)[1] != nf & dim(apairs_f)[2] != k) browser()
  
  # Compute filters for likelihood function
  pairs_seen <- !is.na(apairs_f) # check if NA
  
  if(dim(pairs_seen)[1] != nf & dim(pairs_seen)[2] != k) browser()
  
  # Define filter matrix objects
  pairs_seen_index_list <- list()
  pairs_mask <- matrix(F, nrow = nf, ncol = k-2)
  
  # Check if partnered at time before and after this one (we condition on first/last paired captures for this estimator)
  # Also condition on singles (dummy status nm+1)
  for(t in 2:(k-1)){
    pairs_seen_index_list[[t-1]] <- which(pairs_seen[1:nf,t] & pairs_seen[1:nf,t-1]  & pairs_seen[1:nf,t+1])
    pairs_mask[pairs_seen_index_list[[t-1]] ,t-1] <- (apairs_f[pairs_seen_index_list[[t-1]] ,t] == apairs_f[pairs_seen_index_list[[t-1]] ,t-1]) &
      (apairs_f[pairs_seen_index_list[[t-1]] ,t] == apairs_f[pairs_seen_index_list[[t-1]] ,t+1]) &
      (apairs_f[pairs_seen_index_list[[t-1]] ,t] != nm+1)
  }
  
  block <- 1
  block_matrix <- matrix(0, nrow = nf, ncol = k-2)
  
  # Construct Blocks
  for(i in 1:nf){
    if(!any(pairs_mask[i,1:(k-2)])) next
    
    temp_index <- which(pairs_mask[i,1:(k-2)])
    
    last_partner <- apairs_f[i, temp_index[1]+1]
    
    for(t in temp_index){
      current_partner <- apairs_f[i,t+1]
      if(current_partner != last_partner) block <- block + 1
      block_matrix[i,t] <- block
      last_partner <- apairs_f[i,t+1]
    }
    
    block <- block + 1
  }
  
  # Block Reset
  if(!use_block){
    block_matrix[block_matrix != 0] <- 1:length(which(block_matrix != 0))
  }
  
  # Get viable rows (if conditions remove individual altogether no point in sampling)
  blocks <- pmax(0, max(block_matrix))
  
  if(blocks == 0){
    # Return dummy results if no blocks (shouldnt happen)
    replicates <- lapply(1:iter, function(i) return(list(recap_f_filter = rep(1,100),
                                                         recap_m_filter = rep(1,100))))
  } else {
    # Run iter replicates
    replicates <- lapply(1:iter, function(x) bootstrap_dataset_recapture(recap_f         = recap_f,
                                                                         recap_m         = recap_m,
                                                                         apairs_f        = apairs_f,
                                                                         first_capture_f = first_capture_f,
                                                                         nf              = nf,
                                                                         nm              = nm,
                                                                         k               = k,
                                                                         blocks          = blocks,
                                                                         block_matrix    = block_matrix,
                                                                         PF              = PF,
                                                                         PM              = PM,
                                                                         rho             = rho,
                                                                         parametric      = parametric))
  }
  return(replicates)
}



# Caller for bootstrapping of recapture correaltion
compute_bootstrap_estimates_recapture_correlation <- function(ps_data, 
                                                              iter,
                                                              PF,
                                                              PM,
                                                              rho,
                                                              parametric = F,
                                                              use_block  = F,
                                                              model){
  
  # Make model type lower case
  model <- tolower(trimws(model))
  
  # Check for any bad names
  if(!model %in% (c("likelihood", "full_pearson", "partial_pearson"))) stop("Need to specify model as: likelihood, full_pearson or partial_pearson")
  
  # Correlation Calculations break on boundaries
  if(PF == 0|PM == 0|PF == 1|PM == 1){
    return(rep(0,iter))
  } 
  
  if(rho == 10){
    return(rep(10, iter))
  }
  
  # Generate Bootstrap replicates 
  bootstrap_data_replicates <- generate_bootstrap_replicates_recapture(ps_data, iter, PF, PM, rho, parametric, use_block)
  
  # Get correlation bounds
  if(model == "likelihood"){
    rl               <- compute_jbin_param_cjs(PF,PM)$cor_lower_bound + 1e-6
    ru               <- compute_jbin_param_cjs(PF,PM)$cor_upper_bound - 1e-6
    rho_bs           <- sapply(1:iter, function(i) nlminb(start            = runif(1, min = rl,max = ru),
                                                          objective        = partial_likelihood_recapture_correlation,
                                                          PF               = PF,
                                                          PM               = PM,
                                                          recap_f_filter   = bootstrap_data_replicates[[i]]$recap_f_filter,
                                                          recap_m_filter   = bootstrap_data_replicates[[i]]$recap_m_filter,
                                                          lower            = rl,
                                                          upper            = ru)$par)
  } else if(model == "partial_pearson"){
    rho_bs <- sapply(1:iter, function(i) sum((bootstrap_data_replicates[[i]]$recap_f_filter - PF) * (bootstrap_data_replicates[[i]]$recap_m_filter - PM))/sqrt(sum((bootstrap_data_replicates[[i]]$recap_f_filter-PF)^2) * sum((bootstrap_data_replicates[[i]]$recap_m_filter-PM)^2)))
    
  } else if(model == "full_pearson"){
    # Run full pearson model
    rho_bs <- sapply(1:iter, function(i) cor(bootstrap_data_replicates[[i]]$recap_f_filter, bootstrap_data_replicates[[i]]$recap_m_filter))
  }
  
  # Set unknown to zero (occurs when recap_f or recap_m have no variance)
  rho_bs[is.na(rho_bs)] <- 0
  return(rho_bs)
}

# Survival Correlation Functions ----------------------------------------------------------------------------------------------------------------------------------

# Estimate Survival Correlation for One Replicate 
compute_surv_cor <- function(x, PMF, PhiF, PhiM){
  
  if(PhiF == 1|PhiM == 1) return(0)
  if(is.nan(x)|is.na(x)) return(10)
  
  # Compute all relevant parameters 
  params_cjs <- compute_jbin_param_cjs(PhiF, PhiM)
  
  # Standard Deviation of PhiF, PhiM
  sigF <- params_cjs$sig.prob.f
  sigM <- params_cjs$sig.prob.m
  
  # Upper and Lower bounds of correlation 
  gl <- params_cjs$cor_lower_bound+1e-5
  gu <- params_cjs$cor_upper_bound-1e-5
  
  y <- ((x/PMF) - (PhiF * PhiM))/(sigF * sigM)
  
  return(pmax(pmin(gu, y),gl))
}

# Compute Survival Correlation for One Replicate Using Derived Formula
compute_survival_correlation <- function(ps_data,
                                         PFM,
                                         PhiF,
                                         PhiM){
  
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
  pairs_taken[is.na(pairs_taken)] <- FALSE
  partner_grid <- pairs_known * pairs_taken
  trial_matrix <- matrix(0, nrow = nf, ncol = k)
  success_matrix <- matrix(0, nrow = nf, ncol = k)
  
  # Count Trials
  for(i in 1:nf){
    if(first_capture_f[i] == k) next
    for(t in first_capture_f[i]:k){
      if(partner_grid[i,t] == 1){
        if(recap_f[i,t] == 1 & recap_m[apairs_f[i,t], t] == 1){
          trial_matrix[i,t] <- 1
        }
      }
    }
  }
  
  # Count Success
  for(i in 1:nf){
    if(first_capture_f[i] == k) next
    for(t in first_capture_f[i]:(k-1)){
      if(trial_matrix[i,t] == 1 & trial_matrix[i,t+1] == 1){
        if(apairs_f[i,t] == apairs_f[i,t+1]){
          success_matrix[i,t] <- 1
        }
      }
    }
  }
  
  # Summarize results (dont count k as we cannot jump from there)
  n <- sum(success_matrix[,1:(k-1)])
  N <- sum(trial_matrix[,1:(k-1)])
  
  if(N == 0){
    gamma <- 10
  } else {
    # Use invariance property of MLE to compute gamma
    gamma <- compute_surv_cor(n/N, PFM, PhiF, PhiM)
  }
  
  return(list(gamma       = gamma,
              n_eff_gamma = N,
              n_success   = n,
              ybar        = n/N))
}

# generate_bootstrap_replicates_surv2 <- function(ps_data,
#                                                 PFM,
#                                                 PhiF,
#                                                 PhiM,
#                                                 gamma,
#                                                 iter){
#   
#   
#   # Unpack data
#   nm              <- ps_data$nm
#   nf              <- ps_data$nf
#   k               <- ps_data$k
#   first_capture_f <- ps_data$first_capture_f
#   recap_f         <- ps_data$recap_f
#   recap_m         <- ps_data$recap_m
#   apairs_f        <- ps_data$apairs_f
#   
#   # Filter on known pairs that are not single
#   pairs_known <- !is.na(apairs_f[1:nf,1:k])
#   pairs_taken <- (apairs_f[1:nf,1:k] != (nm+1))
#   pairs_taken[is.na(pairs_taken)] <- FALSE
#   partner_grid <- pairs_known * pairs_taken
#   trial_matrix <- matrix(0, nrow = nf, ncol = k)
#   success_matrix <- matrix(0, nrow = nf, ncol = k)
#   
#   # Count Trials
#   for(i in 1:nf){
#     if(first_capture_f[i] == k) next
#     for(t in first_capture_f[i]:k){
#       if(partner_grid[i,t] == 1){
#         if(recap_f[i,t] == 1 & recap_m[apairs_f[i,t], t] == 1){
#           trial_matrix[i,t] <- 1
#         }
#       }
#     }
#   }
#   
#   # Count Success
#   for(i in 1:nf){
#     if(first_capture_f[i] == k) next
#     for(t in first_capture_f[i]:(k-1)){
#       if(trial_matrix[i,t] == 1 & trial_matrix[i,t+1] == 1){
#         if(apairs_f[i,t] == apairs_f[i,t+1]){
#           success_matrix[i,t] <- 1
#         }
#       }
#     }
#   }
#   
#   # Summarize results (dont count k as we cannot jump from there)
#   n <- sum(success_matrix[,1:(k-1)])
#   N <- sum(trial_matrix[,1:(k-1)])
#   
#   PhiFM <- compute_jbin_cjs(PhiF, PhiM, gamma)[[1]]
#   
#   # Get bootstraps
#   gamma_bs <- replicate(iter, compute_surv_cor(rbinom(1, N, PhiFM * PFM)/N,PFM, PhiF, PhiM)) 
#   return(gamma_bs)
#   
# }

# Run bootstrap to get standard error estimates of survival
generate_bootstrap_replicates_surv <- function(ps_data, prob_prod, parametric, use_block, iter){
  
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
  pairs_taken[is.na(pairs_taken)] <- FALSE
  partner_grid <- pairs_known * pairs_taken
  trial_matrix <- matrix(0, nrow = nf, ncol = k)
  success_matrix <- matrix(0, nrow = nf, ncol = k)
  block_matrix <- matrix(0, nrow = nf, ncol = k)
  block <- 1
  
  # Count how many trials
  for(i in 1:nf){
    if(first_capture_f[i] == k) next
    for(t in first_capture_f[i]:k){
      if(partner_grid[i,t] == 1){
        if(recap_f[i,t] == 1 & recap_m[apairs_f[i,t], t] == 1){
          trial_matrix[i,t] <- 1
        }
      }
    }
  }
  
  # Count how many successes 
  for(i in 1:nf){
    if(first_capture_f[i] == k) next
    for(t in first_capture_f[i]:(k-1)){
      if(trial_matrix[i,t] == 1 & trial_matrix[i,t+1] == 1){
        if(apairs_f[i,t] == apairs_f[i,t+1]){
          success_matrix[i,t] <- 1
        }
      }
    }
  }
  
  # Construct Blocks (block is defined temporal-pair combos)
  for(i in 1:nf){
    if(first_capture_f[i] == k) next
    if(sum(trial_matrix[i,1:(k-1)]) == 0) next
    temp_index <- which(trial_matrix[i,1:(k-1)] == 1)
    last_partner <- apairs_f[i, temp_index[1]]
    for(t in temp_index){
      current_partner <- apairs_f[i, t]
      if(current_partner != last_partner) block <- block + 1
      block_matrix[i,t] <- block
      last_partner <- apairs_f[i,t]
    }
    block <- block + 1
  }
  
  # Block Reset
  if(!use_block){
    block_matrix[block_matrix != 0] <- 1:length(which(block_matrix != 0))
  }

  # Collection of Y/N
  outcome_vector <- unlist(lapply(1:nf, function(i) success_matrix[i, ][block_matrix[i,] != 0]))
  # Block Label for Outcomes
  block_vector   <- unlist(lapply(1:nf, function(i) block_matrix[i, ][block_matrix[i,] != 0]))
  # Number of replicates 
  replicates <- rep(0, iter)
  # Generate Bootstrapped Estimates of Gamma using Semi-Parametric and Non-Parametric
  if(parametric){
    for(i in 1:iter){
      # Sample blocks of outcomes based on pairs 
      block_replicate <- sample(max(block_vector), max(block_vector), replace = T)
      outcome_replicate <- unlist(lapply(1:length(block_replicate), function(t) rbinom(n    = length(which(block_vector == block_replicate[t])), 
                                                                                       size = 1, 
                                                                                       prob = prob_prod)))
      replicates[i] <- mean(outcome_replicate)
    }
    
  } else {
    # Non-Parametric Approach
    for(i in 1:iter){
      # Sample blocks of outcomes based on pairs 
      block_replicate <- sample(max(block_vector), max(block_vector), replace = T)
      # Grab outcomes from the original dataset 
      outcome_replicate <- unlist(lapply(1:length(block_replicate), function(t) outcome_vector[which(block_vector == block_replicate[t])]))
      # Convert to MLE 
      replicates[i] <- mean(outcome_replicate)
    }
  }
  
  return(replicates)
}

# Call bootstrap gamma
compute_bootstrap_estimates_survival_correlation <- function(ps_data, 
                                                             iter,
                                                             rho,
                                                             PF,
                                                             PM,
                                                             gamma,
                                                             PhiF,
                                                             PhiM,
                                                             parametric,
                                                             use_block){
  # In case of empty data
  if(gamma == 10|rho == 10) return(rep(10, iter))
  # Correlation must be zero here
  if(PhiF == 0|PhiF == 1|PhiM == 0|PhiM == 1) return(rep(0,iter))
  
  # Joint Probabilties
  PFM                       <- compute_jbin_cjs(PF, PM, rho)$prob.mf
  PhiFM                     <- compute_jbin_cjs(PhiF, PhiM, gamma)$prob.mf
  bootstrap_data_replicates <- generate_bootstrap_replicates_surv(ps_data, PhiFM, parametric, use_block, iter)
  gamma_bs <- sapply(1:iter, function(i) compute_surv_cor(x    = PFM * bootstrap_data_replicates[i],
                                                          PMF  =  PFM,
                                                          PhiF = PhiF,
                                                          PhiM = PhiM))
  return(gamma_bs)
}