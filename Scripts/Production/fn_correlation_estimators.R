# Recapture Correlation Functions -----------------------------------------------------------------------------------------

# Partial log-likelihood function 
partial_likelihood_recapture_correlation <- function(pars,
                                                     PF,
                                                     PM,
                                                     recap_f_filter,
                                                     recap_m_filter){
  recap_dist <- compute_jbin_cjs(prob.f = PF, 
                                 prob.m = PM, 
                                 corr   = pars)
  
  lli <- log(c(recap_dist$prob.00,
               recap_dist$prob.m0,
               recap_dist$prob.f0,
               recap_dist$prob.mf))
  
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
  pairs_seen <- !is.na(ps_data$apairs_f)
  
  pairs_seen_index_list <- list()
  pairs_mask <- matrix(F, nrow = nf, ncol = k-2)
  
  for(t in 2:(k-1)){
    pairs_seen_index_list[[t-1]] <- which(pairs_seen[1:nf,t] & pairs_seen[1:nf,t-1] & pairs_seen[1:nf,t+1])
    pairs_mask[pairs_seen_index_list[[t-1]] ,t-1] <- (apairs_f[pairs_seen_index_list[[t-1]] ,t] == apairs_f[pairs_seen_index_list[[t-1]] ,t-1]) & 
      (apairs_f[pairs_seen_index_list[[t-1]] ,t] == apairs_f[pairs_seen_index_list[[t-1]] ,t+1]) & 
      (apairs_f[pairs_seen_index_list[[t-1]] ,t] != nm+1)
  }
  
  
  male_index <- apairs_f[1:nf,2:(k-1)][pairs_mask]
  times <- t((2:(k-1))* t(pairs_mask))[pairs_mask]
  recap_f_filter <- recap_f[1:nf, 2:(k-1)][pairs_mask]
  recap_m_filter <- sapply(1:length(male_index), function(x) recap_m[male_index[x], times[x]])
  
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

generate_bootstrap_replicates_recapture <- function(ps_data, iter){
  
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
  
  pairs_seen_index_list <- list()
  pairs_mask <- matrix(F, nrow = nf, ncol = k-2)
  
  for(t in 2:(k-1)){
    pairs_seen_index_list[[t-1]] <- which(pairs_seen[1:nf,t] & pairs_seen[1:nf,t-1] & pairs_seen[1:nf,t+1])
    pairs_mask[pairs_seen_index_list[[t-1]] ,t-1] <- (apairs_f[pairs_seen_index_list[[t-1]] ,t] == apairs_f[pairs_seen_index_list[[t-1]] ,t-1]) & 
      (apairs_f[pairs_seen_index_list[[t-1]] ,t] == apairs_f[pairs_seen_index_list[[t-1]] ,t+1]) & 
      (apairs_f[pairs_seen_index_list[[t-1]] ,t] != nm+1)
  }
  
  index_viable_f <- which(rowSums(pairs_mask) != 0)
  size <- length(index_viable_f)
  
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

compute_bootstrap_estimates_recapture_correlation <- function(ps_data, 
                                                              iter,
                                                              PF,
                                                              PM){
  
  bootstrap_data_replicates <- generate_bootstrap_replicates_recapture(ps_data, iter)
  rho <- compute_recapture_correlation_simulation(bootstrap_data_replicates, rep(PF, iter), rep(PM, iter))
  return(rho)
}

# SURVIVAL CORRELATION -------------------------------------------------------------------------------------------------

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
  N <- 0
  n <- 0
  
  # Unpack data 
  nm              <- ps_data$nm
  nf              <- ps_data$nf
  k               <- ps_data$k
  first_capture_f <- ps_data$first_capture_f
  recap_f         <- ps_data$recap_f
  recap_m         <- ps_data$recap_m
  apairs_f        <- ps_data$apairs_f
  
  pairs_known <- !is.na(apairs_f[1:nf,1:k])
  pairs_taken <- (apairs_f[1:nf,1:k] != (nm+1)) 
  
  for(i in 1:nf){
    if(first_capture_f[i] == k) next
    
    for(j in (first_capture_f[i]+1):k){
      if(pairs_known[i,j-1]){
        if(pairs_taken[i,j-1]){
          if(recap_f[i,j-1] == 1 & recap_m[apairs_f[i,j-1],j-1] == 1){
            if(pairs_known[i,j]){
              if((apairs_f[i,j] == apairs_f[i,j-1])){
                if(recap_f[i,j] == 1 & recap_m[apairs_f[i,j],j] == 1){
                  n <- n + 1
                  N <- N +1
                } else {
                  N <- N+1
                } 
              } else {
                N <- N + 1
              }
            }else {
              N <- N + 1
            } 
          }
        }
      }
    }
  }
  
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
  
  for(i in 1:length(recapture_correlations)){
    PFM_list[[i]] <- compute_jbin_cjs(PF[i],PM[i],recapture_correlations[i])$prob.mf 
  }
  
  PFM <- unlist(PFM_list)
  
  survival_correlations <- lapply(1:length(ps_data_list), 
                                  function(i) compute_survival_correlation(ps_data = ps_data_list[[i]],
                                                                           PFM     = PFM[i],
                                                                           PhiF    = PhiF[i],
                                                                           PhiM    = PhiM[i]))
  return(unlist(survival_correlations))
} 




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

# Likehood Survival Correlation Approach (EXPERIMENTAL)

# # Compute Survival Correlation for One Replicate Using Derived Formula
# partial_likelihood_survival <- function(pars,
#                                         ps_data,
#                                         PF,
#                                         PM,
#                                         rho,
#                                         PhiF,
#                                         PhiM){
# 
#   # Unpack data
#   nm              <- ps_data$nm
#   nf              <- ps_data$nf
#   k               <- ps_data$k
#   first_capture_f <- ps_data$first_capture_f
#   recap_f         <- ps_data$recap_f
#   recap_m         <- ps_data$recap_m
#   apairs_f        <- ps_data$apairs_f
#   apairs_m        <- ps_data$apairs_m
#   af              <- ps_data$af
#   am              <- ps_data$am
# 
#   # Condtion on a pair having been observed together in the last jump
#   pairs_known <- !is.na(apairs_f[1:nf,1:k])
#   male_pairs_known <- !is.na(apairs_m[1:nm,1:k])
#   pairs_observed <- matrix(F, nrow = nf, ncol = k-1)
#   for(i in 1:nf){
#     for(t in 1:(k-1)){
#       if(pairs_known[i,t]){
#         pairs_observed[i,t] <- (apairs_f[i,t] != (nm+1))# *(recap_f[i,t] == 1) & (recap_m[apairs_f[i,t],t] == 1)
#       }
#     }
#   }
#   
#   # Recapute Probabilities (known)
# 
#   p <- rev(unlist(compute_jbin_cjs(PF,PM, rho)))
#   p_lli   <- log(p)
# 
#   # Survival Probabilties (Proposed)
#   phi <- rev(unlist(compute_jbin_cjs(PhiF,PhiM, pars)))
#   phi_lli  <- log(phi)
# 
#   ll <- matrix(0, nrow = nf, ncol = k)
# 
#   for(i in 1:nf){
#     
#     if(first_capture_f[i] == k) next
#     
#     for(t in (first_capture_f[i]+1):k){
#       if(pairs_observed[i,t-1]){
# 
#         if(apairs_m[apairs_f[i,t-1],t-1] != i) browser()
# 
#         # Scenarios in which pairs at t are unknown
#         if(pairs_known[i,t]){
# 
#           # Pair Seen again and is same
#           if(apairs_f[i,t] == apairs_f[i,t-1]){
# 
#             if(recap_f[i,t] == 1 & recap_m[apairs_f[i,t],t] == 1){
#               ll[i,t] <- phi_lli[4] + p_lli[4]
#               next
#             }
# 
#             if(recap_f[i,t] == 1 & recap_m[apairs_f[i,t],t] == 0){
#               ll[i,t] <- phi_lli[4] + p_lli[3]
#               next
#             }
# 
#             if(recap_f[i,t] == 0 & recap_m[apairs_f[i,t],t] == 1){
#               ll[i,t] <- phi_lli[4] + p_lli[2]
#               next
#             }
# 
#             if(recap_f[i,t] == 0 & recap_m[apairs_f[i,t],t] == 0){
#               ll[i,t] <- phi_lli[4] + p_lli[1]
#               next
#             }
# 
#             browser()
# 
#             # Pair seen but they are different
#           } else {
# 
#             if(apairs_f[i,t] == (nm+1)){
#               # next
#               ll[i,t] <- phi_lli[3] + log(PF)
#               next
#             }
# 
# 
#             if(apairs_f[i,t] != (nm+1)){
# 
#               ll[i,t] <- phi_lli[3] + p_lli[4]
#               next
#             }
# 
#             browser()
#           }
# 
#           # Scenarios in which pairs at t are unknown for female i
#         } else if(male_pairs_known[apairs_f[i,t-1],t]){
# 
#           if(apairs_m[apairs_f[i,t-1],t] == (nf+1)){
#             # next
#             ll[i,t] <- phi_lli[2] + log(PM)
#             next
#           }
# 
# 
#           if(apairs_m[apairs_f[i,t-1],t] != (nf+1)){
# 
#             ll[i,t] <- phi_lli[2] + p_lli[4]
#             next
#           }
# 
# 
#           browser()
# 
#           # If next set of pairs is missing for both male and female
#         } else {
# 
#           # Both known survival states
#           if(!is.na(af[i,t]) & !is.na(am[apairs_f[i,t-1],t])){
# 
#             if(af[i,t] == 1 & am[apairs_f[i,t-1],t] == 1){
#               # Since there is only divorce upon death, the pair should be together, just not fully observed
# 
#               if(recap_f[i,t] == 1 & recap_m[apairs_f[i,t-1],t] == 0){
#                 ll[i,t] <- phi_lli[4] + p_lli[3]
#                 next
#               }
# 
#               if(recap_f[i,t] == 1 & recap_m[apairs_f[i,t-1],t] == 1){
#                 ll[i,t] <- phi_lli[4] + p_lli[4]
#                 next
#               }
# 
#               if(recap_f[i,t] == 0 & recap_m[apairs_f[i,t-1],t] == 1){
#                 ll[i,t] <- phi_lli[4] + p_lli[2]
#                 next
#               }
# 
#               if(recap_f[i,t] == 0 & recap_m[apairs_f[i,t-1],t] == 0){
#                 ll[i,t] <- phi_lli[4] + p_lli[1]
#                 next
#               }
#             }
# 
#             browser()
#           }
# 
#           # Male unknown female known
#           if(is.na(af[i,t]) & !is.na(am[apairs_f[i,t-1],t])){
# 
#             if(am[apairs_f[i,t-1],t] == 1){
#               if(recap_m[apairs_f[i,t-1],t] == 1){
#                 ll[i,t] <- log(phi[4] * p[2] + phi[2] * p[2])
#                 next
# 
#               } else {
#                 ll[i,t] <- log(phi[4] * p[1] + phi[2] * ((1-PM) + p[1] + p[3]))
#                 next
#               }
#             }
# 
#             browser()
#           }
# 
# 
#           # Female unknown male known
#           if(!is.na(af[i,t]) & is.na(am[apairs_f[i,t-1],t])){
# 
# 
#             if(af[i,t] == 1){
#               if(recap_f[i,t] == 1){
#                 ll[i,t] <- log(phi[4] * p[3] + phi[3] * p[3])
#                 next
# 
#               } else {
#                 ll[i,t] <- log(phi[4] * p[1] + phi[3] * ((1-PF) + p[1] + p[2]))
#                 next
#               }
#             }
# 
#             browser()
#           }
# 
#           # Male unknown female unknown
#           if(is.na(af[i,t]) & is.na(am[apairs_f[i,t-1],t])){
#             ll[i,t] <- log(phi[4] * p[1] + phi[2] * ((1-PM) + p[1] + p[3]) + phi[3] * ((1-PF) + p[1] + p[2]) + phi[1])
#             next
#           }
# 
# 
#           browser()
# 
#         }
#       }
#     }
#   }
# 
#   if(any(ll > 0)) browser()
# 
#   return(-sum(ll))
# }
# 
# PhiF <- cjs_out %>% filter(Parameter == "phi.f") %>% pull(Est)
# PhiM <- cjs_out %>% filter(Parameter == "phi.m") %>% pull(Est)
# PM <- cjs_out %>% filter(Parameter == "p.m") %>% pull(Est)
# PF <- cjs_out %>% filter(Parameter == "p.f") %>% pull(Est)
# 
# # Compute all relevant parameters
# params_cjs <- compute_jbin_param_cjs(PhiF, PhiM)
# 
# # Standard Deviation of PhiF, PhiM
# sigF <- params_cjs$sig.prob.f
# sigM <- params_cjs$sig.prob.m
# 
# # Upper and Lower bounds of correlation
# gl <- params_cjs$cor_lower_bound
# gu <- params_cjs$cor_upper_bound
# # #
# # 
# # i <- 5
# # gamma <- nlminb(start            = runif(1, min = gl,max = gu),
# #                 objective        = partial_likelihood_survival,
# #                 ps_data          = ps_data,
# #                 PF               = PF,
# #                 PM               = PM,
# #                 rho              = rho_true,
# #                 PhiF             = PhiF,
# #                 PhiM             = PhiM,
# #                 lower            = gl,
# #                 upper            = gu)$par
# # print(gamma-gam_true)
# # gamma2 <- compute_survival_correlation(ps_data = ps_data,
# #                                       PFM     = compute_jbin_cjs(PF, PM, rho_true)$prob.mf,
# #                                       PhiF    = PhiF,
# #                                       PhiM    = PhiM)
# # print(gamma2-gam_true)
# # 
# # 
# # 
# # rho <- compute_recapture_correlation_simulation(ps_data_list, rep(PF,1000),  rep(PM,1000))
# # plot(density(rho))
# # sd(rho)
# # mean(rho)
# # gamma <- compute_survival_correlation_simulation(ps_data_list, rho, rep(PF,1000),  rep(PM,1000), rep(PhiF,1000),  rep(PhiM,1000))
# # plot(density(gamma))
# # sd(gamma)
# # mean(gamma)
# # 
# # gamma2 <- rep(0,10)
# # 
# # for(i in 1:10){
#   gamma2[i] <- nlminb(start            = runif(1, min = gl,max = gu),
#                       objective        = partial_likelihood_survival,
#                       ps_data          = ps_data_list[[i]],
#                       PF               = PF,
#                       PM               = PM,
#                       rho              = rho[i],
#                       PhiF             = PhiF,
#                       PhiM             = PhiM,
#                       lower            = gl,
#                       upper            = gu)$par
# 
# # }
# # 
# # plot(density(gamma2))
# # c(mean(gamma2),sd(gamma2))
# # lines(density(gamma))
# # abline(v = gam_true, col = "green")
# abline(v = mean(gamma2), col = "blue")
# abline(v = mean(gamma), col = "red")
# 
# abline(v = median(gamma2), lty = 2, col = "blue")
# abline(v = median(gamma), lty = 2, col = "red")
