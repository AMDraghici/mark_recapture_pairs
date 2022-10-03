# Nimble Functions for Running Pair Swap Model ------------------------------------------------------------------------------------------------------------ 

# Produce vector of 1s and 0s to check for matching value in an existing vector
vectorMatch <- nimbleFunction(
  run = function(x= double(1),
                 y = double(0)){
    returnType(double(1))
    output <- 1*(y == x)
    return(output)}
)

# Compute Conditional Partnership Grid
compute_psi_cond_it <- nimbleFunction(
  run = function(psi_it            = double(1),
                 amating_f         = double(0),
                 amating_m         = double(1),
                 arepartner        = double(0),
                 male_taken_jt     = double(1),
                 hist_pairs_i      = double(0),
                 next_partner_i    = double(0),
                 current_partner_i = double(0),
                 nf                = integer(0),
                 nm                = integer(0)){
    
    returnType(double(0))
    # browser()
    # Apply Conditioning on Mating
    
    if(current_partner_i != 0){
      return(current_partner_i)
    }
    
    psi_cond_it <- psi_it
    psi_cond_it[1:nm] <- psi_cond_it[1:nm] * amating_f * amating_m[1:nm]
    psi_cond_it[nm+1] <- psi_cond_it[nm+1] * (1-amating_f)
    
    # If male and female are repartnering 
    if(arepartner == 1){
      psi_cond_it[1:(nm+1)] <- 0
      psi_cond_it[hist_pairs_i] <- 1
      return(hist_pairs_i)
    }
    
    if(arepartner == 0){
      psi_cond_it[1:(nm+1)] <- psi_cond_it[1:(nm+1)] * (1-male_taken_jt[1:(nm+1)])
      
      if(hist_pairs_i != (nm+1)){
        psi_cond_it[hist_pairs_i] <- 0
      }
      
      if(psi_cond_it[next_partner_i] == 1){
        return(next_partner_i) 
      }
    }
    
    # If for some reason all are zero (competing partnerships) force single 
    if(sum(psi_cond_it[1:(nm+1)]) == 0){
      psi_cond_it[nm+1] <- 1
      return(nm+1)
    }
    
    return(nm+1)
  }
)

compute_prob_repartner_it <- nimbleFunction(
  run = function(beta         = double(0),
                 amating_f    = double(0),
                 amating_m    = double(0),
                 hist_pair_i  = double(0),
                 na_repartner = double(0),
                 psi          = double(1),
                 nm           = integer(0)){
    returnType(double(0))
    prob_repartner_it <- beta * amating_f * amating_m * psi[hist_pair_i]
    
    # If previous partner is unknown then arepartner is na
    # If previous partner == current partner then arepartner MUST be 1
    mask1 <- 1 * (sum(psi[1:(nm)]) == 1)
    mask2 <- 1 * (psi[hist_pair_i] == 1)
    mask3 <- 1 * (na_repartner == 1)
    mask4 <- 1* (prob_repartner_it != 0)
    
    if(mask1 & mask2 & mask3 & mask4){
      prob_repartner_it <- 1
    }
    
    if(hist_pair_i == (nm+1)){
      prob_repartner_it <- 0
    }
    
    return(prob_repartner_it)
  }
)

assign_partner <- nimbleFunction(
  run = function(psi = double(1),
                 nm = integer(0)){
    returnType(double(0))
    possible_partners <- 1:(nm+1)
    partner_id <- inprod(psi[1:(nm+1)], possible_partners[1:(nm+1)])
    return(partner_id)
  }
)


# BUGS/JAGS Code
nimble_ps_model <- nimbleCode({
  
  
  # Conditional Partnership/Survival Steps -------------------------------------------------------------------------------------------------
  
  # 1a. Mating at time t-----------------------------------------------------------------------------------------------------------
  # Female Recruitment
  # for(i in 1:nf){
  #   for(t in (first_capture_f[i]):(k+1)){
  #     # amating_f[i,t] ~ dbern(delta * af[i,t])
  #     amating_f[i,t] <- af[i,t]
  #   }
  # }
  # 
  # # Male Recruitment
  # for(j in 1:nm){
  #   for(t in (first_capture_m[j]):(k+1)){
  #     # amating_m[j,t] ~ dbern(delta * am[j,t])
  #     amating_m[j,t] <- am[j,t]
  #   }
  # }
  
  # 1b. Repartnership at time t---------------------------------------------------------------------------------------------------
  for(i in 1:nf){
    prob_repartner_it[i,1:first_capture_f[i]] <- 0
    for(t in (first_capture_f[i]+1):k){
      prob_repartner_it[i,t] <- compute_prob_repartner_it(beta,
                                                          af[i,t],
                                                          am[apairs_f[i,t-1],t],
                                                          apairs_f[i,t-1],
                                                          na_repartner[i,t],
                                                          psi[i,1:(nm+1),t],
                                                          nm)
      
      arepartner[i,t] ~ dbern(prob_repartner_it[i,t])
    }
  }
  
  for(j in 1:nm){
    male_taken_jt[j, 1:first_capture_m[j]] <- 0
    for(t in (first_capture_m[j]+1):k){
      male_taken_jt[j,t] <- sum(vectorMatch(apairs_f[1:nf,t-1],j)*arepartner[1:nf,t])
    }
  }
  
  male_taken_jt[nm+1, 1:k] <- 0
  
  # 1b. Decision to Mate at t-----------------------------------------------------------------------------------------------------------------
  for(i in 1:nf){
    
    # Sample Female Mate 
    apairs_f[i,first_capture_f[i]] <- compute_psi_cond_it(psi[i,1:(nm+1),first_capture_f[i]],
                                                          af[i,first_capture_f[i]],
                                                          am[1:nm,first_capture_f[i]],
                                                          0,
                                                          male_taken_jt[1:(nm+1),first_capture_f[i]],
                                                          (nm+1),
                                                          next_partner_matrix[i,first_capture_f[i]],
                                                          current_partner_matrix[i, first_capture_f[i]],
                                                          nf,
                                                          nm)
    
    single_female[i,first_capture_f[i]] <- equals(apairs_f[i,first_capture_f[i]], nm + 1)
    
    
    for(t in (first_capture_f[i]+1):k){
      
      apairs_f[i,t] <- compute_psi_cond_it(psi[i,1:(nm+1),t],
                                           af[i,t],
                                           am[1:nm,t],
                                           arepartner[i,t],
                                           male_taken_jt[1:(nm+1),t],
                                           apairs_f[i,t-1],
                                           next_partner_matrix[i,t],
                                           current_partner_matrix[i, t],
                                           nf,
                                           nm)
      
      
      single_female[i,t] <- equals(apairs_f[i,t], nm + 1)
    }
  }
  
  # 3. Joint Survival [t-1,t) ---------------------------------------------------------------------------------------------------------------
  
  # Marginal Survival Event for Males in the Population (P[Y^M_T]) 
  for(j in 1:nm){
    for(t in first_capture_m[j]:(k-1)){
      am[j,t+1] ~ dbern(PhiM * am[j,t])
    }
  }
  
  # Draw conditional Survival Event
  for(i in 1:nf){
    for(t in first_capture_f[i]:(k-1)){
      af[i, t+1] ~ dbern((single_female[i,t] * PhiF +
                        (1-single_female[i,t]) * (am[apairs_f[i,t],t+1] * Phif_M1 + (1- am[apairs_f[i,t],t+1]) * Phif_M0)) * af[i,t])
      # af[i, t+1] ~ dbern(PhiF * af[i,t])
    }
  }
  
  # 4. Joint Recapture --------------------------------------------------------------------------------------------------------------------
  
  # Marginal Recapture Event for Males in the Population (P[X^M_T])
  for(j in 1:nm){
    for(t in (first_capture_m[j]+1):(k)){
      recap_m[j,t] ~ dbern(PM * am[j,t])
    }
  }
  
  # Draw Recapture Probability
  for(i in 1:nf){
    for(t in (first_capture_f[i]+1):(k)){
      recap_f[i, t] ~ dbern((single_female[i,t] * PF +
                            (1-single_female[i,t]) * (recap_m[apairs_f[i,t],t] * Pf_M1 + (1- recap_m[apairs_f[i,t],t]) * Pf_M0)) * af[i,t])
      # recap_f[i, t] ~ dbern(PF * af[i,t])
    }
  }
  
  # 5. Prior Distributions-------------------------------------------------------------------------------------------------------------------
  
  # Mating or not
  # delta ~ dbeta(2,2)
  beta ~ dbeta(2,2)
  
  # Survival Terms
  
  # ### Derived Parameters ####
  # Phi_Vector_f[[1]] <- PhiF
  # Phi_Vector_f[[2]] <- Phif_M1
  # Phi_Vector_f[[3]] <- Phif_M0
  
  # Conditional Probabilities from Female Perspective
  Phif_M1 <- Phifm/(PhiM)
  Phif_M0 <- Phif0/(1-PhiM)
  
  #Joint Survival probabilities for paired individuals
  Phi00 <- 1 - PhiF - PhiM + Phifm
  Phif0 <- PhiF - Phifm
  Phim0 <- PhiM - Phifm
  Phifm <- gamma*sig.PhiF*sig.PhiM + PhiF*PhiM
  
  ###Binomial SD for survival 
  sig.PhiF <- sqrt(PhiF*(1-PhiF))
  sig.PhiM <- sqrt(PhiM*(1-PhiM))
  
  ##Correlation (with FH bounds)
  constraint_data[2] ~ dconstraint(gamma <= gu & gamma >= gl)
  gamma <- 2*raw_gamma - 1
  raw_gamma ~ dbeta(2,2)
  # gamma <- (gu-gl) * raw_gamma + gl
  # raw_gamma ~ dbeta(1,1)
  
  # Bounds for Correlation
  
  # Survival Rates (Gamma)
  gu <-  min(sqrt(OR.Phi), 1/sqrt(OR.Phi)) 
  gl <- -min(sqrt(OP.Phi), 1/sqrt(OP.Phi)) 
  
  # Odds Ratio and Product of Survival Rates
  OP.Phi <- odds.PhiF*odds.phiM
  OR.Phi <- odds.PhiF/odds.phiM
  
  ### Odds of Recapture Rates
  odds.PhiF <- PhiF/(1 - PhiF)
  odds.phiM <- PhiM/(1 - PhiM)
  
  ##Survival Rates M/F
  PhiF ~ dbeta(1,1)
  PhiM ~ dbeta(1,1)
  
  # Recapture Terms
  ### Derived Parameters ####
  # # Conditional Probabilities
  # P_Vector_f[[1]] <- PF
  # P_Vector_f[[2]] <- Pf_M1
  # P_Vector_f[[3]] <- Pf_M0
  
  Pf_M1 <- Pfm/(PM)
  Pf_M0 <- Pf0/(1-PM)
  
  #Joint Capture probabilities for paired individuals
  P00 <- 1 - PF - PM + Pfm
  Pf0 <- PF - Pfm
  Pm0 <- PM - Pfm
  Pfm <- rho*sig.PF*sig.PM + PF*PM
  
  ###Binomial SD for recapture 
  sig.PF <- sqrt(PF*(1-PF))
  sig.PM <- sqrt(PM*(1-PM))
  
  ##Correlation using four parameter beta (with FH bounds)
  constraint_data[1] ~ dconstraint(rho <= ru & rho >= rl)
  rho <- 2 * raw_rho - 1
  raw_rho ~ dbeta(2,2)
  # rho <- (ru-rl) * raw_rho + rl
  # raw_rho ~ dbeta(1,1)
  # Bounds for Correlation
  
  # Recapture Rates (Rho)
  ru <-  min(sqrt(OR.P), 1/sqrt(OR.P)) 
  rl <- -min(sqrt(OP.P), 1/sqrt(OP.P)) 
  
  # Odds Ratio and Product of Recapture Rates
  OP.P <- odds.PF*odds.PM
  OR.P <- odds.PF/odds.PM
  
  ### Odds of Survival and Recapture Rates
  odds.PF <- PF/(1 - PF)
  odds.PM <- PM/(1 - PM)
  
  # Recapture Rates M/F
  PF ~ dbeta(1,1)
  PM ~ dbeta(1,1) 
  
})


# Generating Initial Values
generate_nimble_init_pairs <- function(ps_data){
  # set.seed(pi)
  #Unpack Variables -----------------------------------------------------------------
  # Indexes
  k               <- ps_data$k
  nf              <- ps_data$nf
  nm              <- ps_data$nm
  psi             <- ps_data$psi
  first_capture_f <- ps_data$first_capture_f
  first_capture_m <- ps_data$first_capture_m  
  next_partner_matrix <- ps_data$next_partner_matrix
  current_partner_matrix <- ps_data$apairs_f
  current_partner_matrix[is.na(current_partner_matrix)] <- 0 
  # CR data with missing components
  amating_f     <- ps_data$amating_f
  amating_m     <- ps_data$amating_m
  # apairs_f      <- ps_data$apairs_f
  apairs_f      <- matrix(NA, nrow = ps_data$nf, ncol = ps_data$k)
  af            <- ps_data$af
  am            <- ps_data$am
  recap_m       <- ps_data$recap_m
  recap_f       <- ps_data$recap_f
  arepartner    <- ps_data$arepartner
  na_repartner  <- ps_data$na_repartner
  
  # Recapture Prob and Correlation -------------------------------------------------
  PM <- rbeta(1,1,1)
  PF <- rbeta(1,1,1)
  
  ### Odds of Recapture Rates
  odds.PF <- PF/(1 - PF)
  odds.PM <- PM/(1 - PM)
  
  # Odds Ratio and Product of Recapture Rates
  OP.P <- odds.PF*odds.PM
  OR.P <- odds.PF/odds.PM
  
  # Recapture Rates (Rho)
  ru <-  min(sqrt(OR.P), 1/sqrt(OR.P))
  rl <- -min(sqrt(OP.P), 1/sqrt(OP.P))
  
  ##Correlation using four parameter beta (with FH bounds)
  raw_rho <-  -rl/(ru-rl)#rbeta(1,1,1)
  rho <- (ru - rl)*raw_rho + rl
  raw_rho <- (rho+1)/2
  
  ###Binomial SD for recapture
  sig.PF <- sqrt(PF*(1-PF))
  sig.PM <- sqrt(PM*(1-PM))
  
  #Joint Capture probabilities for paired individuals
  Pfm <- rho*sig.PF*sig.PM + PF*PM
  P00 <- 1 - PF - PM + Pfm
  Pf0 <- PF - Pfm
  Pm0 <- PM - Pfm
  
  # Conditional Probabilities from Female Perspective
  Pf_M1 <- Pfm/(PM)
  Pf_M0 <- Pf0/(1-PM)
  
  # Survival Prob and Correlation -------------------------------------------------
  PhiM <- rbeta(1,1,1)
  PhiF <- rbeta(1,1,1)
  
  ### Odds of Survival Rates
  odds.PhiM <- PhiM/(1 - PhiM)
  odds.PhiF <- PhiF/(1 - PhiF)
  
  # Odds Ratio and Product of Recapture Rates
  OP.Phi <- odds.PhiF*odds.PhiM
  OR.Phi <- odds.PhiF/odds.PhiM
  
  # Recapture Rates (Rho)
  gu <-  min(sqrt(OR.Phi), 1/sqrt(OR.Phi))
  gl <- -min(sqrt(OP.Phi), 1/sqrt(OP.Phi))
  
  ##Correlation using four parameter beta (with FH bounds)
  raw_gamma <- -gl/(gu-gl) # rbeta(1,1,1)
  gamma <- (gu - gl)*raw_gamma + gl
  raw_gamma <- (gamma+1)/2
  
  ###Binomial SD for survival
  sig.PhiF <- sqrt(PhiF*(1-PhiF))
  sig.PhiM <- sqrt(PhiM*(1-PhiM))
  
  #Joint Survival probabilities for paired individuals
  Phifm <- gamma*sig.PhiF*sig.PhiM + PhiF*PhiM
  Phi00 <- 1 - PhiF - PhiM + Phifm
  Phif0 <- PhiF - Phifm
  Phim0 <- PhiM - Phifm
  
  # Conditional Probabilities from Female Perspective
  Phif_M1 <- Phifm/(PhiM)
  Phif_M0 <- Phif0/(1-PhiM)
  
  # Partnership params  --------------------------------------------------------------
  delta <- 1 #rbeta(1,1,1)
  beta <- rbeta(1,1,1)
  
  # Missing Data Simulation --------------------------------------------------------
  
  # Intermediate objects defined within NIMBLE
  single_female  <- matrix(NA, nrow = nf, ncol = k)
  psi_cond       <- array(NA, dim = c(nf, nm+1, k))
  male_taken_jt  <- matrix(NA, nrow = (nm+1), ncol = k)
  male_taken_jt[nm+1, 1:k] <- 0
  prob_repartner_it <- matrix(NA, nrow = nf, ncol = k)
  
  
  # Randomly Sample outcomes using initial values
  for(t in 1:k){
    
    # 1a. Decision to Mate ---------------------------------------------------------
    
    # Female Mating
    for(i in 1:nf){
      if(t >= first_capture_f[i]){
        amating_f[i,t] <- ifelse(is.na(amating_f[i,t]),
                                 rbinom(1, 1, delta * af[i,t]), 
                                 amating_f[i,t])
      }
    }
    
    # Male Mating
    for(j in 1:nm){
      if(t >= first_capture_m[j]){
        amating_m[j,t] <- ifelse(is.na(amating_m[j,t]),
                                 rbinom(1, 1, delta * am[j,t]), 
                                 amating_m[j,t])
      }
    }
    
    
    # 1b. Repartnership at time t-------------------------------------------------------------------------------------------------
    for(i in 1:nf){
      if(t <= first_capture_f[i]){
        prob_repartner_it[i,t] <- 0
        arepartner[i,t] <- 0
      }
      
      if(t >= (first_capture_f[i]+1)){
        prob_repartner_it[i,t] <- compute_prob_repartner_it(beta,
                                                            amating_f[i,t],
                                                            amating_m[apairs_f[i,t-1],t],
                                                            apairs_f[i,t-1],
                                                            na_repartner[i,t],
                                                            psi[i,1:(nm+1),t],
                                                            nm)
        
        arepartner[i,t] <- ifelse(is.na(arepartner[i,t]), 
                                  rbinom(1, 1, prob_repartner_it[i,t]),
                                  arepartner[i,t])
      }
    }
    
    for(j in 1:nm){
      if(t <= first_capture_m[j]){
        male_taken_jt[j, t] <- 0
      } else {
        male_taken_jt[j,t] <- sum(vectorMatch(apairs_f[1:nf,t-1],j)*arepartner[1:nf,t])
      }
    }
    
    # 1b. Decision to Mate at t-----------------------------------------------------------------------------------------------------------------
    for(i in 1:nf){
      
      if(t==first_capture_f[i]){

        apairs_f[i,first_capture_f[i]] <- ifelse(is.na(apairs_f[i,first_capture_f[i]]),
                                                 compute_psi_cond_it(psi[i,1:(nm+1),first_capture_f[i]],
                                                                     amating_f[i,first_capture_f[i]],
                                                                     amating_m[1:nm,first_capture_f[i]],
                                                                     0,
                                                                     male_taken_jt[1:(nm+1),first_capture_f[i]],
                                                                     (nm+1),
                                                                     next_partner_matrix[i,first_capture_f[i]],
                                                                     current_partner_matrix[i, first_capture_f[i]],
                                                                     nf,
                                                                     nm),
                                                 apairs_f[i,first_capture_f[i]])
        single_female[i,first_capture_f[i]] <- 1 * (apairs_f[i,first_capture_f[i]] == (nm + 1))
        
      } else if(t > first_capture_f[i]){
        
        apairs_f[i,t] <- ifelse(is.na(apairs_f[i,t]),
                                compute_psi_cond_it(psi[i,1:(nm+1),t],
                                                    amating_f[i,t],
                                                    amating_m[1:nm,t],
                                                    arepartner[i,t],
                                                    male_taken_jt[1:(nm+1),t],
                                                    apairs_f[i,t-1],
                                                    next_partner_matrix[i,t],
                                                    current_partner_matrix[i, t],
                                                    nf,
                                                    nm),
                                apairs_f[i,t])
        single_female[i,t] <- 1 * (apairs_f[i,t] == (nm + 1))
        
      }
    }
    
    
    # Survival Outcome ----------------------------------------------------------
    if(t < k){
      # Marginal Survival Event for Males in the Population (Y^M_T)---------------------------------------------
      for(j in 1:nm){
        if(t >= first_capture_m[j]){
          am[j,t+1] <- ifelse(is.na(am[j,t+1]),
                              rbinom(1,1, PhiM * am[j,t]),
                              am[j,t+1])
        }
      }
      
      # Marginal Recapture Event for Females in the Population ([X^F_T|X^M_T])
      for(i in 1:nf){
        if(t >= first_capture_f[i]){
          af[i, t+1] <- ifelse(is.na(af[i,t+1]),
                               rbinom(1,1, (single_female[i,t] * PhiF + (1-single_female[i,t]) * (am[apairs_f[i,t],t+1] * Phif_M1 + (1-am[apairs_f[i,t],t+1]) * Phif_M0)) * af[i,t]),
                               af[i,t+1])
        }
      }
    }
  }
  
  # Update Initial Values to follow NIMBLE structure -----------------------------------------------------------------
  
  # Fn to Replace known values with NA and NA values with initial values
  build_NA_mat <- function(mat, ps_mat){
    mat_final <- matrix(NA,nrow = dim(mat)[1], ncol = dim(mat)[2])
    mat_final[is.na(ps_mat)] <- mat[is.na(ps_mat)]
    return(mat_final)
  }
  
  build_NA_vec <- function(vec, ps_vec){
    vec_final <- rep(NA, length(ps_vec))
    vec_final[is.na(ps_vec)] <- vec[is.na(ps_vec)]
    return(vec_final)
  }
  
  
  #Female Mating
  amating_f <- build_NA_mat(amating_f, ps_data$amating_f)
  
  # Male Mating
  amating_m <- build_NA_mat(amating_m, ps_data$amating_m)
  
  # Repartnership
  arepartner <- build_NA_mat(arepartner, ps_data$arepartner)
  
  # Partnership
  # apairs_f <- build_NA_mat(apairs_f, ps_data$apairs_f)
  
  # Female Survival
  af <- build_NA_mat(af, ps_data$af)
  
  # Male Survival
  am <- build_NA_mat(am, ps_data$am)
  
  
  
  # Return Results ------------------------------------------------------------------
  
  # Store in object
  ps_inits <- list(
    PF                = PF,
    PM                = PM,
    raw_rho           = raw_rho,
    rho               = rho,
    PhiF              = PhiF,
    PhiM              = PhiM,
    raw_gamma         = raw_gamma,
    gamma             = gamma,
    # delta             = delta,
    beta              = beta,
    # amating_m         = amating_m,
    # amating_f         = amating_f,
    arepartner        = arepartner,
    af                = af,
    am                = am,
    apairs_f          = apairs_f,
    psi_cond          = psi_cond,
    single_female     = single_female,
    male_taken_jt     = male_taken_jt,
    prob_repartner_it =  prob_repartner_it,
    current_partner_matrix = current_partner_matrix 
  )
  
  # Return Initial Values for a single chain
  return(ps_inits)
}

# Compile Model
compile_pair_swap_nimble <- function(ps_data,
                                     params = NULL){
  
  # Generating Initial Values
  cat("Generating Initial Values...", "\n")
  nimble_inits <- generate_nimble_init_pairs(ps_data)
  
  # Construct Nimble Objects 
  cat("Organizing Data for Nimble...", "\n")
  nimble_ps_constants <- list(
    nf              = ps_data$nf,
    nm              = ps_data$nm,
    k               = ps_data$k,
    first_capture_f = ps_data$first_capture_f,
    first_capture_m = ps_data$first_capture_m,
    psi             = ps_data$psi,
    na_repartner    = ps_data$na_repartner,
    next_partner_matrix = ps_data$next_partner_matrix,
    current_partner_matrix = nimble_inits$current_partner_matrix
  )
  
  nimble_ps_dat <- list(
    # amating_f       = ps_data$amating_f,
    # amating_m       = ps_data$amating_m,
    arepartner      = ps_data$arepartner,
    af              = ps_data$af,
    am              = ps_data$am,
    # apairs_f        = ps_data$apairs_f,
    recap_f         = ps_data$recap_f,
    recap_m         = ps_data$recap_m,
    constraint_data = c(1,1)
  )
  
  if(!is.null(params)){
    cat("User-specified Params Detected...","\n")
    cat("Using params := ", "\n")
    cat(params, "\n")
    
    nimble_params <- params
  } else {
    nimble_params <- c("PF","PM","PhiF","PhiM",
                       "gl","gu","gamma", "beta",
                       "ru","rl","rho")
    cat("Params argument is NULL...","\n")
    cat("Using params := ", "\n")
    cat(nimble_params, "\n")
  }
  
  nimble_dims <- list(
    single_female    = c(nimble_ps_constants$nf, nimble_ps_constants$k),
    psi_prob         = c(nimble_ps_constants$nf, nimble_ps_constants$nm+1),
    psi_cond         = c(nimble_ps_constants$nf, nimble_ps_constants$nm+1, nimble_ps_constants$k),
    male_taken_jt    = c(nimble_ps_constants$nm+1, nimble_ps_constants$k),
    prob_repartner_it = c(nimble_ps_constants$nf, nimble_ps_constants$k),
    apairs_f          = c(nimble_ps_constants$nf, nimble_ps_constants$k)
  )
  
  cat("Building Model Nodes in Nimble (SLOW)...", "\n")
  psModel <- nimbleModel(code       = nimble_ps_model, 
                         constants  = nimble_ps_constants, 
                         inits      = nimble_inits,
                         data       = nimble_ps_dat,
                         dimensions = nimble_dims)
  
  
  # psModel$simulate()
  lp_init <- psModel$calculate()
  print(paste0("LP from initial values is ", round(lp_init,3)))
  
  cat("Compiling Graphical Model in C++ (SLOW)...", "\n")
  compile_ps <- compileNimble(psModel, showCompilerOutput = F)
  
  # [Note] SafeDepare.... warnings are annoying so suppress messages 
  cat("Configuring Markov Chain Monte Carlo Process (SLOW)...", "\n")
  psConf  <- suppressMessages(configureMCMC(model = psModel,
                                            print = F,
                                            multivariateNodesAsScalars = T, 
                                            monitors = nimble_params,
                                            onlySlice = F,
                                            useConjugacy = F))
  
  # Display Samplers
  print(psConf)
  
  cat("Adding Monitors and Constructing MCMC...", "\n")
  psConf$addMonitors(nimble_params)
  psMCMC  <- buildMCMC(psConf)
  
  cat("Compiling MCMC Samplers (SLOW)...", "\n")
  CmdlMCMC <- compileNimble(psMCMC, project = psModel)
  
  cat("Project Defined, MCMC and Model are compiled...", "\n")
  
  cat("Returning Model Object...", "\n")
  
  return(list(CmdlMCMC      = CmdlMCMC, 
              nimble_inits = nimble_inits))
}


# Get Samples from Model
run_nimble <- function(CmdlMCMC, 
                       niter,
                       nburnin,
                       nthin,
                       # inits,
                       nchains=3,
                       seed = F){
  
  cat("MCMC Sampling from Model...","\n")
  samples <- runMCMC(mcmc              = CmdlMCMC,
                     niter             = niter,
                     nburnin           = nburnin, 
                     thin              = nthin,
                     # inits             = inits,
                     nchains           = nchains,
                     setSeed           = seed,
                     samplesAsCodaMCMC = TRUE)
  
  cat("Returning Output...","\n")
  return(samples)
}

# Run compilation/initialization/sampling in one call 
# Dummy function intended for parallelization but can be used in serial as well
execute_pair_swap_nimble_pipeline <- function(seed, 
                                              data, 
                                              params, 
                                              niter, 
                                              nthin,
                                              nburnin, 
                                              nchains){
  
  nimble_complied <- compile_pair_swap_nimble(data, params)
  # nim_inits <- lapply(1:nchains, function(i) generate_nimble_init_pairs(data))
  
  # inits <- generate_nimble_init_pairs(data)
  samples <- run_nimble(CmdlMCMC = nimble_complied$CmdlMCMC,
                        niter    = niter,
                        nthin     = nthin,
                        nburnin  = nburnin,
                        nchains  = nchains,
                        # inits    = nim_inits,
                        seed     = seed) 
  return(list(samples = samples,
              inits   = nimble_complied$nimble_inits))
  
}

run_pair_swap_nimble_parallel <- function(data, params, niter, nthin, nburnin, ncores){
  
  # Build cluster
  ncores <- min(ncores, detectCores()-2)
  
  if(ncores == 1){
    cat("Only one core specified, not bothering with parallelization. Set ncores > 1 if so desired.")
    samples <- execute_pair_swap_nimble_pipeline(F, data, params, niter, nthin, nburnin, 1)
    return(samples)
  }
  
  cat(paste0("Building cluster with ",ncores , " sockets...", "\n"))
  cl <- makeCluster(ncores)
  
  # Load packages necessary
  cat(paste0("Loading custom functions onto cluster...", "\n"))
  clusterEvalQ(cl, {
    #  Load Libraries
    libs <- c("boot", "ggplot2", "nimble", "coda", "ggmcmc", "tidyverse")
    lapply(libs,require, character.only = T)
    source(paste0(getwd(), "/Scripts/pair_swap_mod_nimble7.R"))
    `%+%` <- function(a, b) paste0(a, b)
  })
  
  seeds <-  1:ncores #sample(.Machine$integer.max,ncores)
  
  cat(paste0("Running MCMC in Parallel (SLOW) ...", "\n"))
  out_list <- parLapply(cl      = cl,
                        X       = seeds, 
                        fun     = execute_pair_swap_nimble_pipeline, 
                        data    = data,
                        params  = params,
                        niter   = niter,
                        nthin   = nthin,
                        nburnin = nburnin,
                        nchains = 1)
  
  cat(paste0("Success, closing cluster ...", "\n"))
  
  # Pass Nimble Functions 
  stopCluster(cl)
  
  cat(paste0("Formatting and Returning output ...", "\n"))
  samples <- as.mcmc.list(lapply(1:nchains, function(x) out_list[[x]]$samples))
  inits   <- lapply(1:nchains, function(x) out_list[[x]]$inits)
  
  return(list(samples = samples,
              inits   = inits,
              seed    = seeds))
  
}