# Nimble Functions for Running Pair Swap Model ------------------------------------------------------------------------------------------------------------ 

# Produce vector of 1s and 0s to check for matching value in an existing vector
vectorMatch <- nimbleFunction(
  
  run = function(x= double(1),
                 y = double(0)){
    
    returnType(double(1))
    
    output <- 1*(y == x)
    
    return(output)}
)

# # Equal probability of selecting any configuration
# dpaircat <- nimbleFunction(
#   run = function(x = double(1),
#                  available_mates = double(2),
#                  nf = integer(0),
#                  nm = integer(0),
#                  log = integer(0, default = 1)){
# 
#     returnType(double(0))
#     logProb <- 0
#     if(log) return(logProb)
#     else return(exp(logProb))
#   }
# )
# 
# # Sample without replacement 2D
# rpaircat <- nimbleFunction(
#   run = function(n = integer(0),
#                  available_mates = double(2),
#                  nf = integer(0),
#                  nm = integer(0)){
# 
# 
#     returnType(double(1))
#     if(n != 1) print("rpaircat only allows n = 1; using n = 1.")
# 
#     x <- rep(nm+1,nf)
#     available_mates2 <- available_mates
#     n_choices <- sum(available_mates2)
# 
#     while(n_choices != 0){
#       choice_flat <- rcat(n =1 , c(available_mates2[1:(nf+1),1:(nm+1)]))
#       row_choice <-  1 + (choice_flat-1) %% (nf+1)
#       col_choice <- ceiling(choice_flat/(nf+1))
# 
#       if(row_choice == (nf+1)){
#         available_mates2[1:(nf+1),col_choice] <- 0
#         n_choices <- sum(available_mates2)
#       } else if(col_choice == (nm+1)){
#         available_mates2[row_choice,1:(nm+1)] <- 0
#         n_choices <- sum(available_mates2)
#       } else {
#         x[row_choice] <- col_choice
#         available_mates2[row_choice,1:(nm+1)] <- 0
#         available_mates2[1:(nf+1),col_choice] <- 0
#         n_choices <- sum(available_mates2)
#       }
#     }
#     return(x)
#   }
# )

# Probability of survival (or recapture) conditional on partner status
compute_prob_condF <- nimbleFunction(
  run = function(is_single_female = double(1),
                 # amating_m = double(1),
                 current_male_state = double(1),
                 current_pairs_f = double(1),
                 ProbF = double(0),
                 ProbM = double(0),
                 Probfm = double(0),
                 Probf0 = double(0),
                 nf = integer(0),
                 nm = integer(0)){
    
    
    returnType(double(1))
    out <- numeric(nf)
    
    # male_mating <- amating_m[current_pairs_f]
    male_current_alive <- current_male_state[current_pairs_f]
    
    # out <- (is_single_female[1:nf] + (1-is_single_female[1:nf])* (1-male_mating)) * ProbF +
    #   (1-is_single_female) * male_mating * (male_current_alive * (Probfm/ProbM) + (1- male_current_alive) * (Probf0/(1-ProbM)))
    out <- is_single_female[1:nf] * ProbF +
      (1-is_single_female) * (male_current_alive * (Probfm/ProbM) + (1- male_current_alive) * (Probf0/(1-ProbM)))
    
    return(out)
  }
)

# Compute Conditional Partnership Grid
compute_psi_cond <- nimbleFunction(
  run = function(psi_t = double(2),
                 amating_f = double(1),
                 amating_m = double(1),
                 nf = integer(0),
                 nm = integer(0)){
    
    returnType(double(2))
    psi_cond_t <- psi_t
    
    for(i in 1:nf){
      for(j in 1:nm){
        psi_cond_t[i,j] <- psi_t[i,j] *  amating_f[i] * amating_m[j]
      }
    }
    
    return(psi_cond_t)
  }
)

# MH-Step Sampler
sampler_pairs <- nimbleFunction(
  name = 'sampler_pairs',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
    calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    # Vars
    nf <- model$getParam(target, 'nf')
    nm <- model$getParam(target, 'nm')
    psi_cond_t <- model$getParam(target, 'available_mates')
    ## checks
    if(model$getDistribution(target) != 'dpaircat') stop('can only use pair categorical sampler on node with dpaircat distribution')
  },
  run = function() {
    currentLogProb <- model$getLogProb(calcNodes)
    model[[target]] <<- rpaircat(1,psi_cond_t, nf,nm)
    otherLogProbPrior <- model$calculate(target)
    
    if(otherLogProbPrior == -Inf) {
      otherLogProb <- otherLogProbPrior
    } else {
      otherLogProb <- otherLogProbPrior + model$calculate(calcNodesNoSelf)
    }
    
    acceptanceProb <- 1/(exp(currentLogProb - otherLogProb) + 1)
    jump <- (!is.nan(acceptanceProb)) & (runif(1,0,1) < acceptanceProb)
    if(jump) {
      nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
    } else {
      nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
    }
  },
  methods = list(
    reset = function() { }
  )
)


# BUGS/JAGS Code
nimble_ps_model <- nimbleCode({
  
  # 0. Data Augmentation-----------------------------------------------------------------------------------------------------------------
  for(i in 1:nf){
    zf[i] ~ dbern(xi)
  }
  
  for(j in 1:nm){
    zm[j] ~ dbern(xi)
  }
  
  # Number of males and females
  NF <- sum(zf[1:nf])
  NM <- sum(zm[1:nm])
  
  # 1. Recruitment into population-----------------------------------------------------------------------------------------------------------
  # Female Recruitment
  for(i in 1:nf){
    recruit_f[i,1] ~  dbern(eps[1])
    for(t in 2:(k-1)){
      recruit_f[i,t] ~ dbern(recruit_f[i,t-1] + (1-recruit_f[i,t-1]) * eps[t])
    } 
  }
  
  # Male Recruitment
  for(j in 1:nm){
    recruit_m[j,1] ~  dbern(eps[1])
    for(t in 2:(k-1)){
      recruit_m[j,t] ~ dbern(recruit_m[j,t-1] + (1-recruit_m[j,t-1]) * eps[t])
    } 
  }
  
  # Conditional Partnership/Survival Steps -------------------------------------------------------------------------------------------------
  
  # 2a. Mating at time t-----------------------------------------------------------------------------------------------------------
  # Female Recruitment
  for(i in 1:nf){
    for(t in 1:k){
      amating_f[i,t] <- recruit_f[i,t] * af[i,t]  # * zf[i]dbern(delta * recruit_f[i,t] * af[i,t] * zf[i])
    }
  }
  
  # Male Recruitment
  for(j in 1:nm){
    for(t in 1:k){
      amating_m[j,t] <- recruit_m[j,t] * am[j,t]  # * zm[j]dbern(delta * recruit_m[j,t] * am[j,t] * zm[j])
    }
  }
  
  # 2b. Decision to Mate at t-----------------------------------------------------------------------------------------------------------------
  for(t in 1:k){
    # Need to have both been recruited+alive by time t in order to form a pair-bond 
    psi_cond[1:(nf+1),1:(nm+1), t] <- compute_psi_cond(psi[1:(nf+1),1:(nm+1),t],
                                                       amating_f[1:nf,t],
                                                       amating_m[1:nm,t],
                                                       nf,
                                                       nm)
    
    # Assign Pairs using Custom Random Matrix Sampling Distribution 
    apairs_f[1:nf,t] ~ dpaircat(psi_cond[1:(nf+1), 1:(nm+1), t], nf, nm)
    single_female[1:nf,t] <- vectorMatch(apairs_f[1:nf,t], nm + 1)
  }
  
  # 3. Joint Survival [t-1,t) ---------------------------------------------------------------------------------------------------------------
  for(t in 2:k){
    
    # Marginal Survival Event for Males in the Population (P[Y^M_T]) 
    for(j in 1:nm){
      am[j,t] ~ dbern(PhiM * am[j,t-1] * recruit_m[j,t-1]  + (1-recruit_m[j,t-1]))
    }
    
    # Marginal Recapture Event for Females in the Population (P[X^F_T|X^M_T]) given males
    phi.totalF[1:nf,t-1] <- compute_prob_condF(single_female[1:nf,t-1],
                                               #amating_m[1:(nm+1),t-1],
                                               am[1:(nm+1),t],
                                               apairs_f[1:nf,t-1],
                                               PhiF,
                                               PhiM,
                                               Phifm,
                                               Phif0,
                                               nf,
                                               nm)
    
    # Draw conditional Survival Event
    for(i in 1:nf){
      af[i, t] ~ dbern(phi.totalF[i,t-1] * af[i,t-1] * recruit_f[i,t-1]  + (1-recruit_f[i,t-1]))
    }
  }
  
  # 4. Joint Recapture --------------------------------------------------------------------------------------------------------------------
  for(t in 1:k){
    
    # Marginal Recapture Event for Males in the Population (P[X^M_T])
    for(j in 1:nm){
      recap_m[j,t] ~ dbern(PM * am[j,t] * recruit_m[j,t] * zm[j])
    }
    
    # Marginal Recapture Event for females in the Population (P[X^F_T|X^M_T])
    p.totalF[1:nf,t] <- compute_prob_condF(single_female[1:nf,t],
                                           #amating_m[1:(nm+1),t],
                                           recap_m[1:(nm+1),t],
                                           apairs_f[1:nf,t],
                                           PF,
                                           PM,
                                           Pfm,
                                           Pf0,
                                           nf,
                                           nm)
    
    # Draw Recapture Probability
    for(i in 1:nf){
      recap_f[i, t] ~ dbern(p.totalF[i,t] * af[i,t] * recruit_f[i,t] * zf[i])
    }
  }
  
  # 5. Prior Distributions-------------------------------------------------------------------------------------------------------------------
  # Data augmentation
  xi ~ dbeta(0.1,1.0)
  
  # Mating or not
  # delta ~ dbeta(1,1)
  
  # Recruitment 
  for(t in 1:k){
    eps[t] ~ dbeta(1,1)
  }
  
  # Survival Terms
  
  ### Derived Parameters ####
  
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
  
  #Unpack Variables -----------------------------------------------------------------
  # Indexes
  k <- ps_data$k
  nf <- ps_data$nf
  nm <- ps_data$nm
  psi <- ps_data$psi
  
  # CR data with missing components
  zf <- ps_data$zf
  zm <- ps_data$zm
  recruit_f <- ps_data$recruit_f
  recruit_m <- ps_data$recruit_m
  amating_f <- ps_data$amating_f
  amating_m <- ps_data$amating_m
  apairs_f <-  ps_data$apairs_f
  af <- ps_data$af
  am <- ps_data$am
  recap_m <- ps_data$recap_m
  recap_f <- ps_data$recap_f
  single_female <- ps_data$single_female
  
  # Define local fn equals to emulate jags code
  # Equals call (1 if T; 0 if F)
  equals <- function(a, b){
    return(1*(a == b))
  }
  
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
  raw_rho <- rbeta(1,1,1)
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
  raw_gamma <-  rbeta(1,1,1)
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
  
  # Simple Processes --------------------------------------------------------------
  # Recruitment
  eps <- rep(NA, k)
  for(t in 1:k){
    eps[t] <- rbeta(1,1,1)
  }
  
  xi <- rbeta(1,0.3,1)
  
  # Mating or not
  delta <- 1#rbeta(1,1,1)
  
  # Missing Data Simulation --------------------------------------------------------
  
  # Sample augmentation female
  for(i in 1:nf){
    zf[i] <- ifelse(is.na(zf[i]), rbinom(1, 1, xi),zf[i])
  }
  
  # Sample augmentation male
  for(j in 1:nm){
    zm[j] <- ifelse(is.na(zm[j]), rbinom(1, 1, xi),zm[j])
  }
  
  # Recruitment
  # Female Recruitment
  for(i in 1:nf){
    recruit_f[i,1] <- ifelse(is.na(recruit_f[i,1]),rbinom(1, 1, eps[1]),recruit_f[i,1])
    for(t in 2:(k-1)){
      recruit_f[i,t] <- ifelse(is.na(recruit_f[i,t]), 
                               rbinom(1, 1, (recruit_f[i,t-1] + (1-recruit_f[i,t-1]) * eps[t])),
                               recruit_f[i,t])
    }
  }
  
  # Male Recruitment
  for(j in 1:nm){
    recruit_m[j,1] <- ifelse(is.na(recruit_m[j,1]), rbinom(1, 1, eps[1]), recruit_m[j,1])
    for(t in 2:(k-1)){
      recruit_m[j,t] <- ifelse(is.na(recruit_m[j,t]), 
                               rbinom(1, 1, (recruit_m[j,t-1] + (1-recruit_m[j,t-1]) * eps[t])),
                               recruit_m[j,t])
    }
  }
  
  # Intermediate objects defined within NIMBLE
  single_female  <- matrix(NA, nrow = nf, ncol = k)
  phi.totalF     <- matrix(NA, nrow = nf, ncol = k-1)
  p.totalF       <- matrix(NA, nrow = nf, ncol = k)
  psi_cond       <- array(NA, dim = c(nf+1, nm+1, k))
  
  # Time 2 through k initialization
  for(t in 1:k){
    # If past the first occasion then death is possible 
    if(t > 1){
      # Marginal Survival Event for Males in the Population (Y^M_T)---------------------------------------------
      for(j in 1:nm){
        if(is.na(am[j,t])){
          am[j,t] <- rbinom(1,1, PhiM * am[j,t-1] * recruit_m[j,t-1] + (1-recruit_m[j,t-1]))
        }
      }
      
      # Marginal Recapture Probability for Females in the Population (P[X^F_T|X^M_T])
      phi.totalF[1:nf, t-1] <- compute_prob_condF(is_single_female    = single_female[1:nf,t-1],
                                                  #amating_m           = amating_m[1:(nm+1),t-1],
                                                  current_male_state  = am[1:(nm+1), t],
                                                  current_pairs_f     = apairs_f[1:nf, t-1],
                                                  ProbF               = PhiF,
                                                  ProbM               = PhiM,
                                                  Probfm              = Phifm,
                                                  Probf0              = Phif0,
                                                  nf                  = nf,
                                                  nm                  = nm) 
      
      # Marginal Recapture Event for Females in the Population ([X^F_T|X^M_T])
      for(i in 1:nf){
        if(is.na(af[i,t])){
          af[i, t] <- rbinom(1,1, phi.totalF[i,t-1] * af[i,t-1] * recruit_f[i,t-1]  + (1-recruit_f[i,t-1]))
        }
      }
    }
    
    # Partnership ----------------------------------------------------------------------
    # Assign Pairs using Custom Random Matrix Sampling Distribution 
    # Need to have both been recruited+alive by time t in order to form a pair-bond 
    
    # Female Mating
    for(i in 1:nf){
      amating_f[i,t] <- ifelse(is.na(amating_f[i,t]),
                               rbinom(1, 1, delta * recruit_f[i,t] * af[i,t] * zf[i]), 
                               amating_f[i,t])
    }
    
    # Male Mating
    for(j in 1:nm){
      amating_m[j,t] <- ifelse(is.na(amating_m[j,t]),
                               rbinom(1, 1, delta * recruit_m[j,t] * am[j,t]* zm[j]), 
                               amating_m[j,t])
    }
    
    # Need to have both been recruited+alive by time t in order to form a pair-bond 
    psi_cond[1:(nf+1),1:(nm+1), t] <- compute_psi_cond(psi[1:(nf+1),1:(nm+1),t],
                                                       amating_f[1:nf,t],
                                                       amating_m[1:nm,t],
                                                       nf,
                                                       nm)
    
    apairs_f[1:nf,t] <- rpaircat(n = 1, psi_cond[1:(nf+1), 1:(nm+1), t], nf, nm)
    single_female[1:nf,t] <- equals(apairs_f[1:nf,t],nm+1)
    
    # 5. Joint Recapture --------------------------------------------------------------------------------------------------------------------
    
    # Marginal Recapture Probability for Females in the Population (P[X^F_T|X^M_T])
    p.totalF[1:nf, t] <- compute_prob_condF(is_single_female    = single_female[1:nf,t],
                                            #amating_m           = amating_m[1:(nm+1),t],
                                            current_male_state  = recap_m[1:(nm+1), t],
                                            current_pairs_f     = apairs_f[1:nf, t],
                                            ProbF               = PF,
                                            ProbM               = PM,
                                            Probfm              = Pfm,
                                            Probf0              = Pf0,
                                            nf                  = nf,
                                            nm                  = nm) 
    
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
  
  #Female Recruitment
  recruit_f <- build_NA_mat(recruit_f, ps_data$recruit_f)
  
  # Male Recruitment
  recruit_m <- build_NA_mat(recruit_m, ps_data$recruit_m)
  
  #Female Mating
  amating_f <- build_NA_mat(amating_f, ps_data$amating_f)
  
  # Male Mating
  amating_m <- build_NA_mat(amating_m, ps_data$amating_m)
  
  # Female Survival
  af <- build_NA_mat(af, ps_data$af)
  
  # Male Survival
  am <- build_NA_mat(am, ps_data$am)
  
  # Existence blocks
  zf <- build_NA_vec(zf, ps_data$zf)
  zm <- build_NA_vec(zm, ps_data$zm)
  
  # Return Results ------------------------------------------------------------------
  
  # Store in object
  ps_inits <- list(
    PF            = PF,
    PM            = PM,
    raw_rho       = raw_rho,
    rho           = rho,
    PhiF          = PhiF,
    PhiM          = PhiM,
    raw_gamma     = raw_gamma,
    gamma         = gamma,
    eps           = eps,
    delta         = delta,
    xi            = xi,
    zf            = zf,
    zm            = zm,
    recruit_m     = recruit_m,
    recruit_f     = recruit_f,
    # amating_m     = amating_m,
    # amating_f     = amating_f,
    af            = af,
    am            = am,
    apairs_f      = apairs_f,
    psi_cond      = psi_cond,
    single_female = single_female,
    phi.totalF    = phi.totalF,
    p.totalF      = p.totalF
  )
  
  # Return Initial Values for a single chain
  return(ps_inits)
}

# Compile Model
compile_pair_swap_nimble <- function(ps_data,
                                     params = NULL){
  
  # Registering Random Pair-Swap Distribution
  cat("Registering Random Pair-Swap Distribution...", "\n")
  
  registerDistributions(list(
    dpaircat = list(
      BUGSdist = "dpaircat(available_mates, nf, nm)",
      Rdist = "dpaircat(available_mates, nf, nm)",
      discrete = TRUE,
      range = c(1, (ps_data$nm+1)),
      types = c('value = double(1)', 
                'available_mates = double(2)',
                'nf = integer(0)', 
                'nm = integer(0)'),
      pqAvail = FALSE)
  ))
  
  # Generating Initial Values
  cat("Generating Initial Values...", "\n")
  nimble_inits <- generate_nimble_init_pairs(ps_data)
  
  # Construct Nimble Objects 
  cat("Organizing Data for Nimble...", "\n")
  nimble_ps_constants <- list(
    nf           = ps_data$nf,
    nm           = ps_data$nm,
    k            = ps_data$k
  )
  
  nimble_ps_dat <- list(
    zf              = ps_data$zf,
    zm              = ps_data$zm,
    recruit_f       = ps_data$recruit_f,
    recruit_m       = ps_data$recruit_m,
    # amating_f       = ps_data$amating_f,
    # amating_m       = ps_data$amating_m,
    af              = ps_data$af,
    am              = ps_data$am,
    apairs_f        = ps_data$apairs_f,
    recap_f         = ps_data$recap_f,
    recap_m         = ps_data$recap_m,
    psi             = ps_data$psi,
    constraint_data = c(1,1)
  )
  
  if(!is.null(params)){
    cat("User-specified Params Detected...","\n")
    cat("Using params := ", "\n")
    cat(params, "\n")
    
    nimble_params <- params
  } else {
    cat("Params argument is NULL...","\n")
    cat("Using params := ", "\n")
    cat(nimble_params, "\n")
  }
  
  nimble_dims <- list(
    single_female    = c(nimble_ps_constants$nf, nimble_ps_constants$k),
    phi.totalF       = c(nimble_ps_constants$nf, nimble_ps_constants$k-1),
    p.totalF         = c(nimble_ps_constants$nf, nimble_ps_constants$k),
    psi_cond         = c(nimble_ps_constants$nf+1, nimble_ps_constants$nm+1, nimble_ps_constants$k)
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
  
  # Assign dpaircat sampler to the global environment
  assign('sampler_pairs', sampler_pairs, envir = .GlobalEnv)
  
  cat("Compiling Graphical Model in C++ (SLOW)...", "\n")
  compile_ps <- compileNimble(psModel, showCompilerOutput = F)
  
  node_names <- psModel$getNodeNames()[psModel$getNodeType(psModel$getNodeNames())=="stoch"]
  node_names <- node_names[!substr(node_names, 1, nchar("apairs_f")) == "apairs_f"]
  
  # [Note] SafeDepare.... warnings are annoying so suppress messages 
  cat("Configuring Markov Chain Monte Carlo Process (SLOW)...", "\n")
  psConf  <- suppressMessages(configureMCMC(model = psModel,
                                            print = F,
                                            nodes = node_names,
                                            multivariateNodesAsScalars = T, 
                                            monitors = nimble_params,
                                            onlySlice = F,
                                            useConjugacy = T))
  
  cat("Adding custom random pair sampler to apairs_f...", "\n")
  
  # Default sampler wont work here....
  psConf$removeSampler("apairs_f", print = F)
  
  for(t in 1:ps_data$k){
    temp_name_cat <- "apairs_f[1:" %+% ps_data$nf %+% "," %+% t %+% "]"
    psConf$addSampler(target = temp_name_cat, type = "sampler_pairs", print = T)
  }
  
  # NIMBLE thinks recap_m is partially observed. Remove this and posterior predictive sampler on F.
  psConf$removeSampler("recap_m",print = F)
  psConf$removeSampler("recap_f",print = F)
  
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
                       thin,
                       # inits,
                       nchains=3,
                       seed = F){
  
  cat("MCMC Sampling from Model...","\n")
  samples <- runMCMC(mcmc              = CmdlMCMC,
                     niter             = niter,
                     nburnin           = nburnin, 
                     thin              = thin,
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
  
  # inits <- generate_nimble_init_pairs(data)
  samples <- run_nimble(CmdlMCMC = nimble_complied$CmdlMCMC,
                        niter    = niter,
                        thin     = nthin,
                        nburnin  = nburnin,
                        nchains  = nchains,
                        # inits    = nimble_complied$nimble_inits,
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
    source(paste0(getwd(), "/Scripts/pair_swap_mod_nimble4.R"))
    `%+%` <- function(a, b) paste0(a, b)
  })
  
  seeds <-  sample(.Machine$integer.max,ncores)
  
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