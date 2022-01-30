# Nimble Functions for Running Pair Swap Model ------------------------------------------------------------------------------------------------------------ 

# Check if any entries in a vector are equal to the provided float
anyMatch <- nimbleFunction(
  
  run = function(x= double(1),
                 y = double(0)){
    
    returnType(double(0))
    
    output <- 1*any(y == x)
    
    return(output)}
)

# Produce vector of 1s and 0s to check for matching value in an existing vector
vectorMatch <- nimbleFunction(
  
  run = function(x= double(1),
                 y = double(0)){
    
    returnType(double(1))
    
    output <- 1*(y == x)
    
    return(output)}
)

# Pair-Swap random assignment distribution function
dpaircat <- nimbleFunction(
  run = function(x = double(1), 
                 psi_cond = double(2), 
                 nf = integer(0),
                 nm = integer(0), 
                 log = integer(0, default = 1)){
    
    possible_mates <- rep(0,nf)
    psi_cond2 <- matrix(value = 0, nrow = nf, ncol = nm + 1)
    
    returnType(double(0))
    
    psi_cond2[1, 1:(nm+1)] <- c(psi_cond[1,1:nm],equals(sum(psi_cond[1,1:nm]),0))
    possible_mates[1] <- sum(psi_cond2[1,])
    
    for(i in 2:nf){
      
      # Remove Formed Pairs
      for(j in 1:nm){
        psi_cond2[i,j] <- psi_cond[i,j]*(1-1*any(j == x[1:(i-1)]))
      }
      
      #Add case in which no pairs are available 
      psi_cond2[i,(nm+1)] <- equals(sum(psi_cond2[i,1:nm]),0)
      
      possible_mates[i] <- sum(psi_cond2[i,])
    }
    
    logProb <- sum(log(1.0/possible_mates))
    
    if(log) return(logProb)
    else return(exp(logProb))
  }
)

# Pair-Swap random assignment random value generator function
rpaircat <- nimbleFunction(
  run = function(n = integer(0),psi_cond = double(2), nf = integer(0), nm = integer(0)){
    
    x <- rep(0,nf)
    psi_cond2 <- matrix(value = 0, nrow = nf, ncol = nm + 1)
    
    returnType(double(1))
    if(n != 1) print("rpaircat only allows n = 1; using n = 1.")
    
    psi_cond2[1, 1:(nm+1)] <- c(psi_cond[1,1:nm],equals(sum(psi_cond[1,1:nm]),0))
    x[1] <- rcat(n=1, prob = psi_cond2[1, 1:(nm+1)])
    
    for(i in 2:nf){
      
      # Remove Formed Pairs
      for(j in 1:nm){
        psi_cond2[i,j] <- psi_cond[i,j]*(1-1*any(j == x[1:(i-1)]))
      }
      
      #Add case in which no pairs are available 
      psi_cond2[i,(nm+1)] <- equals(sum(psi_cond2[i,1:nm]),0)
      
      # Find mate for i 
      x[i] <- rcat(n=1, prob = psi_cond2[i, 1:(nm+1)])
      
    }
    
    return(x)
  }
)

# BUGS/JAGS Code
nimble_ps_model <- nimbleCode({
  
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
  
  # Initialize History Array (All Zero at time 1)
  for(i in 1:(nf)){
    for(j in 1:(nm+1)){
      histories[i, j, 1] <- 0
    }
  }
  
  # Conditional Partnership/Survival/Recapture Steps ----------------------------------------------------------------------------------------
  
  # Model Events from t=1 to k --------------------------------------------------------------------------------------------------------------
  for(t in 1:k){
    
    # 2. Decision to Mate -------------------------------------------------------------------------------------------------------------------
    
    # Female Mating Choice at time t
    for(i in 1:nf){
      amating_f[i,t] ~ dbern(af[i,t] * recruit_f[i,t] * delta)
    }
    
    # # Male Mating Choice at time t
    for(j in 1:nm){
      amating_m[j,t] ~ dbern(am[j,t] * recruit_m[j,t] * delta)
    }
    
    # Choose to re-form pairs
    for(i in 1:nf){
      prob_repartner[i,t] <- ilogit(beta0 + beta1*histories[i, apairs_f[i,t] , t]) * psi[i, apairs_f[i,t], t]
      arepartner[i,t] ~ dbern(prob_repartner[i,t] * amating_f[i,t] * amating_m[apairs_f[i,t],t])
    }
    
    # Is Male j taken at time t based on re-partnership? 
    # we need Exclude Males who are now unavailable from the catalog of non-repairing individuals
    for(j in 1:nm){
      male_taken_jt[j,t] <- sum(vectorMatch(apairs_f[1:nf,t],j)*arepartner[1:nf,t])
    }
    
    # 3. Mate Selection -------------------------------------------------------------------------------------------------------------------
    # Use Categorical Distribution to classify mates
    
    # Build Homogeneous Partnership probabilities 
    for(i in 1:nf){
      # Flat likelihood of mating conditional on decision to mate
      # If repairing then force partner to be only choice
      # If not repairing then exclude past partner plus any non-mating males
      for(j in 1:nm){
        psi_cond[i, j, t] <- (psi[i,j,t] * amating_f[i,t] * amating_m[j,t] * (1-equals(apairs_f[i,t],j)) * (1-arepartner[i,t]) * (1-male_taken_jt[j,t]) + 
                                arepartner[i,t] * equals(apairs_f[i,t],j))  
      }
    }
    
    # Assign Pairs using Custom Random Matrix Sampling Distribution 
    apairs_f[1:nf,t+1] ~ dpaircat(psi_cond[1:nf,1:nm,t], nf, nm)
    single_female[1:nf,t] <- vectorMatch(apairs_f[1:nf,t+1], nm + 1) 
    
    # Update Total History for Next Time Step
    for(i in  1:nf){
      for(j in  1:(nm+1)){
        histories[i, j, t+1] <- histories[i, j, t] + equals(apairs_f[i,t+1],j) * (1 - single_female[i,t])
      }
    }
    
    # 4. Joint Survival ---------------------------------------------------------------------------------------------------------------------
    
    # Marginal Survival Event for Males in the Population (P[Y^M_T])
    for(j in 1:nm){
      am[j,t+1] ~ dbern(PhiM * am[j,t] * recruit_m[j,t]  + (1-recruit_m[j,t]))
    }
    # Marginal Recapture Event for Males in the Population (P[X^M_T|X^F_T])
    for(i in 1:nf){
      # Probability of female surviving given partnership and partner recapture status
      phi.totalF[i, t] <- single_female[i,t] * PhiF + # female was single
        (1 - single_female[i,t]) * (am[apairs_f[i,t+1],t+1] * (Phifm/PhiM) + # Male mated and female surived
                                      (1 - am[apairs_f[i,t+1],t+1]) * (Phif0/(1-PhiM))) # Male mated and female perished

      # Draw Survival Event
      af[i, t+1] ~ dbern(phi.totalF[i,t] * af[i,t] * recruit_f[i,t] + (1-recruit_f[i,t]))
    }
    
    # 5. Joint Recapture --------------------------------------------------------------------------------------------------------------------
    
    # Marginal Recapture Event for Males in the Population (P[X^M_T])
    for(j in 1:nm){
      recap_m[j,t] ~ dbern(PM * am[j,t+1] * recruit_m[j,t])
    }
    
    # Marginal Recapture Event for females in the Population (P[X^F_T|X^M_T])
    for(i in 1:nf){
      
      # Probability of female being captured given partnership and partner recapture status
      p.totalF[i, t] <- single_female[i,t] * PF + # Female was single
        (1 - single_female[i,t]) * (recap_m[apairs_f[i,t+1],t] * (Pfm/PM) + # Male mated and female captured
                                      (1 - recap_m[apairs_f[i,t+1],t]) * (Pf0/(1-PM))) # Male mated and female not captured
      
      # Draw Recapture Probability
      recap_f[i, t] ~ dbern(p.totalF[i,t] * af[i,t+1] * recruit_f[i,t])
    }
  }
  
  # 6. Prior Distributions-------------------------------------------------------------------------------------------------------------------
  
  # Recruitment 
  for(t in 1:k){
    eps[t] ~ dbeta(1,1)
  }
  
  # Attempt to Mate 
  delta ~ dbeta(1,1)
  
  # Pairs reforming
  beta0 ~ dnorm(0, 1)
  beta1 ~ dnorm(0, 1)
  
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
  gamma <- (gu - gl)*gamma_raw + gl 
  gamma_raw ~ dbeta(1,1)

  # Bounds for Correlation
  
  # Survival Rates (Gamma)
  gu <-  min(sqrt(OR.Phi), 1/sqrt(OR.Phi)) 
  gl <- -min(sqrt(OP.Phi), 1/sqrt(OP.Phi)) 
  
  # Odds Ratio and Product of Survival Rates
  OP.Phi <- odds.PhiF*odds.phiM
  OR.Phi <- odds.PhiF/odds.phiM
  
  ### Odds of Recapture Rates
  odds.phiM <- PhiM/(1 - PhiM)
  odds.PhiF <- PhiF/(1 - PhiF)
  
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
  rho <- (ru - rl)*rho_raw + rl 
  rho_raw ~ dbeta(1,1)
  
  # Bounds for Correlation
  
  # Recapture Rates (Rho)
  ru <-  min(sqrt(OR.P), 1/sqrt(OR.P)) 
  rl <- -min(sqrt(OP.P), 1/sqrt(OP.P)) 
  
  # Odds Ratio and Product of Recapture Rates
  OP.P <- odds.PF*odds.PM
  OR.P <- odds.PF/odds.PM
  
  ### Odds of Survival and Recapture Rates
  odds.PM <- PM/(1 - PM)
  odds.PF <- PF/(1 - PF)
  
  ### Prior Parameters ####
  
  # Recapture Rates M/F
  PF ~ dbeta(1,1)
  PM ~ dbeta(1,1)
})

# Generating Initial Values
generate_nimble_init_pairs <- function(jags_data){
  
  #Unpack Variables -----------------------------------------------------------------
  # Indexes
  k <- jags_data$k
  nf <- jags_data$nf
  nm <- jags_data$nm
  psi <- jags_data$psi # index who is taken
  
  # CR data with missing components
  recruit_f <- jags_data$recruit_f
  recruit_m <- jags_data$recruit_m
  amating_f <- jags_data$amating_f
  amating_m <- jags_data$amating_m
  arepartner <- jags_data$arepartner
  apairs_f <-  jags_data$apairs_f
  af <- jags_data$af
  am <- jags_data$am
  
  # Define local fn equals to emulate jags code
  # Equals call (1 if T; 0 if F)
  equals <- function(a, b){
    return(1*(a == b))
  }
  
  #Randomly Sample categorical Distribution
  rcat <- function(prob){
    return(which(rmultinom(1,1,prob)==1))
  }
  
  # Recapture Prob and Correlation -------------------------------------------------
  PF <- rbeta(1,2,2)
  PM <- rbeta(1,2,2)
  
  ### Odds of Recapture Rates
  odds.PM <- PM/(1 - PM)
  odds.PF <- PF/(1 - PF)
  
  # Odds Ratio and Product of Recapture Rates
  OP.P <- odds.PF*odds.PM
  OR.P <- odds.PF/odds.PM
  
  # Recapture Rates (Rho)
  ru <-  min(sqrt(OR.P), 1/sqrt(OR.P))
  rl <- -min(sqrt(OP.P), 1/sqrt(OP.P))
  
  ##Correlation using four parameter beta (with FH bounds)
  rho_raw <- rbeta(1,3,3)
  rho <- (ru - rl)*rho_raw + rl
  
  ###Binomial SD for recapture
  sig.PF <- sqrt(PF*(1-PF))
  sig.PM <- sqrt(PM*(1-PM))
  
  #Joint Capture probabilities for paired individuals
  Pfm <- rho*sig.PF*sig.PM + PF*PM
  P00 <- 1 - PF - PM + Pfm
  Pf0 <- PF - Pfm
  Pm0 <- PM - Pfm
  
  # Survival Prob and Correlation -------------------------------------------------
  PhiF <- rbeta(1,2,2)
  PhiM <- rbeta(1,2,2)
  
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
  gamma_raw <-  rbeta(1,3,3)
  gamma <- (gu - gl)*gamma_raw + gl
  
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
  
  # Attempt to Mate
  delta <- rbeta(1, 3, 2)
  
  # Pairs reforming
  beta0 <- rnorm(1, 0, 1/4)
  beta1 <- rnorm(1, 0, 1/4)
  
  # Missing Data Simulation --------------------------------------------------------
  
  # amating f/m
  # arepartner f/m
  # af/am
  # apairs f/m
  
  # package up the mating stuff into a few functions its a little unweildy
  
  # Recruitment
  # Female Recruitment
  for(i in 1:nf){
    recruit_f[i,1] <- ifelse(is.na(recruit_f[i,1]), rbinom(1, 1, eps[1]),recruit_f[i,1])
    for(t in 2:(k-1)){
      recruit_f[i,t] <- ifelse(is.na(recruit_f[i,t]), rbinom(1, 1, (recruit_f[i,t-1] + (1-recruit_f[i,t-1]) * eps[t])),recruit_f[i,t])
    }
  }
  
  # Male Recruitment
  for(j in 1:nm){
    recruit_m[j,1] <- ifelse(is.na(recruit_m[j,1]), rbinom(1, 1, eps[1]), recruit_m[j,1])
    for(t in 2:(k-1)){
      recruit_m[j,t] <- ifelse(is.na(recruit_m[j,t]), rbinom(1, 1, (recruit_m[j,t-1] + (1-recruit_m[j,t-1]) * eps[t])), recruit_m[j,t])
    }
  }
  
  # Intermediate objects defined within JAGS
  histories <- array(0, dim = c(nf, nm+1, k+1))
  single_female <- matrix(NA, nrow = nf, ncol = k)
  prob_repartner <- matrix(NA, nrow = nf, ncol = k)
  phi.totalF <- matrix(NA, nrow = nf, ncol = k)
  male_taken_jt <- matrix(NA, nrow = nm, ncol = k)
  psi_raw <- array(NA, dim = c(nf, nm, k))   
  psi_cond <- psi_raw
  psi_cond2 <- array(NA, dim = c(nf, nm + 1, k))
  
  # Time 2 through k initialization
  for(t in 1:k){
    # Female Mating Choice at time t
    for(i in 1:nf){
      amating_f[i,t] <- ifelse(is.na(amating_f[i,t]), rbinom(1, 1, af[i,t] * recruit_f[i,t] * delta), amating_f[i,t])
    }
    
    # Male Mating Choice at time t
    for(j in 1:nm){
      amating_m[j,t] <- ifelse(is.na( amating_m[j,t]),rbinom(1, 1, am[j,t] * recruit_m[j,t] * delta),  amating_m[j,t])
    }
    
    # Choose to re-form pairs
    for(i in 1:nf){
      prob_repartner[i,t] <- inv.logit(beta0 + beta1*histories[i, apairs_f[i,t] , t]) * psi[i, apairs_f[i,t], t]
      arepartner[i,t] <- ifelse(is.na(arepartner[i,t]),
                                rbinom(1,1,prob_repartner[i,t] * amating_f[i,t] * amating_m[apairs_f[i,t],t]),
                                arepartner[i,t])
    }
    
    for(i in 1:nf){
      # Flat likelihood of mating conditional on decision to mate
      # If repairing then force partner to be only choice
      # If not repairing then exclude past partner plus any non-mating males
      for(j in 1:nm){
        psi_raw[i, j, t] <- psi[i,j,t] * amating_f[i,t] * amating_m[j,t] * (1 - equals(apairs_f[i,t],j)) * (1 - arepartner[i,t]) +
          arepartner[i,t] * equals(apairs_f[i,t],j)
      }
    }
    
    #  Exclude Males who are now unavailable
    for(j in 1:nm){
      
      # Is male j available at time t (****based on repartner structure***)
      male_taken_jt[j,t] <- sum(equals(apairs_f[1:nf,t],j)*arepartner[1:nf,t])
      
      # Add Exclusion
      for(i in 1:nf){
        # Remove all possible pairings with females who aren't repairing at t+1 for repairing males
        # (repairing females already have the correct exclusion applied)
        psi_cond[i,j,t] <- psi_raw[i, j, t] * (arepartner[i,t] + (1- arepartner[i,t])*(1-male_taken_jt[j,t]))
      }
    }
    
    # some error handling
    if(!all(sort(which(male_taken_jt[,t] == 1)) == sort(apairs_f[which(arepartner[,t] == 1),t]))) stop("Male taken not working at time:" %+% t)
    
    # Initialize choice selection
    psi_cond2[1, 1:(nm+1), t] <- c(psi_cond[1,1:nm,t],equals(sum(psi_cond[1,1:nm,t]),0))
    # Find mate for 1
    if(is.na(apairs_f[1,t+1])){
      if(is.na(sum(psi_cond2[1, 1:(nm+1), t]))) browser()
      apairs_f[1,t+1] <- rcat(psi_cond2[1, 1:(nm+1), t])
    }
    
    single_female[1,t] <- equals(apairs_f[1,t+1],nm+1)
    
    # Mate Selection
    for(i in 2:nf){
      # Remove Formed Pairs
      for(j in 1:nm){
        
        taken <- (1-sum(equals(apairs_f[1:(i-1),t+1],j)))
        
        if(taken < 0|taken > 1) stop("psi error")
        # something about how psi is indexing is off
        
        psi_cond2[i,j,t] <- psi_cond[i,j,t]*taken
      }
      
      #Add case in which no pairs are available
      psi_cond2[i,(nm+1),t] <- equals(sum(psi_cond2[i,1:nm,t]),0)
      
      # Find mate for i
      if(is.na(apairs_f[i,t+1])){
        if(is.na(sum(psi_cond2[i, 1:(nm+1), t]))) browser()
        apairs_f[i,t+1] <- rcat(psi_cond2[i, 1:(nm+1), t])
      }
      single_female[i,t] <- equals(apairs_f[i,t+1],nm+1)
    }
    
    # some error handling
    if(any(sort(table(apairs_f[,t]))[-length(table(apairs_f[,t]))] > 1)) stop("Illegal apairing_f at " %+% t)
    
    for(i in  1:nf){
      for(j in  1:(nm+1)){
        histories[i, j, t+1] <- histories[i, j, t] + equals(apairs_f[i,t+1],j)*(1-single_female[i,t])
      }
    }
    
    # Marginal Survival Event for Males in the Population (P[Y^M_T])---------------------------------------------
    for(j in 1:nm){
      if(is.na(am[j,t+1])){
        am[j,t+1] <- rbinom(1,1,PhiM * am[j,t] * recruit_m[j,t] + (1-recruit_m[j,t]))
      }
    }
    
    # Marginal Recapture Event for Males in the Population (P[X^M_T|X^F_T])
    for(i in 1:nf){
      
      # Probability of female surviving given partnership and partner recapture status
      phi.totalF[i, t] <- single_female[i,t] * PhiF + # female was single
        (1 - single_female[i,t]) * (am[apairs_f[i,t+1],t+1] * (Phifm/PhiM) + # Male mated and female surived
                                      (1 - am[apairs_f[i,t+1],t+1]) * (Phif0/(1-PhiM))) # Male mated and female perished
      
      
      # print("female:" %+% i %+% "; time:" %+% t %+% "; phi.totalF:" %+% phi.totalF[i,t] %+% "; recruit_f" %+% recruit_f[i,t])
      
      
      # Draw Survival Event
      if(is.na(af[i,t+1])){
        af[i, t+1] <- rbinom(1,1, phi.totalF[i,t] * af[i,t] * recruit_f[i,t] + (1-recruit_f[i,t]))
      }
    }
  }
  
  # Update Initial Values to follow JAGS structure -----------------------------------------------------------------
  
  # Fn to Replace known values with NA and NA values with initial values
  build_NA_mat <- function(mat, jags_mat){
    mat_final <- matrix(NA,nrow = dim(mat)[1], ncol = dim(mat)[2])
    mat_final[is.na(jags_mat)] <- mat[is.na(jags_mat)]
    return(mat_final)
  }
  
  #Female Recruitment
  recruit_f <- build_NA_mat(recruit_f, jags_data$recruit_f)
  
  # Male Recruitment
  recruit_m <- build_NA_mat(recruit_m, jags_data$recruit_m)
  
  # Female Survival
  af <- build_NA_mat(af, jags_data$af)
  
  # Female Mating Status
  amating_f <- build_NA_mat(amating_f, jags_data$amating_f)
  
  # Male Survival
  am <- build_NA_mat(am, jags_data$am)
  
  # Male Mating Status
  amating_m <- build_NA_mat(amating_m, jags_data$amating_m)
  
  # Pair index (female perspective)
  apairs_f <- build_NA_mat(apairs_f, jags_data$apairs_f)
  
  # Repartner index (female perspective)
  arepartner <- build_NA_mat(arepartner, jags_data$arepartner)
  #
  # Return Results ------------------------------------------------------------------
  
  # Store in object
  jags_inits <- list(
    PF = PF,
    PM = PM,
    rho_raw = rho_raw,
    PhiF = PhiF,
    PhiM = PhiM,
    gamma_raw = gamma_raw,
    eps = eps,
    delta = delta,
    beta0 = beta0,
    beta1 = beta1,
    recruit_m = recruit_m,
    recruit_f = recruit_f,
    amating_f = amating_f,
    amating_m = amating_m,
    arepartner = arepartner,
    apairs_f =  apairs_f,
    af = af,
    am = am,
    histories = histories,
    single_female = single_female,
    phi.totalF = phi.totalF,
    male_taken_jt = male_taken_jt,
    prob_repartner = prob_repartner,
    psi_cond = psi_cond
  )
  
  # Return Initial Values for a single chain
  return(jags_inits)
}

# Compile Model
compile_pair_swap_nimble <- function(jags_data,
                                     params = NULL){
  
  # Registering Random Pair-Swap Distribution
  cat("Registering Random Pair-Swap Distribution...", "\n")
  
  registerDistributions(list(
      dpaircat = list(
        BUGSdist = "dpaircat(psi_cond, nf, nm)",
        Rdist = "dpaircat(psi_cond, nf, nm)",
        discrete = T,
        types = c('value = double(1)', 'psi_cond = double(2)', 'nf = integer(0)', 'nm = integer(0)'),
        pqAvail = FALSE)))
  
  # Generating Initial Values
  cat("Generating Initial Values...", "\n")
  nimble_inits <- generate_nimble_init_pairs(jags_data)
  
  # Construct Nimble Objects 
  cat("Organizing Data for Nimble...", "\n")
  nimble_ps_constants <- list(
    nf = jags_data$nf,
    nm = jags_data$nm,
    k = jags_data$k
  )
  
  nimble_ps_dat <- list(
    recruit_f = jags_data$recruit_f,
    recruit_m = jags_data$recruit_m,
    amating_f = jags_data$amating_f,
    amating_m = jags_data$amating_m,
    psi = jags_data$psi,
    af = jags_data$af,
    am = jags_data$am,
    apairs_f = jags_data$apairs_f,
    arepartner = jags_data$arepartner, 
    recap_f = jags_data$recap_f,
    recap_m = jags_data$recap_m
  )
  
  if(!is.null(params)){
    cat("User-specified Params Detected...","\n")
    cat("Using params := " %+% params)
    nimble_params <- params
  } else {
    cat("Params argument is NULL...","\n")
    nimble_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps", "gl", "gu", "ru", "rl")
    cat("Using params := " %+% nimble_params)
  }
  
  
  nimble_dims <- list(histories       = c(nimble_ps_constants$nf, nimble_ps_constants$nm+1, nimble_ps_constants$k+1),
                      prob_repartner  = c(nimble_ps_constants$nf, nimble_ps_constants$k),
                      male_taken_jt   = c(nimble_ps_constants$nm, nimble_ps_constants$k),
                      psi_cond        = c(nimble_ps_constants$nf, nimble_ps_constants$nm, nimble_ps_constants$k),
                      single_female   = c(nimble_ps_constants$nf, nimble_ps_constants$k),
                      phi.totalF      = c(nimble_ps_constants$nf, nimble_ps_constants$k),
                      p.totalF        = c(nimble_ps_constants$nf, nimble_ps_constants$k)
  )
  
  cat("Building Model Nodes in Nimble (SLOW)...", "\n")
  psModel <- nimbleModel(code       = nimble_ps_model, 
                         constants  = nimble_ps_constants, 
                         inits      = nimble_inits,
                         data       = nimble_ps_dat,
                         dimensions = nimble_dims)
  
  cat("Compiling Graphical Model in C++ (SLOW)...", "\n")
  compile_ps <- compileNimble(psModel, showCompilerOutput = TRUE)
  
  
  cat("Configuring Markov Chain Monte Carlo Process (SLOW)...", "\n")
  psConf  <- configureMCMC(psModel, print = TRUE)
  psConf$addMonitors(nimble_params)
  psMCMC  <- buildMCMC(psConf)
  
  cat("Compiling MCMC Samplers (SLOW)...", "\n")
  CpsMCMC <- compileNimble(psMCMC, project = psModel)
  
  cat("Project Defined, MCMC and Model are compiled...", "\n")
  
  cat("Returning Model Object...", "\n")
  
  return(CpsMCMC)
  
  
}

# Get Samples from Model
run_nimble <- function(CpsMCMC, 
                       niter,
                       nburnin,
                       thin){
  
  cat("MCMC Sampling from Model...","\n")
  samples <- runMCMC(CpsMCMC, niter = niter, nburnin = nburnin, thin = thin, setSeed = TRUE, samplesAsCodaMCMC = TRUE)
  cat("Converting to CODA samples...","\n")
  coda.samples <- as.mcmc(samples)
  cat("Returning Output...","\n")
  return(coda.samples)
}