# Nimble Functions for Running Pair Swap Model ------------------------------------------------------------------------------------------------------------ 

# Produce vector of 1s and 0s to check for matching value in an existing vector
vectorMatch <- nimbleFunction(
  
  run = function(x= double(1),
                 y = double(0)){
    
    returnType(double(1))
    
    output <- 1*(y == x)
    
    return(output)}
)
# Arrival dpaircat

# Pair-Swap random assignment distribution function
dpaircat <- nimbleFunction(
  run = function(x = double(1), 
                 available_mates = double(2), 
                 mating_f = double(1),
                 mating_m = double(1),
                 nf = integer(0),
                 nm = integer(0), 
                 log = integer(0, default = 1)){
    
    
    available_mates2 <- matrix( 0, nrow = nf+1, ncol = nm + 1)
    
    returnType(double(0))
    
    # How many available males and females
    mating_females <- sum(mating_f) + 1
    mating_males  <- sum(mating_m) + 1
    osr <- mating_females/mating_males
    
    if(osr <= 1){
      ll <- rep(0, nf)
      # possible_mates <- rep(0,nf)
      available_mates2[1, 1:(nm+1)] <- c(available_mates[1,1:nm],equals(sum(available_mates[1,1:nm]),0))
      ll[1] <- dcat(x[1],available_mates2[1,1:(nm+1)],log = TRUE)
      
      for(i in 2:nf){
        
        # Remove Formed Pairs
        for(j in 1:nm){
          available_mates2[i,j] <- available_mates[i,j]*(1-1*any(j == x[1:(i-1)]))
        }
        
        #Add case in which no pairs are available 
        available_mates2[i,(nm+1)] <- equals(sum(available_mates2[i,1:nm]),0)
        
        ll[i] <- dcat(x[i],available_mates2[i,1:(nm+1)],log = TRUE)
      }
    } else {
      y <- rep(0, nm)
      ll <- rep(0,nm)
      
      # Build vector of male pairs
      # Need base likelihood on these values
      for(j in 1:nm){
        # Find which is your mate
        for(i in 1:nf){
          if(x[i] == j){
            y[j] <- i
          } 
        }
        
        # If no mate set to single state
        if(y[j]==0){
          y[j] <- nf+1
        }
        
      }
      
      
      # Apply mating selection for first male
      available_mates2[1:(nf+1), 1] <- c(available_mates[1:nf,1],equals(sum(available_mates[1:nf,1]),0))
      ll[1] <- dcat(y[1],available_mates2[1:(nf+1),1],log = TRUE)
      
      # Iterate over remaining males
      for(j in 2:nm){
        
        # Remove Formed Pairs
        for(i in 1:nf){
          available_mates2[i,j] <- available_mates[i,j]*(1-1*any(i == y[1:(j-1)]))
        }
        
        #Add case in which no pairs are available 
        available_mates2[(nf+1),j] <- equals(sum(available_mates2[1:nf,j]),0)
        
        # Find mate for i 
        ll[j] <- dcat(y[j],available_mates2[1:(nf+1),j],log = TRUE)
        
      }
    }
    
    logProb <- sum(ll)
    
    if(log) return(logProb)
    else return(exp(logProb))
  }
)


# Pair-Swap random assignment random value generator function
rpaircat <- nimbleFunction(
  run = function(n = integer(0),
                 available_mates = double(2),
                 mating_f = double(1),
                 mating_m = double(1),
                 nf = integer(0),
                 nm = integer(0)){
    
    x <- rep(0,nf)
    available_mates2 <- matrix(value = 0, nrow = nf+1, ncol = nm + 1)
    
    returnType(double(1))
    if(n != 1) print("rpaircat only allows n = 1; using n = 1.")
    
    # How many available males and females
    mating_females <- sum(mating_f) + 1
    mating_males  <- sum(mating_m) + 1
    osr <- mating_females/mating_males
    
    # If operating sex ratio is balanced or in favor of females 
    # Iterate across female
    if(osr <= 1){
      
      available_mates2[1, 1:(nm+1)] <- c(available_mates[1,1:nm],equals(sum(available_mates[1,1:nm]),0))
      x[1] <- rcat(n=1, prob = available_mates2[1, 1:(nm+1)])
      
      for(i in 2:nf){
        
        # Remove Formed Pairs
        for(j in 1:nm){
          available_mates2[i,j] <- available_mates[i,j]*(1-1*any(j == x[1:(i-1)]))
        }
        
        #Add case in which no pairs are available 
        available_mates2[i,(nm+1)] <- equals(sum(available_mates2[i,1:nm]),0)
        
        # Find mate for i 
        x[i] <- rcat(n=1, prob = available_mates2[i, 1:(nm+1)])
        
      }
      
    } else {
      # Otherwise OSR favors males so we iterate across their mating pop instead
      y <- rep(0,nm)
      
      # Apply mating selection for first male
      available_mates2[1:(nf+1), 1] <- c(available_mates[1:nf,1],equals(sum(available_mates[1:nf,1]),0))
      y[1] <- rcat(n = 1, prob = available_mates2[1:(nf+1),1])
      
      # Iterate over remaining males
      for(j in 2:nm){
        
        # Remove Formed Pairs
        for(i in 1:nf){
          available_mates2[i,j] <- available_mates[i,j]*(1-1*any(i == y[1:(j-1)]))
        }
        
        #Add case in which no pairs are available 
        available_mates2[(nf+1),j] <- equals(sum(available_mates2[1:nf,j]),0)
        
        # Find mate for i 
        y[j] <- rcat(n= 1, prob = available_mates2[1:(nf+1), j])
        
      }
      
      # Mapping male choices back to the female vector
      for(i in 1:nf){
        for(j in 1:nm){
          if(y[j] == i){
            x[i] <- j
          } 
        }
        
        #If not mates set to single
        if(x[i]==0){
          x[i] <- nm+1
        }
      }
      
    }
    
    return(x)
  }
)

# rng vectorized bernoulli trials
rmvbern <- nimbleFunction(
  run = function(n = integer(0, default = 1),
                 prob = double(1)){
    
    
    returnType(double(1))
    if(n != 1) print("rmvbern only allows n = 1; using n = 1.")
    
    return(rbinom(length(prob), size = 1, prob = prob))
  }
)

# distn vectorized bernoulli trials
dmvbern <- nimbleFunction(
  run = function(x = double(1),
                 prob = double(1),
                 log  = integer(0, default = 1)){
    
    
    returnType(double(0))
    
    logProb <- sum(dbinom(x, size = 1, prob = prob, log = TRUE))
    
    if(log) return(logProb)
    else return(exp(logProb))
  }
)
# 
compute_pr_repartner <- nimbleFunction(
  run = function(intercept = double(0),
                 slope = double(0),
                 history = double(2),
                 psi_uncond = double(2),
                 mating_f = double(1),
                 mating_m = double(1),
                 former_pairs_f  = double(1),
                 na_repartner = double(1),
                 nf = integer(0),
                 nm = integer(0)){
    
    returnType(double(1))
    
    out <- rep(0,nf)
    
    for(i in 1:nf){
      
      if(former_pairs_f[i]>=(nm+1)){
        out[i] <- 0
      } else if(sum(psi_uncond[i,1:nm]) == 1.0 & psi_uncond[i,former_pairs_f[i]] == 1.0 & na_repartner[i] == 1.0){
        # if only 1 partner available, force nibmle to make sure they pair
        out[i] <- 1
      } else {
        out[i] <- ilogit(intercept + slope * history[i,former_pairs_f[i]]) *
          psi_uncond[i,former_pairs_f[i]] * mating_f[i] * mating_m[former_pairs_f[i]]
      }
    }
    
    return(out)
  }
)

compute_prob_condF <- nimbleFunction(
  run = function(is_single_female = double(1),
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
    
    for(i in 1:nf){
      if(is_single_female[i]==1|current_pairs_f[i]>=(nm+1)|ProbM==0|ProbM == 1){
        out[i] <- ProbF
      } else {
        out[i] <- (current_male_state[current_pairs_f[i]] * (Probfm/ProbM) + # Male mated and female surived
                     (1 - current_male_state[current_pairs_f[i]]) * (Probf0/(1-ProbM)))
      }
    }
    
    return(out)
  }
)

# Custom Sampler for dpaircat 
# Variation on categorical sampler and MH Proposal 
# sampler_pairs <- nimbleFunction(
#   name = 'sampler_pairs',
#   contains = sampler_BASE,
#   setup = function(model, mvSaved, target, control) {
#     ## node list generation
#     targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
#     calcNodes <- model$getDependencies(target)
#     calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
#     isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
#     calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
#     calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
#     # Vars
#     nf <- model$getParam(target, 'nf')
#     nm <- model$getParam(target, 'nm')
#     amating_m <- model$getParam(target, 'mating_m')
#     amating_f <- model$getParam(target, 'mating_f')
#     psi_cond_t <- model$getParam(target, 'available_mates')
#     ## checks
#     if(model$getDistribution(target) != 'dpaircat') stop('can only use pair categorical sampler on node with dpaircat distribution')
#   },
#   run = function() {
#     # Simulate new partners
#     model[[target]] <<- rpaircat(1,psi_cond_t, amating_f, amating_m, nf, nm) # accept target 
#     model$calculate(calcNodes) # calculate logprobs 
#     nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
#     nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
#     nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
#   },
#   methods = list(
#     reset = function() { }
#   )
# )

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
    amating_m <- model$getParam(target, 'mating_m')
    amating_f <- model$getParam(target, 'mating_f')
    psi_cond_t <- model$getParam(target, 'available_mates')
    logProbs <- numeric(2)
    ## checks
    if(model$getDistribution(target) != 'dpaircat') stop('can only use pair categorical sampler on node with dpaircat distribution')
  },
  run = function() {
    current_pairs <- model[[target]]
    logProbs[1] <<- model$getLogProb(calcNodes)
    if(is.nan(logProbs[1])) logProbs[1] <<- -Inf
    # Simulate new partners
    model[[target]] <<- rpaircat(1,psi_cond_t, amating_f, amating_m, nf, nm) # accept target 
    logProbs[2] <<- model$calculate(calcNodes) # calculate logprobs 
    if(is.nan(logProbs[2])) logProbs[2] <<- -Inf
    
    logProbs <<- logProbs - max(logProbs)
    probs <<- exp(logProbs)
    jump <- rcat(1, probs)
    
    if(jump == 2) {
      model$calculate(calcNodes)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
    } else {
      model[[target]] <<- current_pairs
      model$calculate(calcNodes)
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
  
  # # Initialize History Array (All Zero at time 1)
  # for(i in 1:(nf)){
  #   for(j in 1:(nm+1)){
  #     histories[i, j, 1] <- 0
  #   }
  # }
  
  # Conditional Partnership/Survival Steps ----------------------------------------------------------------------------------------
  
  # Model Events from t=1 to k --------------------------------------------------------------------------------------------------------------
  
  # Time 1 (all animals must be alive/no animals are re-partnering yet)
  
  # If they've been recruited do they mate at this time step
  for(i in 1:nf){
    amating_f[i,1] ~ dbern(recruit_f[i,1] * delta *  zf[i])
  }
  
  for(j in 1:nm){
    amating_m[j,1] ~ dbern(recruit_m[j,1] * delta * zm[j])
  }
  
  # Build Homogeneous Partnership probabilities 
  for(i in 1:nf){
    # Flat likelihood of mating conditional on decision to mate
    for(j in 1:nm){
      psi_cond[i, j, 1] <- psi[i,j,1] * amating_f[i,1] * amating_m[j,1] 
    }
  }
  
  # Assign Pairs using Custom Random Matrix Sampling Distribution 
  apairs_f[1:nf,1] ~ dpaircat(psi_cond[1:nf, 1:nm, 1], amating_f[1:nf,1], amating_m[1:nm,1], nf, nm)
  single_female[1:nf,1] <- vectorMatch(apairs_f[1:nf,1], nm + 1) 
  
  # # Update Total History (at start of time 2, what is the history data)
  # for(i in  1:nf){
  #   for(j in  1:(nm+1)){
  #     histories[i, j, 2] <- histories[i, j, 1] + equals(apairs_f[i,1],j) * (1 - single_female[i,1])
  #   }
  # }
  
  # Time 2-k
  for(t in 2:k){
    
    # 2. Decision to Mate -------------------------------------------------------------------------------------------------------------------
    
    # Female Mating Choice at time t
    for(i in 1:nf){
      amating_f[i,t] ~ dbern(af[i,t-1] * recruit_f[i,t] * delta* zf[i])
    }
    
    # Male Mating Choice at time t
    for(j in 1:nm){
      amating_m[j,t] ~ dbern(am[j,t-1] * recruit_m[j,t] * delta* zm[j])
    }
    
    # prob_repartner[1:nf,t-1] <- compute_pr_repartner(beta0,
    #                                                  beta1,
    #                                                  histories[1:nf,1:(nm+1),t],
    #                                                  psi[1:nf,1:(nm+1),t],
    #                                                  amating_f[1:nf,t],
    #                                                  amating_m[1:nm,t],
    #                                                  apairs_f[1:nf,t-1],
    #                                                  na_repartner[1:nf,t-1],
    #                                                  nf,
    #                                                  nm)
    # 
    # # Choose to re-form pairs
    # for(i in 1:nf){
    #   arepartner[i,t-1] ~ dbern(prob_repartner[i,t-1])
    # }
    # 
    # # arepartner[1:nf,t-1] ~ dmvbern(prob_repartner[1:nf,t-1])
    # 
    # # Is Male j from t-1taken at time t based on re-partnership? 
    # # we need Exclude Males who are now unavailable from the catalog of non-repairing individuals
    # for(j in 1:nm){
    #   male_taken_jt[j,t-1] <- sum(vectorMatch(apairs_f[1:nf,t-1],j)*arepartner[1:nf,t-1])
    # }
    
    # 3. Mate Selection -------------------------------------------------------------------------------------------------------------------
    # Use Categorical Distribution to classify mates
    
    # Build Homogeneous Partnership probabilities 
    for(i in 1:nf){
      # Flat likelihood of mating conditional on decision to mate
      # If repairing then force partner to be only choice
      # If not repairing then exclude past partner plus any non-mating males
      for(j in 1:nm){
        psi_cond[i, j, t] <- psi[i,j,t] * amating_f[i,t] * amating_m[j,t] #*
        # (1-equals(apairs_f[i,t-1],j)) * (1-arepartner[i,t-1]) * (1-male_taken_jt[j,t-1]) +
        # arepartner[i,t-1] * equals(apairs_f[i,t-1],j)) 
      }
    }
    
    # Assign Pairs using Custom Random Matrix Sampling Distribution 
    apairs_f[1:nf,t] ~ dpaircat(psi_cond[1:nf,1:nm,t], 
                                amating_f[1:nf,t], 
                                amating_m[1:nm,t], 
                                nf, 
                                nm)
    
    single_female[1:nf,t] <- vectorMatch(apairs_f[1:nf,t], nm + 1) 
    
    # Update Total History for Next Time Step
    # for(i in  1:nf){
    #   for(j in  1:(nm+1)){
    #     histories[i, j, t+1] <- histories[i, j, t] + equals(apairs_f[i,t],j) * (1 - single_female[i,t])
    #   }
    # }
    # 
    # 4. Joint Survival ---------------------------------------------------------------------------------------------------------------------
    
    # Marginal Survival Event for Males in the Population (P[Y^M_T])
    for(j in 1:nm){
      am[j,t] ~ dbern(PhiM * am[j,t-1] * recruit_m[j,t]   + (1-recruit_m[j,t]))
    }
    
    # am[1:nm,t] ~ dmvbern(PhiM * am[1:nm,t-1] * recruit_m[1:nm,t]   + (1-recruit_m[1:nm,t]))
    
    # Marginal Recapture Event for Females in the Population (P[X^F_T|X^M_T]) given males 
    phi.totalF[1:nf,t-1] <- compute_prob_condF(single_female[1:nf,t],
                                               am[1:nm,t],
                                               apairs_f[1:nf,t],
                                               PhiF,
                                               PhiM,
                                               Phifm,
                                               Phif0,
                                               nf,
                                               nm)
    
    # Draw conditional Survival Event
    for(i in 1:nf){
      af[i, t] ~ dbern(phi.totalF[i,t-1] * af[i,t-1] * recruit_f[i,t] + (1-recruit_f[i,t]))
    }
    
    # af[1:nf,t] ~ dmvbern(phi.totalF[1:nf,t-1] * af[1:nf,t-1] * recruit_f[1:nf,t] + (1-recruit_f[1:nf,t]))
  }
  
  # 5. Joint Recapture --------------------------------------------------------------------------------------------------------------------
  for(t in 1:k){
    
    # Marginal Recapture Event for Males in the Population (P[X^M_T])
    for(j in 1:nm){
      recap_m[j,t] ~ dbern(PM * am[j,t] * recruit_m[j,t] * zm[j])
    }
    
    # recap_m[1:nm,t] ~ dmvbern(PM * am[1:nm,t] * recruit_m[1:nm,t] * zm[1:nm])
    
    # Marginal Recapture Event for females in the Population (P[X^F_T|X^M_T])
    p.totalF[1:nf,t] <- compute_prob_condF(single_female[1:nf,t],
                                           recap_m[1:nm,t],
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
    
    # recap_f[1:nf, t] ~ dmvbern(p.totalF[1:nf,t] * af[1:nf,t] * recruit_f[1:nf,t] * zf[1:nf])
  }
  
  
  # 6. Prior Distributions-------------------------------------------------------------------------------------------------------------------
  # Data augmentation
  xi ~ dbeta(0.1,1.0)
  
  # Recruitment 
  for(t in 1:k){
    eps[t] ~ dbeta(1,1)
  }
  
  # Attempt to Mate 
  delta ~ dbeta(1,1)
  
  # # Pairs reforming
  # beta0 ~ dnorm(0, 1)
  # beta1 ~ dnorm(0, 1)
  
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
  gamma <- 0#(gu - gl)*gamma_raw + gl 
  gamma_raw ~ dbeta(1,1)
  
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
  
  # Recapture Rates M/F
  # PhiF ~ dbeta(alpha.phif,beta.phif)
  # # alpha.phif <- beta.phif * odds.phiM
  # # beta.phif ~ dgamma(3,1)
  # #
  # alpha.phif <- PhiM * v.phif
  # beta.phif  <- (1-PhiM) * v.phif
  # v.phif ~ dgamma(10,1)
  # PhiM ~ dbeta(2,2)
  # 
  # 
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
  rho <- 0#(ru - rl)*rho_raw + rl 
  rho_raw ~ dbeta(1,1)
  
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
  
  ### Prior Parameters ####
  
  # Recapture Rates M/F
  # PF ~ dbeta(alpha.pf,beta.pf)
  # # alpha.pf <- beta.pf * odds.PM
  # # beta.pf  ~ dgamma(3,1)
  # alpha.pf <- PM * v.pf
  # beta.pf  <- (1-PM) * v.pf
  # v.pf ~ dgamma(10,1)
  # PM ~ dbeta(2,2)
  
  # 
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
  psi <- ps_data$psi # index who is taken
  
  # CR data with missing components
  zf <- ps_data$zf
  zm <- ps_data$zm
  recruit_f <- ps_data$recruit_f
  recruit_m <- ps_data$recruit_m
  amating_f <- ps_data$amating_f
  amating_m <- ps_data$amating_m
  arepartner <- ps_data$arepartner
  apairs_f <-  ps_data$apairs_f
  af <- ps_data$af
  am <- ps_data$am
  recap_m <- ps_data$recap_m
  na_repartner <- ps_data$na_repartner
  
  # Define local fn equals to emulate jags code
  # Equals call (1 if T; 0 if F)
  equals <- function(a, b){
    return(1*(a == b))
  }
  
  # Recapture Prob and Correlation -------------------------------------------------
  PM <- rbeta(1,2,2)
  v.pf <- rgamma(1,10,1)
  alpha.pf <- PM * v.pf
  beta.pf <- (1-PM) * v.pf
  PF <- rbeta(1,alpha.pf,beta.pf)
  
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
  PhiM <- rbeta(1,2,2)
  v.phif <- rgamma(1,10,1)
  alpha.phif <- PhiM * v.phif
  beta.phif <- (1-PhiM) * v.phif
  PhiF <- rbeta(1,alpha.phif,beta.phif)
  
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
  
  xi <- rbeta(1,0.3,1)
  
  # Missing Data Simulation --------------------------------------------------------
  
  # amating f/m
  # arepartner f/m
  # af/am
  # apairs f/m
  
  # package up the mating stuff into a few functions its a little unweildy
  
  # Sample augmentation female
  for(i in 1:nf){
    zf[i] <- ifelse(is.na(zf[i]), rbinom(1, 1, xi), zf[i])
  }
  
  
  # Sample augmentation male
  for(j in 1:nm){
    zm[j] <- ifelse(is.na(zm[j]), rbinom(1, 1, xi), zm[j])
  }
  
  
  # Recruitment
  # Female Recruitment
  for(i in 1:nf){
    recruit_f[i,1] <- ifelse(is.na(recruit_f[i,1]),
                             rbinom(1, 1, eps[1]),
                             recruit_f[i,1])
    for(t in 2:(k-1)){
      recruit_f[i,t] <- ifelse(is.na(recruit_f[i,t]), 
                               rbinom(1, 1, (recruit_f[i,t-1] + (1-recruit_f[i,t-1]) * eps[t])),
                               recruit_f[i,t])
    }
  }
  
  # Male Recruitment
  for(j in 1:nm){
    recruit_m[j,1] <- ifelse(is.na(recruit_m[j,1]), 
                             rbinom(1, 1, eps[1]), 
                             recruit_m[j,1])
    for(t in 2:(k-1)){
      recruit_m[j,t] <- ifelse(is.na(recruit_m[j,t]), 
                               rbinom(1, 1, (recruit_m[j,t-1] + (1-recruit_m[j,t-1]) * eps[t])),
                               recruit_m[j,t])
    }
  }
  
  # Intermediate objects defined within NIMBLE
  histories <- array(0, dim = c(nf, nm+1, k+1))
  single_female <- matrix(NA, nrow = nf, ncol = k)
  prob_repartner <- matrix(NA, nrow = nf, ncol = k-1)
  phi.totalF <- matrix(NA, nrow = nf, ncol = k-1)
  p.totalF <- matrix(NA, nrow = nf, ncol = k)
  male_taken_jt <- matrix(NA, nrow = nm, ncol = k-1)
  psi_cond <- array(NA, dim = c(nf, nm, k))   
  
  # Time 2 through k initialization
  for(t in 1:k){
    # Female Mating Choice at time t
    for(i in 1:nf){
      if(t == 1){
        amating_f[i,t] <- ifelse(is.na(amating_f[i,t]), 
                                 rbinom(1, 1, recruit_f[i,t] * delta* zf[i]), 
                                 amating_f[i,t])
      } else {
        amating_f[i,t] <- ifelse(is.na(amating_f[i,t]),
                                 rbinom(1, 1, af[i,t-1] * recruit_f[i,t] * delta * zf[i]), 
                                 amating_f[i,t])
      }
      
    }
    
    # Male Mating Choice at time t
    for(j in 1:nm){
      if(t == 1){
        amating_m[j,t] <- ifelse(is.na( amating_m[j,t]),
                                 rbinom(1, 1, recruit_m[j,t] * delta * zm[j]), 
                                 amating_m[j,t])
      } else {
        amating_m[j,t] <- ifelse(is.na( amating_m[j,t]),
                                 rbinom(1, 1, am[j,t-1] * recruit_m[j,t] * delta * zm[j]),  
                                 amating_m[j,t])
      }
      
    }
    
    # Re-partnership happens after time  1
    # arepartner is zero at time  1
    if(t > 1){
      
      # Probability of re-forming 
      prob_repartner[1:nf, t-1] <- compute_pr_repartner(intercept      = beta0,
                                                        slope          = beta1,
                                                        history        = histories[1:nf, 1:(nm+1), t],
                                                        psi_uncond     = psi[1:nf, 1:nm, t],
                                                        mating_f       = amating_f[1:nf,t],
                                                        mating_m       = amating_m[1:nm,t],
                                                        former_pairs_f = apairs_f[1:nf,t-1],
                                                        na_repartner   = na_repartner[1:nf, t-1],
                                                        nf             = nf,
                                                        nm             = nm)
      # Choose to re-form pairs
      for(i in 1:nf){
        arepartner[i,t-1] <- ifelse(is.na(arepartner[i,t-1]), 
                                    rbinom(1,1,prob_repartner[i,t-1]),
                                    arepartner[i,t-1])
        
        lp <- dbinom(arepartner[i,t-1],1,prob_repartner[i,t-1],log=T)
        if(lp == -Inf) browser()
        
      }
      
      # Is Male j taken at time t based on re-partnership? 
      # we need Exclude Males who are now unavailable from the catalog of non-repairing individuals
      for(j in 1:nm){
        male_taken_jt[j,t-1] <- sum(vectorMatch(apairs_f[1:nf,t-1],j)*arepartner[1:nf,t-1])
      }
    }
    
    # Build Homogeneous Partnership probabilities 
    for(i in 1:nf){
      # Flat likelihood of mating conditional on decision to mate
      # If repairing then force partner to be only choice
      # If not repairing then exclude past partner plus any non-mating males
      for(j in 1:nm){
        if(t == 1){
          psi_cond[i, j, t] <- (psi[i,j,t] * amating_f[i,t] * amating_m[j,t])
        } else {
          psi_cond[i, j, t] <- (psi[i,j,t] * amating_f[i,t] * amating_m[j,t] * 
                                  (1-equals(apairs_f[i,t-1],j)) * (1-arepartner[i,t-1]) * (1-male_taken_jt[j,t-1]) + 
                                  arepartner[i,t-1] * equals(apairs_f[i,t-1],j)) 
        }
      }
    }
    
    # Sample from possible partnerships
    apairs_f[1:nf, t] <- rpaircat(n=1,
                                  available_mates = psi_cond[1:nf,1:nm,t],
                                  mating_f = amating_f[1:nf,t],
                                  mating_m = amating_m[1:nm,t],
                                  nf = nf,
                                  nm = nm)
    
    lp <- dpaircat(apairs_f[1:nf,t], 
                   available_mates = psi_cond[1:nf,1:nm,t],
                   mating_f = amating_f[1:nf,t],
                   mating_m = amating_m[1:nm,t],
                   nf = nf,
                   nm = nm,
                   log = 1)
    
    
    # If sampling is going wrong
    if(lp == -Inf) browser()
    
    # Assign single females
    single_female[1:nf,t] <- equals(apairs_f[1:nf,t],nm+1)
    
    # some error handling
    if(any(sort(table(apairs_f[,t]))[-length(table(apairs_f[,t]))] > 1)) stop("Illegal apairing_f at " %+% t)
    
    # Update histories 
    for(i in  1:nf){
      for(j in  1:(nm+1)){
        histories[i, j, t+1] <- histories[i, j, t] + 
          equals(apairs_f[i,t],j)*(1-single_female[i,t])
      }
    }
    
    # If past the first occasion then death is possible 
    if(t > 1){
      # Marginal Survival Event for Males in the Population (P[Y^M_T])---------------------------------------------
      for(j in 1:nm){
        if(is.na(am[j,t])){
          am[j,t] <- rbinom(1,1, PhiM * am[j,t-1] * recruit_m[j,t] + (1-recruit_m[j,t]))
        }
      }
      
      # Marginal Recapture Event for Females in the Population (P[X^F_T|X^M_T])
      phi.totalF[1:nf, t-1] <- compute_prob_condF(is_single_female   = single_female[1:nf,t],
                                                  current_male_state = am[1:nm, t],
                                                  current_pairs_f    = apairs_f[1:nf, t],
                                                  ProbF              = PhiF,
                                                  ProbM              = PhiM,
                                                  Probfm             = Phifm,
                                                  Probf0             = Phif0,
                                                  nf                 = nf,
                                                  nm                 = nm) 
      
      for(i in 1:nf){
        # Draw Survival Event
        if(is.na(af[i,t])){
          af[i, t] <- rbinom(1,1, phi.totalF[i,t-1] * af[i,t-1] * recruit_f[i,t]  + (1-recruit_f[i,t]))
        }
      }
    }
    
    # 5. Joint Recapture --------------------------------------------------------------------------------------------------------------------
    
    # Marginal Recapture Event for females in the Population (P[X^F_T|X^M_T])
    
    p.totalF[1:nf, t] <- compute_prob_condF(is_single_female   = single_female[1:nf,t],
                                            current_male_state = recap_m[1:nm, t],
                                            current_pairs_f    = apairs_f[1:nf, t],
                                            ProbF              = PF,
                                            ProbM              = PM,
                                            Probfm             = Pfm,
                                            Probf0             = Pf0,
                                            nf                 = nf,
                                            nm                 = nm) 
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
  
  # Female Survival
  af <- build_NA_mat(af, ps_data$af)
  
  # Female Mating Status
  amating_f <- build_NA_mat(amating_f, ps_data$amating_f)
  
  # Male Survival
  am <- build_NA_mat(am, ps_data$am)
  
  # Male Mating Status
  amating_m <- build_NA_mat(amating_m, ps_data$amating_m)
  
  # Pair index (female perspective)
  apairs_f <- build_NA_mat(apairs_f, ps_data$apairs_f)
  
  # Repartner index (female perspective)
  arepartner <- build_NA_mat(arepartner, ps_data$arepartner)
  
  zf <- build_NA_vec(zf, ps_data$zf)
  zm <- build_NA_vec(zm, ps_data$zm)
  # Return Results ------------------------------------------------------------------
  
  # Store in object
  ps_inits <- list(
    PF = PF,
    PM = PM,
    # v.pf   = v.pf,
    # v.phif = v.phif,
    rho_raw = rho_raw,
    PhiF = PhiF,
    PhiM = PhiM,
    gamma_raw = gamma_raw,
    eps = eps,
    xi = xi,
    delta = delta,
    beta0 = beta0,
    beta1 = beta1,
    zf = zf,
    zm = zm,
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
    p.totalF = p.totalF,
    male_taken_jt = male_taken_jt,
    prob_repartner = prob_repartner,
    psi_cond = psi_cond
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
      BUGSdist = "dpaircat(available_mates, mating_f, mating_m, nf, nm)",
      Rdist = "dpaircat(available_mates, mating_f, mating_m, nf, nm)",
      discrete = TRUE,
      range = c(1, (ps_data$nm+1)),
      types = c('value = double(1)', 'available_mates = double(2)','mating_f = double(1)','mating_m = double(1)','nf = integer(0)', 'nm = integer(0)'),
      pqAvail = FALSE),
    dmvbern = list(
      BUGSdist = "dmvbern(prob)",
      Rdist = "dmvbern(prob)",
      discrete = TRUE,
      range = c(0, 1),
      types = c('value = double(1)', 'prob = double(1)'),
      pqAvail = FALSE)
  ))
  
  # Generating Initial Values
  cat("Generating Initial Values...", "\n")
  nimble_inits <- generate_nimble_init_pairs(ps_data)
  
  # Construct Nimble Objects 
  cat("Organizing Data for Nimble...", "\n")
  nimble_ps_constants <- list(
    nf = ps_data$nf,
    nm = ps_data$nm,
    k = ps_data$k
  )
  
  nimble_ps_dat <- list(
    zf = ps_data$zf,
    zm = ps_data$zm,
    recruit_f = ps_data$recruit_f,
    recruit_m = ps_data$recruit_m,
    amating_f = ps_data$amating_f,
    amating_m = ps_data$amating_m,
    psi = ps_data$psi,
    af = ps_data$af,
    am = ps_data$am,
    apairs_f = ps_data$apairs_f,
    arepartner = ps_data$arepartner, 
    recap_f = ps_data$recap_f,
    recap_m = ps_data$recap_m,
    na_repartner = ps_data$na_repartner
  )
  
  if(!is.null(params)){
    cat("User-specified Params Detected...","\n")
    cat("Using params := ", "\n")
    cat(params, "\n")
    
    nimble_params <- params
  } else {
    cat("Params argument is NULL...","\n")
    nimble_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps", "gl", "gu", "ru", "rl","Nf","Nm","xi")
    cat("Using params := ", "\n")
    cat(nimble_params, "\n")
  }
  
  
  nimble_dims <- list(histories       = c(nimble_ps_constants$nf, nimble_ps_constants$nm+1, nimble_ps_constants$k+1),
                      prob_repartner  = c(nimble_ps_constants$nf, nimble_ps_constants$k-1),
                      male_taken_jt   = c(nimble_ps_constants$nm, nimble_ps_constants$k-1),
                      psi_cond        = c(nimble_ps_constants$nf, nimble_ps_constants$nm, nimble_ps_constants$k),
                      single_female   = c(nimble_ps_constants$nf, nimble_ps_constants$k),
                      phi.totalF      = c(nimble_ps_constants$nf, nimble_ps_constants$k-1),
                      p.totalF        = c(nimble_ps_constants$nf, nimble_ps_constants$k)
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
  
  # Assign dpaircat sampler
  assign('sampler_pairs', sampler_pairs, envir = .GlobalEnv)
  
  cat("Compiling Graphical Model in C++ (SLOW)...", "\n")
  compile_ps <- compileNimble(psModel, showCompilerOutput = F)
  
  node_names <- psModel$getNodeNames()[psModel$getNodeType(psModel$getNodeNames())=="stoch"]
  node_names <- node_names[!substr(node_names, 1, nchar("apairs_f")) == "apairs_f"]
  # [Note] SafeDepare.... warnings are annoying so suppress messages 
  cat("Configuring Markov Chain Monte Carlo Process (SLOW)...", "\n")
  psConf  <- suppressMessages(configureMCMC(psModel,
                                            print = F,
                                            nodes = node_names,
                                            multivariateNodesAsScalars = T, 
                                            onlySlice = F))
  
  cat("Adding custom random pair sampler to apairs_f...", "\n")
  
  
  psConf$removeSampler("apairs_f", print = F)
  
  for(t in 1:ps_data$k){
    temp_name_cat <- "apairs_f[1:" %+% ps_data$nf %+% "," %+% t %+% "]"
    psConf$addSampler(target = temp_name_cat, type = "sampler_pairs", print = T)
  }
  
  # cat("..also changing recap_f and recap_m be binary...", "\n")
  psConf$removeSampler("recap_m",print = F)
  psConf$removeSampler("recap_f",print = F)
  
  print(psConf)
  
  cat("Adding Monitors and Constructing MCMC...", "\n")
  psConf$addMonitors(nimble_params)
  psMCMC  <- buildMCMC(psConf)
  
  cat("Compiling MCMC Samplers (SLOW)...", "\n")
  CpsMCMC <- compileNimble(psMCMC, project = psModel)
  
  cat("Project Defined, MCMC and Model are compiled...", "\n")
  
  cat("Returning Model Object...", "\n")
  
  return(list(CpsMCMC = CpsMCMC, 
              nimble_inits = nimble_inits))
  
  
}

# Get Samples from Model
run_nimble <- function(CmdlMCMC, 
                       niter,
                       nburnin,
                       thin,
                       inits = NULL,
                       nchains=3,
                       seed = F){
  
  cat("MCMC Sampling from Model...","\n")
  samples <- runMCMC(CmdlMCMC,
                     niter = niter,
                     nburnin = nburnin, 
                     thin = thin,
                     inits = inits,
                     nchains = nchains,
                     setSeed = seed,
                     samplesAsCodaMCMC = TRUE)
  
  cat("Returning Output...","\n")
  return(samples)
}
