# Nimble Functions for Running Pair Swap Model ------------------------------------------------------------------------------------------------------------ 

# Produce vector of 1s and 0s to check for matching value in an existing vector
vectorMatch <- nimbleFunction(
  
  run = function(x= double(1),
                 y = double(0)){
    
    returnType(double(1))
    
    output <- 1*(y == x)
    
    return(output)}
)


# Randomly sample order of distribution for arrival at mating site
darrivalcat <- nimbleFunction(
  run = function(x = double(1),
                 mating_f = double(1),
                 nf = integer(0),
                 log = integer(0, default = 1)){
    
    returnType(double(0))
    
    n_females_mating <- sum(mating_f[1:nf]) + 1
    
    indicies <- seq(0,n_females_mating-1,1)
    logProb <- -sum(log(n_females_mating-indicies)) 
    
    if(log) return(logProb)
    else return(exp(logProb))
  }
)


rarrivalcat <- nimbleFunction(
  run = function(n = integer(0),
                 mating_f = double(1),
                 nf = integer(0)){
    
    returnType(double(1))
    if(n != 1) print("rarrivalcat only allows n = 1; using n = 1.")
    
    prob <- rep(1, nf)
    x <- rep(0, nf)
    for(i in 1:nf){
      x[i] <- rcat(n = 1, prob = prob)
      prob[x[i]] <- 0
    }
    
    
    return(x)
  }
  
)

dpaircat <- nimbleFunction(
  run = function(x = double(1),
                 available_mates = double(2),
                 arrivals = double(1),
                 nf = integer(0),
                 nm = integer(0),
                 log = integer(0, default = 1)){
    
    
    available_mates2 <- available_mates
    returnType(double(0))
    
    ll <- rep(0,nf)
    i_t <- arrivals[1]
    
    if(x[i_t] != (nm+1)){
      ll[1] <- dcat(x[i_t],available_mates2[i_t,1:nm],log = TRUE)
      available_mates2[1:nf, x[i_t]] <- 0
    }
    
    for(i in 2:nf){
      
      i_t <- arrivals[i]
      
      if(x[i_t] != (nm+1)){
        ll[i] <- dcat(x[i_t],available_mates2[i_t,1:nm],log = TRUE)
        available_mates2[1:nf, x[i_t]] <- 0
      }
    }
    
    
    logProb <- sum(ll)
    if(log) return(logProb)
    else return(exp(logProb))
  }
)

rpaircat <- nimbleFunction(
  run = function(n = integer(0),
                 available_mates = double(2),
                 arrivals = double(1),
                 nf = integer(0),
                 nm = integer(0)){
    
    
    returnType(double(1))
    if(n != 1) print("rpaircat only allows n = 1; using n = 1.")
    
    x <- rep(nm+1,nf)
    available_mates2 <- available_mates

    i_t <- arrivals[1]
    
    if(sum(available_mates2[i_t,1:nm]) != 0){
      x[i_t] <- rcat(n=1, prob = available_mates2[i_t, 1:nm])
      available_mates2[1:nf, x[i_t]] <- 0
    }
   
    for(i in 2:nf){
      
      i_t <- arrivals[i]
      
      if(sum(available_mates2[i_t,1:nm]) != 0){
        x[i_t] <- rcat(n=1, prob = available_mates2[i_t, 1:nm])
        available_mates2[1:nf, x[i_t]] <- 0
      }
    }
    
    return(x)
  }
)

compute_pr_repartner <- nimbleFunction(
  run = function(intercept = double(0),
                 slope = double(0),
                 history = double(2),
                 psi_uncond = double(2),
                 mating_f = double(1),
                 mating_m = double(1),
                 single_female = double(1),
                 former_pairs_f  = double(1),
                 na_repartner = double(1),
                 nf = integer(0),
                 nm = integer(0)){
    
    returnType(double(1))
    
    out <- rep(0,nf)
    for(i in 1:nf){
      forced_repartner <- equals(sum(psi_uncond[i,1:nm]),1) * psi_uncond[i,former_pairs_f[i]]* na_repartner[i]
      
      out[i] <- (1-single_female[i]) * mating_f[i] * mating_m[former_pairs_f[i]] * psi_uncond[i, former_pairs_f[i]] * 
        (ilogit(intercept + slope * history[i,former_pairs_f[i]]) * (1-forced_repartner) + forced_repartner)
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
    
    out <- is_single_female     * ProbF + 
      (1-is_single_female) * (current_male_state[current_pairs_f]        * (Probfm/ProbM) +
                                (1- current_male_state[current_pairs_f]) * (Probf0/(1-ProbM)))
    return(out)
  }
)

# Custom Sampler for dpaircat
# Variation on categorical sampler and MH Proposal

sampler_arrivals <- nimbleFunction(
  name = 'sampler_arrivals',
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
    amating_f <- model$getParam(target, 'mating_f')
    logProbs <- numeric(2)
    probs <- numeric(2)
    ## checks
    if(model$getDistribution(target) != 'darrivalcat') stop('can only use pair categorical sampler on node with darrivalcat distribution')
  },
  run = function() {
    current_pairs <- model[[target]]
    logProbs[1] <<- model$getLogProb(calcNodes)
    if(is.nan(logProbs[1])) logProbs[1] <<- -Inf
    # Simulate new partners
    model[[target]] <<- rarrivalcat(1,amating_f, nf) # accept target
    logProbs[2] <<- model$calculate(calcNodes) # calculate logprobs
    if(is.nan(logProbs[2])) logProbs[2] <<- -Inf


    acceptanceProb <- 1/(exp(logProbs[1] - logProbs[2]) + 1)
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
    arrivals <- model$getParam(target, "arrivals")
    psi_cond_t <- model$getParam(target, 'available_mates')
    logProbs <- numeric(2)
    probs <- numeric(2)
    ## checks
    if(model$getDistribution(target) != 'dpaircat') stop('can only use pair categorical sampler on node with dpaircat distribution')
  },
  run = function() {
    current_pairs <- model[[target]]
    logProbs[1] <<- model$getLogProb(calcNodes)
    if(is.nan(logProbs[1])) logProbs[1] <<- -Inf
    # Simulate new partners
    model[[target]] <<- rpaircat(1,psi_cond_t,arrivals, nf, nm) # accept target
    logProbs[2] <<- model$calculate(calcNodes) # calculate logprobs
    if(is.nan(logProbs[2])) logProbs[2] <<- -Inf


    acceptanceProb <- 1/(exp(logProbs[1] - logProbs[2]) + 1)
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


# sampler_arrivals <- nimbleFunction(
#   name = 'sampler_arrivals',
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
#     amating_f <- model$getParam(target, 'mating_f')
#     ## checks
#     if(model$getDistribution(target) != 'darrivalcat') stop('can only use pair categorical sampler on node with darrivalcat distribution')
#   },
#   run = function() {
#     model[[target]] <<- rarrivalcat(1,amating_f, nf) # accept target
#     nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
#     nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
#     nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
#   },
#   methods = list(
#     reset = function() { }
#   )
# )
# 
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
#     arrivals <- model$getParam(target, "arrivals")
#     psi_cond_t <- model$getParam(target, 'available_mates')
#     ## checks
#     if(model$getDistribution(target) != 'dpaircat') stop('can only use pair categorical sampler on node with dpaircat distribution')
#   },
#   run = function() {
#     model[[target]] <<- rpaircat(1,psi_cond_t,arrivals, nf, nm) # accept target
#     nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
#     nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
#     nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
#   },
#   methods = list(
#     reset = function() { }
#   )
# )

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
  
  # # # Initialize History Array (All Zero at time 1)
  for(i in 1:(nf)){
    for(j in 1:(nm+1)){
      histories[i, j, 1] <- 0
    }
  }
  
  # Conditional Partnership/Survival Steps ----------------------------------------------------------------------------------------
  
  # Model Events from t=1 to k --------------------------------------------------------------------------------------------------------------
  
  # Time 1 (all animals must be alive/no animals are re-partnering yet)
  
  # If they've been recruited do they mate at this time step
  for(i in 1:nf){
    amating_f[i,1] ~ dbern(delta * recruit_f[i,1]* zf[i])
  }
  
  for(j in 1:nm){
    amating_m[j,1] ~ dbern(delta * recruit_m[j,1] * zm[j])
  }
  
  # Build Homogeneous Partnership probabilities 
  for(i in 1:nf){
    # Flat likelihood of mating conditional on decision to mate
    for(j in 1:nm){
      psi_cond[i, j, 1] <- psi[i,j,1] * amating_f[i,1] * amating_m[j,1]
    }
  }
  
  # Assign Pairs using Custom Random Matrix Sampling Distribution 
  arrivals[1:nf,1] ~ darrivalcat(amating_f[1:nf,1],
                                 nf)
  
  apairs_f[1:nf,1] ~ dpaircat(psi_cond[1:nf, 1:nm, 1], 
                              arrivals[1:nf,1], 
                              nf, 
                              nm)
  
  single_female[1:nf,1] <- vectorMatch(apairs_f[1:nf,1], nm + 1)
  
  #Update Total History (at start of time 2, what is the history data)
  for(i in  1:nf){
    for(j in  1:(nm+1)){
      histories[i, j, 2] <- histories[i, j, 1] + equals(apairs_f[i,1],j) * (1 - single_female[i,1])
    }
  }
  
  # Time 2-k
  for(t in 2:k){
    
    # 2. Joint Survival ---------------------------------------------------------------------------------------------------------------------
    
    # Marginal Survival Event for Males in the Population (P[Y^M_T])
    for(j in 1:nm){
      am[j,t] ~ dbern(PhiM * am[j,t-1] * recruit_m[j,t-1]  + (1-recruit_m[j,t-1]))
    }
    
    # Marginal Recapture Event for Females in the Population (P[X^F_T|X^M_T]) given males
    phi.totalF[1:nf,t-1] <- compute_prob_condF(single_female[1:nf,t-1],
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
      af[i, t] ~ dbern(phi.totalF[i,t-1] * af[i,t-1] * recruit_f[i,t-1] + (1-recruit_f[i,t-1]))
    }
    
    
    # 3. Decision to Mate -------------------------------------------------------------------------------------------------------------------
    
    # Female Mating Choice at time t
    for(i in 1:nf){
      amating_f[i,t] ~ dbern(delta * af[i,t] * recruit_f[i,t] * zf[i])
    }
    
    # Male Mating Choice at time t
    for(j in 1:nm){
      amating_m[j,t] ~ dbern(delta * am[j,t] * recruit_m[j,t] * zm[j])
    }
    
    prob_repartner[1:nf,t-1] <- compute_pr_repartner(beta0,
                                                     beta1,
                                                     histories[1:nf,1:(nm+1),t],
                                                     psi[1:nf,1:(nm+1),t],
                                                     amating_f[1:nf,t],
                                                     amating_m[1:(nm+1),t],
                                                     single_female[1:nf,t-1],
                                                     apairs_f[1:nf,t-1],
                                                     na_repartner[1:nf,t-1],
                                                     nf,
                                                     nm)

    # Choose to re-form pairs
    for(i in 1:nf){
      arepartner[i,t-1] ~ dbern(prob_repartner[i,t-1])
    }

    #Is Male j from t-1taken at time t based on re-partnership?
    #we need Exclude Males who are now unavailable from the catalog of non-repairing individuals
    for(j in 1:nm){
      male_taken_jt[j,t-1] <- sum(vectorMatch(apairs_f[1:nf,t-1],j)*arepartner[1:nf,t-1])
    }
    
    # 4. Mate Selection -------------------------------------------------------------------------------------------------------------------
    # Use Categorical Distribution to classify mates
    
    # Build Homogeneous Partnership probabilities 
    for(i in 1:nf){
      # Flat likelihood of mating conditional on decision to mate
      # If repairing then force partner to be only choice
      # If not repairing then exclude past partner plus any non-mating males
      for(j in 1:nm){
        psi_cond[i, j, t] <- (psi[i,j,t] * amating_f[i,t] * amating_m[j,t]) * (1-equals(apairs_f[i,t-1],j)) * (1-arepartner[i,t-1]) * (1-male_taken_jt[j,t-1]) +
          arepartner[i,t-1] * equals(apairs_f[i,t-1],j) * (male_taken_jt[j,t-1]) 
      }
    }
    
    # Assign Pairs using Custom Random Matrix Sampling Distribution 
    arrivals[1:nf,t] ~ darrivalcat(amating_f[1:nf,t],
                                   nf)
    
    apairs_f[1:nf,t] ~ dpaircat(psi_cond[1:nf, 1:nm, t], 
                                arrivals[1:nf,t], 
                                nf, 
                                nm)
    
    single_female[1:nf,t] <- vectorMatch(apairs_f[1:nf,t], nm + 1)
    
    # Update Total History for Next Time Step
    for(i in  1:nf){
      for(j in  1:(nm+1)){
        histories[i, j, t+1] <- histories[i, j, t] + equals(apairs_f[i,t],j) * (1 - single_female[i,t])
      }
    }
  }
  
  # 5. Joint Recapture --------------------------------------------------------------------------------------------------------------------
  for(t in 1:k){
    
    # Marginal Recapture Event for Males in the Population (P[X^M_T])
    for(j in 1:nm){
      recap_m[j,t] ~ dbern(PM * am[j,t] * recruit_m[j,t] * zm[j])
    }
    
    # Marginal Recapture Event for females in the Population (P[X^F_T|X^M_T])
    p.totalF[1:nf,t] <- compute_prob_condF(single_female[1:nf,t],
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
  # gamma <- (gu - gl)*raw_gamma + gl #2*raw_gamma - 1
  # raw_gamma ~ dbeta(1,1) 
  
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
  # rho <- (ru - rl)*raw_rho + rl # 
  # raw_rho ~ dbeta(1,1)#dbeta(3,3)
  
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
  arrivals <- ps_data$arrivals 
  
  
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
  
  delta <- rbeta(1, 1, 1)
  
  # delta <- c(ps_data$amating_f[!is.na(ps_data$af) & ps_data$af == 1],ps_data$amating_m[!is.na(ps_data$am) & ps_data$am == 1]) %>% table() %>% prop.table()
  # delta <- delta[2]
  
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
    zf[i] <- ifelse(is.na(zf[i]), 
                    rbinom(1, 1, xi), 
                    zf[i])
  }
  
  
  # Sample augmentation male
  for(j in 1:nm){
    zm[j] <- ifelse(is.na(zm[j]), 
                    rbinom(1, 1, xi), 
                    zm[j])
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
  histories      <- array(0,   dim  = c(nf, nm+1, k+1))
  psi_cond       <- array(NA,  dim = c(nf, nm, k))   
  single_female  <- matrix(NA, nrow = nf, ncol = k)
  prob_repartner <- matrix(NA, nrow = nf, ncol = k-1)
  phi.totalF     <- matrix(NA, nrow = nf, ncol = k-1)
  p.totalF       <- matrix(NA, nrow = nf, ncol = k)
  male_taken_jt  <- matrix(NA, nrow = nm, ncol = k-1)
  forced_repartner <- matrix(NA, nrow = nf, ncol = k-1)
  lp_cat <- list()
  
  # Time 2 through k initialization
  for(t in 1:k){
    
    # If past the first occasion then death is possible 
    if(t > 1){
      # Marginal Survival Event for Males in the Population (P[Y^M_T])---------------------------------------------
      for(j in 1:nm){
        if(is.na(am[j,t])){
          am[j,t] <- rbinom(1,1, PhiM * am[j,t-1] * recruit_m[j,t-1] + (1-recruit_m[j,t-1]))
        }
      }
      
      # Marginal Recapture Event for Females in the Population (P[X^F_T|X^M_T])
      phi.totalF[1:nf, t-1] <- compute_prob_condF(is_single_female   = single_female[1:nf,t-1],
                                                  current_male_state = am[1:(nm+1), t],
                                                  current_pairs_f    = apairs_f[1:nf, t-1],
                                                  ProbF              = PhiF,
                                                  ProbM              = PhiM,
                                                  Probfm             = Phifm,
                                                  Probf0             = Phif0,
                                                  nf                 = nf,
                                                  nm                 = nm) 
      
      for(i in 1:nf){
        # Draw Survival Event
        if(is.na(af[i,t])){
          af[i, t] <- rbinom(1,1, phi.totalF[i,t-1] * af[i,t-1] * recruit_f[i,t-1]  + (1-recruit_f[i,t-1]))
        }
      }
    }
    
    # Female Mating Choice at time t
    for(i in 1:nf){
      if(t == 1){
        amating_f[i,t] <- ifelse(is.na(amating_f[i,t]), 
                                 rbinom(1, 1, recruit_f[i,t] * delta* zf[i]),
                                 amating_f[i,t])
      } else {
        amating_f[i,t] <- ifelse(is.na(amating_f[i,t]),
                                 rbinom(1, 1, af[i,t] * recruit_f[i,t] * delta* zf[i]),
                                 amating_f[i,t])
      }
      
    }
    
    # Male Mating Choice at time t
    for(j in 1:nm){
      if(t == 1){
        amating_m[j,t] <- ifelse(is.na( amating_m[j,t]),
                                 rbinom(1, 1, recruit_m[j,t] * delta* zm[j]),
                                 amating_m[j,t])
      } else {
        amating_m[j,t] <- ifelse(is.na( amating_m[j,t]),
                                 rbinom(1, 1, am[j,t] * recruit_m[j,t] * delta * zm[j]),
                                 amating_m[j,t])
      }
      
    }
    
    # Re-partnership happens after time  1
    # arepartner is zero at time  1
    if(t > 1){
      
      # # Probability of re-forming
      prob_repartner[1:nf, t-1] <- compute_pr_repartner(intercept      = beta0,
                                                        slope          = beta1,
                                                        history        = histories[1:nf, 1:(nm+1), t],
                                                        psi_uncond     = psi[1:nf, 1:(nm+1), t],
                                                        mating_f       = amating_f[1:nf,t],
                                                        mating_m       = amating_m[1:(nm+1),t],
                                                        former_pairs_f = apairs_f[1:nf,t-1],
                                                        single_female  = single_female[1:nf, t-1],
                                                        na_repartner   = na_repartner[1:nf, t-1],
                                                        nf             = nf,
                                                        nm             = nm)
      
      # Choose to re-form pairs
      for(i in 1:nf){
        arepartner[i,t-1] <- ifelse(is.na(arepartner[i,t-1]),
                                    rbinom(1,1,prob_repartner[i,t-1]),
                                    arepartner[i,t-1])
        
        lp <- dbinom(arepartner[i,t-1],1,prob_repartner[i,t-1],log=T)
        if(lp == -Inf|lp == Inf|is.nan(lp)) browser()
        
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
                                  arepartner[i,t-1] * equals(apairs_f[i,t-1],j) * (male_taken_jt[j,t-1])) 
        }
      }
    }
    
    
    arrivals[1:nf,t] <- rarrivalcat(n        = 1, 
                                    mating_f = amating_f[1:nf,t],
                                    nf       = nf)
    
    apairs_f[1:nf,t]<- rpaircat(n               = 1,
                                available_mates = psi_cond[1:nf,1:nm,t],
                                arrivals        = arrivals[1:nf,t], 
                                nf              = nf,
                                nm              = nm)
    
    
    lp_cat[[t]] <- dpaircat(x               = apairs_f[1:nf,t], 
                            available_mates = psi_cond[1:nf,1:nm,t],
                            arrivals        = arrivals[1:nf,t],
                            nf              = nf,
                            nm              = nm,
                            log             = 1)
    
    # If sampling is going wrong
    if(lp_cat[[t]] == -Inf|lp_cat[[t]] == Inf|is.nan(lp_cat[[t]])) browser()
    
    # Assign single females
    single_female[1:nf,t] <- equals(apairs_f[1:nf,t],nm+1)
    
    # some error handling
    if(any(sort(table(apairs_f[,t]))[-length(table(apairs_f[,t]))] > 1)) stop("Illegal apairing_f at " %+% t)
    
    # Update histories 
    for(i in  1:nf){
      for(j in  1:(nm+1)){
        histories[i, j, t+1] <- histories[i, j, t] + equals(apairs_f[i,t],j)*(1-single_female[i,t])
      }
    }
    
    # 5. Joint Recapture --------------------------------------------------------------------------------------------------------------------
    
    # Marginal Recapture Event for females in the Population (P[X^F_T|X^M_T])
    
    p.totalF[1:nf, t] <- compute_prob_condF(is_single_female   = single_female[1:nf,t],
                                            current_male_state = recap_m[1:(nm+1), t],
                                            current_pairs_f    = apairs_f[1:nf, t],
                                            ProbF              = PF,
                                            ProbM              = PM,
                                            Probfm             = Pfm,
                                            Probf0             = Pf0,
                                            nf                 = nf,
                                            nm                 = nm) 
  }
  sum(unlist(lp_cat))
  
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
    raw_rho = raw_rho,
    rho = rho,
    PhiF = PhiF,
    PhiM = PhiM,
    raw_gamma = raw_gamma,
    gamma = gamma,
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
    arrivals = arrivals,
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
      BUGSdist = "dpaircat(available_mates, arrivals, nf, nm)",
      Rdist = "dpaircat(available_mates, arrivals, nf, nm)",
      discrete = TRUE,
      range = c(1, (ps_data$nm+1)),
      types = c('value = double(1)', 'available_mates = double(2)','arrivals = double(1)','nf = integer(0)', 'nm = integer(0)'),
      pqAvail = FALSE),
    darrivalcat = list(
      BUGSdist = "darrivalcat(mating_f, nf)",
      Rdist = "darrivalcat(mating_f, nf)",
      discrete = TRUE,
      range = c(1, (ps_data$nf)),
      types = c('value = double(1)', 'mating_f = double(1)', 'nf = integer(0)'),
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
    amating_f       = ps_data$amating_f,
    amating_m       = ps_data$amating_m,
    af              = ps_data$af,
    am              = ps_data$am,
    apairs_f        = ps_data$apairs_f,
    arepartner      = ps_data$arepartner, 
    recap_f         = ps_data$recap_f,
    recap_m         = ps_data$recap_m,
    psi             = ps_data$psi,
    na_repartner    = ps_data$na_repartner,
    constraint_data = c(1,1),
    arrivals        = ps_data$arrivals
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
  
  
  nimble_dims <- list(histories        = c(nimble_ps_constants$nf, nimble_ps_constants$nm+1, nimble_ps_constants$k+1),
                      prob_repartner   = c(nimble_ps_constants$nf, nimble_ps_constants$k-1),
                      male_taken_jt    = c(nimble_ps_constants$nm, nimble_ps_constants$k-1),
                      psi_cond         = c(nimble_ps_constants$nf, nimble_ps_constants$nm, nimble_ps_constants$k),
                      single_female    = c(nimble_ps_constants$nf, nimble_ps_constants$k),
                      forced_repartner = c(nimble_ps_constants$nf, nimble_ps_constants$k-1),
                      p.totalF         = c(nimble_ps_constants$nf, nimble_ps_constants$k)
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
  assign('sampler_arrivals', sampler_arrivals, envir = .GlobalEnv)
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
  
  
  cat("Adding custom random pair sampler to arrivals", "\n")
  
  # Default sampler wont work here....
  psConf$removeSampler("arrivals", print = F)
  
  for(t in 1:ps_data$k){
    temp_name_arrival <- "arrivals[1:" %+% ps_data$nf %+% "," %+% t %+% "]"
    psConf$addSampler(target = temp_name_arrival, type = "sampler_arrivals", print = T)
  }
  
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
    source(paste0(getwd(), "/Scripts/pair_swap_mod_nimble.R"))
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




# # Pair-Swap random assignment distribution function
# dpaircat <- nimbleFunction(
#   run = function(x = double(1), 
#                  available_mates = double(2), 
#                  nf = integer(0),
#                  nm = integer(0), 
#                  log = integer(0, default = 1)){
#     
#     returnType(double(0))
#     
#     # prob_grid <- available_mates#/(sum(available_mates) + equals(sum(available_mates),0))
#     
#     # females <- rep(1:nf,nm)
#     # males <- numeric(nf * nm)
#     # for(j in 1:nm){
#     #   males[(1:nf) + nf*(j-1)] <- rep(j, nf)
#     # }
#     
#     # n_outcomes <- sum(prob_grid)#sum(1 * (prob_grid > 0))
#     
#     ll <- numeric(nf)
#     
#     
#     # for(i in 1:nf){
#     #   if(x[i] != (nm+1)){
#     #     lli <- -log((sum(prob_grid) + equals(sum(prob_grid),0)))
#     #     ll[i] <- lli 
#     #     prob_grid[i,1:nm] <- 0
#     #     prob_grid[1:nf,x[i]] <- 0
#     #   } else {
#     #     potential_mates <- available_mates[i,1:nm]
#     #     N_potential <- sum(potential_mates) 
#     #     if(N_potential > 1.0){
#     #     #   
#     #     #   index_potential <- (1:nm) * potential_mates
#     #     #   index_potential <- index_potential[index_potential != 0]
#     #     #   
#     #     #   xlli <- 1
#     #     #   
#     #     #   for(j in 1:N_potential){
#     #     #     Nj_potential <- sum(available_mates[1:nf,index_potential[j]])
#     #     #     xlli <- xlli * ((N_potential)/(N_potential+1)) * ((Nj_potential)/(Nj_potential+1))  
#     #     #   }
#     #     #   
#     #     #   ll[i] <- log(pmin(pmax(0, 1 - xlli), 1))
#     #       ll[i] <- N_potential * log(1 - 1/N_potential)
#     #         #dbinom(0,size = N_potential, 1/N_potential)
#     #     }
#     #   }
#     # }
#     
#     logProb <- 0#sum(ll)
#     if(log) return(logProb)
#     else return(exp(logProb))
#   }
# )
# 
# 
# 
# # Pair-Swap random assignment random value generator function
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
#     
#     prob_grid <- available_mates#/(sum(available_mates) + equals(sum(available_mates),0))
#     
#     females <- rep(1:nf,nm)
#     males <- numeric(nf * nm)
#     for(j in 1:nm){
#       males[(1:nf) + nf*(j-1)] <- rep(j, nf)
#     }
#     
#     n_outcomes <- sum(prob_grid)#sum(1 * (prob_grid > 0))
#     
#     y <- rep(0,nf*nm)
#     ij <- 0
#     probVect <- numeric(nf * nm)
#       
#     while(n_outcomes > 0){
#       for(j in 1:nm){
#         probVect[1:nf + nf * (j-1)] <- prob_grid[1:nf,j] 
#       }
#       ij <- rcat(n = 1,prob = probVect)
#       prob_grid[females[ij],1:nm] <- 0
#       prob_grid[1:nf,males[ij]] <- 0
#       y[ij] <- 1
#       # prob_grid <- prob_grid/(sum(prob_grid) + equals(sum(prob_grid),0))
#       n_outcomes <- sum(prob_grid)#sum(1 * (prob_grid > 0))
#     }
#     
#     x[females[y==1]] <- males[y==1]
#     
#     
#     return(x)
#   }
# )
# 
# # Pair-Swap random assignment distribution function
# dpaircat <- nimbleFunction(
#   run = function(x = double(1),
#                  available_mates = double(2),
#                  mating_f = double(1),
#                  mating_m = double(1),
#                  nf = integer(0),
#                  nm = integer(0),
#                  log = integer(0, default = 1)){
#     
#     
#     available_mates2 <- matrix( 0, nrow = nf+1, ncol = nm + 1)
#     returnType(double(0))
#     
#     # How many available males and females
#     mating_females <- sum(mating_f) + 1
#     mating_males  <- sum(mating_m) + 1
#     osr <- mating_females/mating_males
#     
#     if(osr <= 1){
#       ll <- rep(0, nf)
#       # possible_mates <- rep(0,nf)
#       available_mates2[1, 1:(nm+1)] <- c(available_mates[1,1:nm], equals(sum(available_mates[1,1:nm]),0)) #+ 1e-6
#       ll[1] <- dcat(x[1],available_mates2[1,1:(nm+1)],log = TRUE)
#       
#       for(i in 2:nf){
#         
#         # Remove Formed Pairs
#         for(j in 1:nm){
#           available_mates2[i,j] <- available_mates[i,j]*(1-1*any(j == x[1:(i-1)]))
#         }
#         
#         #Add case in which no pairs are available
#         available_mates2[i,(nm+1)] <-  equals(sum(available_mates2[i,1:nm]),0) #+ 1e-6
#         
#         ll[i] <- dcat(x[i],available_mates2[i,1:(nm+1)],log = TRUE)
#         
#         
#       }
#     } else {
#       y <- rep(nf+1, nm)
#       ll <- rep(0,nm)
#       
#       # Build vector of male pairs
#       # Need base likelihood on these values
#       for(j in 1:nm){
#         # Find which is your mate
#         for(i in 1:nf){
#           # if(is.nan(x[i])|is.na(x[i])) x[i] <- nm+1
#           if(x[i] == j){
#             y[j] <- i
#           }
#         }
#       }
#       
#       
#       # Apply mating selection for first male
#       available_mates2[1:(nf+1), 1] <- c(available_mates[1:nf,1], equals(sum(available_mates[1:nf,1]),0)) #+ 1e-6
#       ll[1] <- dcat(y[1],available_mates2[1:(nf+1),1],log = TRUE)
#       
#       # Iterate over remaining males
#       for(j in 2:nm){
#         
#         # Remove Formed Pairs
#         for(i in 1:nf){
#           available_mates2[i,j] <- available_mates[i,j]*(1-1*any(i == y[1:(j-1)]))
#         }
#         
#         #Add case in which no pairs are available
#         available_mates2[(nf+1),j] <- equals(sum(available_mates2[1:nf,j]),0) #+ 1e-6
#         
#         # Find mate for i
#         ll[j] <- dcat(y[j],available_mates2[1:(nf+1),j],log = TRUE)
#         
#       }
#     }
#     
#     logProb <- sum(ll)
#     if(log) return(logProb)
#     else return(exp(logProb))
#   }
# )
# 
# # 
# # # Pair-Swap random assignment random value generator function
# rpaircat <- nimbleFunction(
#   run = function(n = integer(0),
#                  available_mates = double(2),
#                  mating_f = double(1),
#                  mating_m = double(1),
#                  nf = integer(0),
#                  nm = integer(0)){
#     
#     x <- rep(nm+1,nf)
#     
#     available_mates2 <- matrix(value = 0, nrow = nf+1, ncol = nm + 1)
#     
#     returnType(double(1))
#     if(n != 1) print("rpaircat only allows n = 1; using n = 1.")
#     
#     # How many available males and females
#     mating_females <- sum(mating_f) + 1
#     mating_males  <- sum(mating_m) + 1
#     osr <- mating_females/mating_males
#     
#     # If operating sex ratio is balanced or in favor of females
#     # Iterate across female
#     if(osr <= 1){
#       
#       available_mates2[1, 1:(nm+1)] <- c(available_mates[1,1:nm],equals(sum(available_mates[1,1:nm]),0)) #+ 1e-6
#       x[1] <- rcat(n=1, prob = available_mates2[1, 1:(nm+1)])
#       
#       for(i in 2:nf){
#         
#         # Remove Formed Pairs
#         for(j in 1:nm){
#           available_mates2[i,j] <- available_mates[i,j]*(1-1*any(j == x[1:(i-1)]))
#         }
#         
#         #Add case in which no pairs are available
#         available_mates2[i,(nm+1)] <- equals(sum(available_mates2[i,1:nm]),0) + 1e-6
#         
#         # Find mate for i
#         x[i] <- rcat(n=1, prob = available_mates2[i, 1:(nm+1)])
#         
#       }
#       
#     } else {
#       # Otherwise OSR favors males so we iterate across their mating pop instead
#       y <- rep(nf+1,nm)
#       
#       # Apply mating selection for first male
#       available_mates2[1:(nf+1), 1] <- c(available_mates[1:nf,1], equals(sum(available_mates[1:nf,1]),0)) #+ 1e-6
#       y[1] <- rcat(n = 1, prob = available_mates2[1:(nf+1),1])
#       
#       # Iterate over remaining males
#       for(j in 2:nm){
#         
#         # Remove Formed Pairs
#         for(i in 1:nf){
#           available_mates2[i,j] <- available_mates[i,j]*(1-1*any(i == y[1:(j-1)]))
#         }
#         
#         #Add case in which no pairs are available
#         available_mates2[(nf+1),j] <- equals(sum(available_mates2[1:nf,j]),0) #+ 1e-6
#         
#         # Find mate for i
#         y[j] <- rcat(n= 1, prob = available_mates2[1:(nf+1), j])
#         
#       }
#       
#       # Mapping male choices back to the female vector
#       for(i in 1:nf){
#         for(j in 1:nm){
#           if(y[j] == i){
#             x[i] <- j
#           }
#         }
#       }
#     }
#     
#     return(x)
#   }
# )