# Nimble Functions for Running Pair Swap Model ------------------------------------------------------------------------------------------------------------ 

# Produce vector of 1s and 0s to check for matching value in an existing vector
vectorMatch <- nimbleFunction(
  run = function(x= double(1),
                 y = double(0)){
    returnType(double(1))
    output <- 1*(y == x)
    return(output)}
)

dpair <- nimbleFunction(
  run = function(x = double(1),
                 prob = double(0),
                 male_state = double(1),
                 male_state_previous = double(1),
                 female_state_previous = double(1),
                 known_states = double(1),
                 partner_state_recap = double(1),
                 nf = integer(0),
                 log = integer(0, default = 1)){
    returnType(double(0))
    ll <- 0
    if(!log) ll <- exp(ll)
    return(ll)
  }
)

rpair <- nimbleFunction(
  run = function(n = integer(0),
                 prob = double(0),
                 male_state = double(1),
                 male_state_previous = double(1),
                 female_state_previous = double(1),
                 known_states = double(1),
                 partner_state_recap = double(1),
                 nf = integer(0)){
    returnType(double(1))
    if(n != 1) print("rpair only allows n = 1; using n = 1.")
    
    
    # Initialize Counts for +1/0 partners 
    x <- known_states
    
    # Number of females who have been observed with a state 2 or 3 partner
    num_2 <- sum(vectorMatch(x,2))
    num_3 <- sum(vectorMatch(x,3))
    
    # Number of Males in State 2 & state 3, conditional on survival
    NM <- sum(male_state_previous)
    N3 <- sum(male_state) 
    N2 <- NM - N3
    
    # Number of males who are in state 2/3 but have NOT been observed with a female
    N3 <- max(0, N3 - num_3)
    N2 <- max(0, N2 - num_2)
    
    # Randomly order females 
    prob_arrival <- rep(1, nf) 
    arrivals <- rep(0, nf)
    for(i in 1:nf){
      arrivals[i] <- rcat(n = 1, prob = prob_arrival)
      prob_arrival[arrivals[i]] <- 0
    }
    
    # Go through random order and assign possible states
    for(i in 1:nf){
      
      # Which female arrived at i? 
      it <- arrivals[i]
      
      # If X[it] is unknown then assign partner status
      if(x[it] == 0){
        
        # If only state 2 males exist 
        if(N3 == 0 & N2 != 0){
          x[it] <- 1 + female_state_previous[it]
        }
        
        # If only state 3 males exist
        if(N3 != 0 & N2 == 0){
          x[it] <- 1 +  2 * female_state_previous[it]
        }
        
        # If both state 2 and state 3 males are available
        if(N3 > 0 & N2 > 0){
          x[it] <- 1 + female_state_previous[it] * (1 + rbinom(1, 1, prob))
        }
        
        
        # If no males are left then sample single
        if(N3 == 0 & N2 == 0){
          x[it] <- 1
        }
        
        # If during a recapture draw, the survival status of f's partner was 1/2 then she must be single
        if(partner_state_recap[it] == 1){
          x[it] <- 1
        }
        
        # If during a recapture draw, the survival status of f's partner was 3 then she must be paired 
        if(partner_state_recap[it] > 1 & x[it] == 1){
          x[it] <- 1 + female_state_previous[it] * (1 + rbinom(1, 1, prob))
        }
        
        # Update availability 
        # For females who are known this is done outside of the loop
        if(x[it] == 3){
          N3 <- max(0, N3 - 1)
        }
        
        if(x[it] == 2){
          N2 <- max(0, N2 - 1)
        }
      }
      
    }
    
    return(x)
    
  }
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
    prob <- model$getParam(target, 'prob')
    male_state <- model$getParam(target, 'male_state')
    male_state_previous <- model$getParam(target, 'male_state_previous')
    female_state_previous <- model$getParam(target, 'female_state_previous')
    known_states  <- model$getParam(target, 'known_states')
    partner_state_recap <- model$getParam(target, "partner_state_recap")
    nf <- model$getParam(target, 'nf')
    ## checks
    if(model$getDistribution(target) != 'dpair') stop('can only use pair categorical sampler on node with dpair distribution')
  },
  run = function() {
    currentLogProb <- model$getLogProb(calcNodes)
    model[[target]] <<- rpair(1,prob, male_state, male_state_previous, female_state_previous, known_states, partner_state_recap, nf)
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
  
  # 1. Joint Survival [t-1,t) ---------------------------------------------------------------------------------------------------------------
  
  # Marginal Survival Event for Males in the Population (P[Y^M_T]) 
  for(j in 1:nm){
    for(t in first_capture_m[j]:(k-1)){
      am[j,t+1] ~ dbern(PhiM * am[j,t])
    }
  }
  
  # Apply potential partner status
  for(t in 1:(k-1)){
    az[1:nf, t] ~ dpair(0, am[1:nm,t+1], am[1:nm,t],af[1:nf,t], az_known[1:nf,t], ar[1:nf,t], nf) 
  }
  
  # Draw conditional Survival Event
  for(i in 1:nf){
    for(t in first_capture_f[i]:(k-1)){
      af[i, t+1] ~ dbern(Phi_Vector_f[az[i,t]] * af[i,t])
    }
  }
  
  # 4. Joint Recapture --------------------------------------------------------------------------------------------------------------------
  
  # Marginal Recapture Event for Males in the Population (P[X^M_T])
  for(j in 1:nm){
    for(t in (first_capture_m[j]+1):(k)){
      recap_m[j,t] ~ dbern(PM * am[j,t])
    }
  }
  
  
  # Apply potential partner status
  for(t in 1:k){
    ar[1:nf, t] ~ dpair(0, recap_m[1:nm,t], am[1:nm,t],af[1:nf,t],ar_known[1:nf,t], ar_known[1:nf,t] * 0 , nf) 
  }
  
  
  # Draw Recapture Probability
  for(i in 1:nf){
    for(t in (first_capture_f[i]+1):(k)){
      # ar[i, t-1] ~ dcat(r_xi[1:3])
      recap_f[i, t] ~ dbern(P_Vector_f[ar[i, t]] * af[i,t])
    }
  }
  
  # 5. Prior Distributions-------------------------------------------------------------------------------------------------------------------
  # Survival Terms
  
  ### Derived Parameters ####
  Phi_Vector_f[1] <- PhiF
  Phi_Vector_f[2] <- Phif_M0
  Phi_Vector_f[3] <- Phif_M1
  
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
  P_Vector_f[1] <- PF
  P_Vector_f[2] <- Pf_M0
  P_Vector_f[3] <- Pf_M1
  
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

  #Unpack Variables -----------------------------------------------------------------
  # Indexes
  k               <- ps_data$k
  nf              <- ps_data$nf
  nm              <- ps_data$nm
  first_capture_f <- ps_data$first_capture_f
  first_capture_m <- ps_data$first_capture_m  
  recap_f         <- ps_data$recap_f
  recap_m         <- ps_data$recap_m 
  af              <- ps_data$af
  am              <- ps_data$am
  az_known        <- ps_data$az_known
  ar_known        <- ps_data$ar_known
  ar              <- matrix(NA, nrow = nf, ncol = k)
  az              <- matrix(NA, nrow = nf, ncol = k) 
  
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
  P_Vector_f <- c(PF, Pf_M0, Pf_M1)
  
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
  Phi_Vector_f <- c(PhiF, Phif_M0, Phif_M1)

  # Randomly Sample outcomes using initial values
  for(t in 1:k){
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
      
      # Error check males
      if(any(is.na(am[j,t+1]))) browser()
      
      # Assign male partners via states to females
      ar[1:nf, t] <- rpair(n = 1, 0, recap_m[1:nm,t], am[1:nm,t],af[1:nf,t], ar_known[1:nf,t], rep(0,nf), nf) 
      az[1:nf,t] <- rpair(n = 1, 0, am[1:nm,t+1], am[1:nm,t],af[1:nf,t], az_known[1:nf,t], ar[1:nf,t], nf) 
    
      
      # Conditional Survival Event for Females in the Population ([X^F_T|X^M_T])
      for(i in 1:nf){
        if(t >= first_capture_f[i]){
          af[i, t+1] <- ifelse(is.na(af[i,t+1]),
                               rbinom(1,1, (Phi_Vector_f[az[i,t]]) * af[i,t]),
                               af[i,t+1])
        }
        
        # Error check females
        if(any(is.na(af[i,t+1]))) browser()
      }
    } else {
      ar[1:nf, t] <- rpair(n = 1, 0, recap_m[1:nm,t], am[1:nm,t],af[1:nf,t], ar_known[1:nf,t], rep(0,nf), nf) 
    }
  }
  
  # # Add male partners via observations to females
  # for(t in 2:k){
  #   ar[1:nf, t] <- rpair(n = 1, 0.4, recap_m[1:nm,t], am[1:nm,t],af[1:nf,t], ar_known[1:nf,t], nf) 
  # }
  # 
  
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
    af                = af,
    am                = am,
    az                = az,
    ar                = ar
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
    dpair = list(
      BUGSdist = "dpair(prob, male_state, male_state_previous, female_state_previous, known_states, partner_state_recap, nf)",
      Rdist = "dpair(prob, male_state, male_state_previous, female_state_previous, known_states, partner_state_recap, nf)",
      discrete = TRUE,
      range = c(1, 3),
      types = c('value = double(1)', 
                'prob = double(0)', 
                'male_state = double(1)',
                'male_state_previous = double(1)',
                'female_state_previous = double(1)',
                'known_states = double(1)',
                'partner_state_recap = double(1)',
                'nf = integer(0)'),
      pqAvail = FALSE)
  ))
  
  
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
    az_known        = ps_data$az_known,
    ar_known        = ps_data$ar_known
  )
  
  nimble_ps_dat <- list(
    af              = ps_data$af,
    am              = ps_data$am,
    recap_f         = ps_data$recap_f,
    recap_m         = ps_data$recap_m,
    az              = matrix(NA, nrow = ps_data$nf, ncol = ps_data$k), 
    ar              = matrix(NA, nrow = ps_data$nf, ncol = ps_data$k), 
    constraint_data = c(1,1)
  )
  
  if(!is.null(params)){
    cat("User-specified Params Detected...","\n")
    cat("Using params := ", "\n")
    cat(params, "\n")
    
    nimble_params <- params
  } else {
    nimble_params <- c("PF","PM","PhiF","PhiM",
                       "gl","gu","gamma", 
                       "ru","rl","rho")
    cat("Params argument is NULL...","\n")
    cat("Using params := ", "\n")
    cat(nimble_params, "\n")
  }
  
  nimble_dims <- list(
    P_Vector_f    = c(3),
    Phi_Vector_f  = c(3)
  )
  
  cat("Building Model Nodes in Nimble (SLOW)...", "\n")
  psModel <- nimbleModel(code       = nimble_ps_model, 
                         constants  = nimble_ps_constants, 
                         inits      = nimble_inits,
                         data       = nimble_ps_dat,
                         dimensions = nimble_dims)
  
  
  lp_init <- psModel$calculate()
  print(paste0("LP from initial values is ", round(lp_init,3)))
  
  cat("Compiling Graphical Model in C++ (SLOW)...", "\n")
  compile_ps <- compileNimble(psModel, showCompilerOutput = F)
  
  # [Note] SafeDepare.... warnings are annoying so suppress messages 
  # Conjugacy is slow to detect and not useful here so turn off
  cat("Configuring Markov Chain Monte Carlo Process (SLOW)...", "\n")
  psConf  <- suppressMessages(configureMCMC(model = psModel,
                                            print = F,
                                            multivariateNodesAsScalars = T, 
                                            monitors = nimble_params,
                                            onlySlice = F,
                                            useConjugacy = F))
  
  #Assign dpair sampler to the global environment
  assign('sampler_pairs', sampler_pairs, envir = .GlobalEnv)
  
  cat("Adding custom random pair sampler to az and ar...", "\n")

  # Default sampler wont work here....
  psConf$removeSampler("az", print = F)
  psConf$removeSampler("ar", print = F)
  
  for(t in 1:(ps_data$k-1)){
    temp_name_cat <- "az[1:" %+% ps_data$nf %+% "," %+% t %+% "]"
    psConf$addSampler(target = temp_name_cat, type = "sampler_pairs", print = F)
    rm(temp_name_cat)
  }
  
  for(t in 2:ps_data$k){
    temp_name_cat <- "ar[1:" %+% ps_data$nf %+% "," %+% t %+% "]"
    psConf$addSampler(target = temp_name_cat, type = "sampler_pairs", print = F)
    rm(temp_name_cat)
  }
  
  cat("Adding AF-Slice Sampler to Recapture and Survival Parameters")
  psConf$removeSampler(c("PhiF","PhiM","raw_gamma"), print = T)
  psConf$removeSampler(c("PF","PM","raw_rho"), print = T)
  psConf$addSampler(target = c("PhiF","PhiM","raw_gamma"), type = 'AF_slice')
  psConf$addSampler(target = c("PF","PM","raw_rho"), type = 'AF_slice')
  
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
    source(paste0(getwd(), "/Scripts/pair_swap_mod_nimble8.R"))
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