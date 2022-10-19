# Nimble Functions for Running Pair Swap Model ------------------------------------------------------------------------------------------------------------ 

# Produce vector of 1s and 0s to check for matching value in an existing vector
vectorMatch <- nimbleFunction(
  run = function(x= double(1),
                 y = double(0)){
    returnType(double(1))
    output <- 1*(y == x)
    return(output)}
)

compute_pair_probs <- nimbleFunction(
  run = function(psi_slice = double(1),
                 af_i = double(0),
                 am_slice = double(1),
                 hist_slice = double(1),
                 beta0 = double(0),
                 beta1 = double(0),
                 nm = integer(0)){
    
    returnType(double(1))
    
    probs <- c(ilogit(beta0 + beta1 * hist_slice[1:nm]) *  af_i * am_slice[1:nm],0)
    
    if(sum(probs)==0|sum(psi_slice[1:nm]) == 0){
      probs[nm+1] <- 1
    }
    
    return(probs[1:(nm+1)]/sum(probs[1:(nm+1)]))
  }
)

compute_partner_status <- nimbleFunction(
  run = function(apairs_f_i    = double(0),
                 partner_slice = double(1),
                 nm            = integer(0)){
    returnType(double(0))
    if(apairs_f_i == (nm+1)){
      status <- 1
    } else {
      status <- 2 + partner_slice[apairs_f_i]
    } 
    
    return(status)
    
  }
)

update_pair_history <- nimbleFunction(
  run = function(hist_last  = double(1),
                 apairs_f_i = double(0),
                 nm         = integer(0)){
    returnType(double(1))
    
    hist_i <- hist_last
    
    if(apairs_f_i != (nm+1)){
      hist_i[apairs_f_i] <- hist_i[apairs_f_i] + 1
    }
    
    return(hist_i)
  }
)

# BUGS/JAGS Code
nimble_ps_model <- nimbleCode({
  
  #0. Assign Partners at t-------------------------------------------------------------------------------------------------------------------
  
  # Initialize Histories
  for(i in 1:nf){
    for(j in 1:nm){
      hist[i,j,first_capture_f[i]] <- 0
    }
  }
  
  for(i in 1:nf){
    for(t in first_capture_f[i]:(k-1)){
      pair_probs[i, 1:(nm+1) ,t] <- compute_pair_probs(psi[i,1:nm,t],
                                                       af[i,t],
                                                       am[1:nm,t],
                                                       hist[i,1:nm,t],
                                                       beta0,
                                                       beta1,
                                                       nm)
      
      apairs_f[i,t] ~ dcat(pair_probs[i,1:(nm+1),t])
      partner_status[i,t]      <- compute_partner_status(apairs_f[i,t],  am[1:nm,t+1],    nm)
      partner_observation[i,t] <- compute_partner_status(apairs_f[i,t],  recap_m[1:nm,t], nm)
      hist[i,1:nm,t+1]         <- update_pair_history(hist[i,1:nm,t], apairs_f[i,t],   nm)
    }
    
    pair_probs[i, 1:(nm+1) ,k] <- compute_pair_probs(psi[i,1:nm,k],
                                                     af[i,k],
                                                     am[1:nm,k],
                                                     hist[i,1:nm,k],
                                                     beta0,
                                                     beta1,
                                                     nm)
    apairs_f[i,k] ~ dcat(pair_probs[i,1:(nm+1),k])
    partner_observation[i,k] <- compute_partner_status(apairs_f[i,k],  recap_m[1:nm,k], nm)
  }
  
  # 1. Joint Survival [t-1,t) ---------------------------------------------------------------------------------------------------------------
  
  # Marginal Survival Event for Males in the Population (P[Y^M_T]) 
  for(j in 1:nm){
    for(t in first_capture_m[j]:(k-1)){
      am[j,t+1] ~ dbern(PhiM * am[j,t])
    }
  }
  
  # Draw conditional Survival Event
  for(i in 1:nf){
    for(t in first_capture_f[i]:(k-1)){
      af[i, t+1] ~ dbern(Phi_Vector_f[partner_status[i,t]] * af[i,t])
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
      recap_f[i, t] ~ dbern(P_Vector_f[partner_observation[i,t]] * af[i,t])
    }
  }
  
  # 5. Prior Distributions-------------------------------------------------------------------------------------------------------------------
  
  # Partnership Terms
  beta0 ~ dnorm(0,1)
  beta1 ~ dnorm(0,1)
  
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
  apairs_f        <- ps_data$apairs_f
  psi             <- ps_data$psi
  
  # Partner Prob -------------------------------------------------------------------
  beta0 <- rnorm(1,0,1)
  beta1 <- rnorm(1,0,1)
  
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
  
  
  # Build Derived Nodes
  pair_probs          <- array(NA, dim = c(nf, nm+1, k))
  partner_status      <- matrix(NA, nrow = nf, ncol = k)
  partner_observation <- matrix(NA, nrow = nf, ncol = k)
  hist                <- array(0,   dim = c(nf, nm, k))
  
  # Randomly Sample outcomes using initial values
  for(t in 1:k){
    # Survival Outcome ----------------------------------------------------------
    if(t < k){
      
      
      # Assign Partner Status/Observation
      for(i in 1:nf){
        if(t >= first_capture_f[i]){
          pair_probs[i,1:(nm+1),t] <- compute_pair_probs(psi[i,1:nm,t],
                                                         af[i,t],
                                                         am[1:nm,t],
                                                         hist[i,1:nm,t],
                                                         beta0,
                                                         beta1,
                                                         nm)
          
          if(is.na(apairs_f[i,t])){
            apairs_f[i,t] <- rcat(n = 1, pair_probs[i,1:(nm+1),t])
          }
          
          hist[i,1:nm,t+1]           <- update_pair_history(hist[i,1:nm,t],apairs_f[i,t],nm)
          partner_observation[i,t] <- compute_partner_status(apairs_f[i,t],  recap_m[1:nm,t], nm)
        }
      }
      
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
      
      # Conditional Survival Event for Females in the Population ([X^F_T|X^M_T])
      for(i in 1:nf){
        if(t >= first_capture_f[i]){
          partner_status[i,t] <- compute_partner_status(apairs_f[i,t],  am[1:nm,t+1],    nm)
          af[i, t+1] <- ifelse(is.na(af[i,t+1]),
                               rbinom(1,1, (Phi_Vector_f[partner_status[i,t]]) * af[i,t]),
                               af[i,t+1])
        }
        
        # Error check females
        if(any(is.na(af[i,t+1]))) browser()
      }
    } else {
      for(i in 1:nf){
        if(t >= first_capture_f[i]){
          pair_probs[i,1:(nm+1),t] <- compute_pair_probs(psi[i,1:nm,t],
                                                         af[i,t],
                                                         am[1:nm,t],
                                                         hist[i,1:nm,t],
                                                         beta0,
                                                         beta1,
                                                         nm)
          
          if(is.na(apairs_f[i,t])){
            apairs_f[i,t] <- rcat(n = 1, pair_probs[i,1:(nm+1),t])
          }
          
          partner_observation[i,t] <- compute_partner_status(apairs_f[i,t],  recap_m[1:nm,t], nm)
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
  
  # Female Survival
  af <- build_NA_mat(af, ps_data$af)
  
  # Male Survival
  am <- build_NA_mat(am, ps_data$am)
  
  # Pair Status
  apairs_f <- build_NA_mat(apairs_f, ps_data$apairs_f)
  
  # Return Results ------------------------------------------------------------------
  
  # Store in object
  ps_inits <- list(
    beta0             = beta0,
    beta1             = beta1,
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
    apairs_f          = apairs_f
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
    psi             = ps_data$psi
  )
  
  nimble_ps_dat <- list(
    af              = ps_data$af,
    am              = ps_data$am,
    recap_f         = ps_data$recap_f,
    recap_m         = ps_data$recap_m,
    apairs_f        = ps_data$apairs_f,
    constraint_data = c(1,1)
  )
  
  if(!is.null(params)){
    cat("User-specified Params Detected...","\n")
    cat("Using params := ", "\n")
    cat(params, "\n")
    
    nimble_params <- params
  } else {
    nimble_params <- c("PF","PM","PhiF","PhiM", "beta0", "beta1",
                       "gl","gu","gamma", 
                       "ru","rl","rho")
    cat("Params argument is NULL...","\n")
    cat("Using params := ", "\n")
    cat(nimble_params, "\n")
  }
  
  nimble_dims <- list(
    P_Vector_f          = c(3),
    Phi_Vector_f        = c(3),
    hist                = c(ps_data$nf, ps_data$nm, ps_data$k),
    pair_probs          = c(ps_data$nf, ps_data$nm+1, ps_data$k),
    partner_status      = c(ps_data$nf, ps_data$k),
    partner_observation = c(ps_data$nf, ps_data$k)
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
  
  cat("Adding AF-Slice Sampler to Recapture and Survival Parameters")
  psConf$removeSampler(c("PhiF","PhiM","raw_gamma"), print = F)
  psConf$removeSampler(c("PF","PM","raw_rho"), print = F)
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
    source(paste0(getwd(), "/Scripts/pair_swap_mod_nimble9.R"))
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