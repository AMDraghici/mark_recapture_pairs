# Nimble Functions for Running Cormack Jolly Seber Model ------------------------------------------------------------------------------------------------------------ 

# BUGS/JAGS Code
nimble_cjs_model <- nimbleCode({
  
  # CJS Likelihood -----------------------------------------------------------------------------------
  for(i in 1:n){
    for(t in initial_entry[i]:(k-1)){
      ## 1) Survival
      a[i,t+1] ~ dbern(phi[i] * a[i,t])
    }
  }
  
  for(i in 1:n){
    for(t in (initial_entry[i]+1):k){
      x[i,t] ~ dbern(p[i] * a[i,t])
    }
  }
  
  # Priors--------------------------------------------------------------------------------------------
  
  # Split probs by sex 
  p[1:n]   <- female[1:n] * PF   + (1-female[1:n]) * PM
  phi[1:n] <- female[1:n] * PhiF + (1-female[1:n])* PhiM
  
  # Survival by sex
  PhiF ~ dbeta(1,1)
  PhiM ~ dbeta(1,1)
  
  # Recapture by sex
  PF ~ dbeta(1,1)
  PM ~ dbeta(1,1)
})


generate_init_cjs <- function(cjs_data){
  
  #Unpack Variables -----------------------------------------------------------------
  
  # Known data and indices 
  k             <- cjs_data$k # Number of capture occasions
  n             <- cjs_data$n # Number of animals 
  female        <- cjs_data$female # Sex
  x             <- cjs_data$x # Recapture 
  initial_entry <- cjs_data$initial_entry
  # CR data with missing components
  a             <- cjs_data$a  # Survival
  
  # Recapture Prob and Survival Prob -------------------------------------------------
  PF   <- rbeta(1,1,1)
  PM   <- rbeta(1,1,1)
  PhiF <- rbeta(1,1,1)
  PhiM <- rbeta(1,1,1)
  phi  <- PhiF * female + PhiM * (1-female) 
  p    <- PF   * female + PM   * (1-female)
  
  # Sample Survival
  for(i in 1:n){
    for(t in initial_entry[i]:(k-1)){
      if(is.na(a[i, t+1])){
        a[i, t+1] <- rbinom(1, 1, phi[i] * a[i, t])
      }
    }
  }
  
  # Add unknown status 
  build_NA_mat <- function(mat, cjs_mat){
    mat_final <- matrix(NA,nrow = dim(mat)[1], ncol = dim(mat)[2])
    mat_final[is.na(cjs_mat)] <- mat[is.na(cjs_mat)]
    return(mat_final)
  }
  
   # Survival Init
  a <- build_NA_mat(a, cjs_data$a)
  
  # Return Results ------------------------------------------------------------------
  
  # Store in object
  cjs_inits <- list(
    PF      = PF,
    PM      = PM,
    PhiF    = PhiF,
    PhiM    = PhiM,
    a       = a,
    p       = p,
    phi     = phi
  )
  
  # Return Initial Values for a single chain
  return(cjs_inits)
  
}

# Compile Model
compile_cjs_nimble <- function(cjs_data,
                               params = NULL){
  
  
  # Generating Initial Values
  cat("Generating Initial Values...", "\n")
  nimble_inits <- generate_init_cjs(cjs_data)
  
  # Construct Nimble Objects 
  cat("Organizing Data for Nimble...", "\n")
  
  nimble_cjs_constants <- list(
    n             = cjs_data$n,
    k             = cjs_data$k,
    female        = cjs_data$female,
    initial_entry = cjs_data$initial_entry
  )
  
  
  nimble_cjs_dat <- list(
    a       = cjs_data$a,
    x       = cjs_data$x
  )
  
  if(!is.null(params)){
    cat("User-specified Params Detected...","\n")
    cat("Using params := ", "\n")
    cat(params, "\n")
    
    nimble_params <- params
  } else {
    cat("Params argument is NULL...","\n")
    nimble_params <- c("PF","PM","PhiF","PhiM")
    cat("Using params := ", "\n")
    cat(nimble_params, "\n")
  }
  
  nimble_dims <- list(p       = c(nimble_cjs_constants$n),
                      phi     = c(nimble_cjs_constants$n))
  
  cat("Building Model Nodes in Nimble (SLOW)...", "\n")
  cjsModel <- nimbleModel(code       = nimble_cjs_model, 
                          constants  = nimble_cjs_constants, 
                          inits      = nimble_inits,
                          data       = nimble_cjs_dat,
                          dimensions = nimble_dims)
  
  
  # jsModel$simulate()
  lp_init <- cjsModel$calculate()
  print(paste0("LP from initial values is ", round(lp_init,3)))
  
  cat("Compiling Graphical Model in C++ (SLOW)...", "\n")
  compile_cjs <- compileNimble(cjsModel, showCompilerOutput = F)
  
  # [Note] SafeDepare.... warnings are annoying so suppress messages 
  cat("Configuring Markov Chain Monte Carlo Process (SLOW)...", "\n")
  cjsConf  <- suppressMessages(configureMCMC(cjsModel,
                                             print = F,
                                             useConjugacy = F,
                                             onlySlice = F))
  
  print(cjsConf)
  
  cat("Adding Monitors and Constructing MCMC...", "\n")
  cjsConf$addMonitors(nimble_params)
  cjsMCMC  <- buildMCMC(cjsConf)
  
  cat("Compiling MCMC Samplers (SLOW)...", "\n")
  CmdlMCMC <- compileNimble(cjsMCMC, project = cjsModel)
  
  cat("Project Defined, MCMC and Model are compiled...", "\n")
  
  cat("Returning Model Object...", "\n")
  
  return(list(CmdlMCMC = CmdlMCMC, 
              nimble_inits = nimble_inits))
  
  
}


# Get Samples from Model
run_nimble_cjs <- function(CmdlMCMC,
                           niter,
                           nburnin,
                           thin,
                           inits,
                           nchains=3,
                           seed = F){
  
  cat("MCMC Sampling from Model...","\n")
  samples <- runMCMC(mcmc              = CmdlMCMC,
                     niter             = niter,
                     nburnin           = nburnin,
                     thin              = thin,
                     inits             = inits,
                     nchains           = nchains,
                     setSeed           = seed,
                     samplesAsCodaMCMC = TRUE)
  
  cat("Returning Output...","\n")
  return(samples)
}

# Run compilation/initialization/sampling in one call 
# Dummy function intended for parallelization but can be used in serial as well
execute_cjs_nimble_pipeline <- function(seed, 
                                        data, 
                                        params, 
                                        niter, 
                                        nthin,
                                        nburnin, 
                                        nchains){
  
  nimble_complied <- compile_cjs_nimble(data, params)
  nim_inits <- lapply(1:nchains, function(i) generate_init_cjs(data))
  samples <- run_nimble_cjs(CmdlMCMC = nimble_complied$CmdlMCMC,
                            niter    = niter,
                            thin     = nthin,
                            nburnin  = nburnin,
                            nchains  = nchains,
                            inits    = nim_inits,
                            seed     = seed) 
  return(list(samples = samples,
              inits   = nim_inits))
  
}

run_cjs_nimble_parallel <- function(data, params, niter, nthin, nburnin, ncores){
  
  # Build cluster
  ncores <- min(ncores, detectCores()-2)
  
  if(ncores == 1){
    cat("Only one core specified, not bothering with parallelization...", "\n")
    cat("Set ncores > 1 if so desired...","\n")
    samples <- execute_cjs_nimble_pipeline(F, data, params, niter, nthin, nburnin, 1)
    return(samples)
  }
  
  cat(paste0("Building cluster with ", ncores , " sockets...", "\n"))
  cl <- makeCluster(ncores)
  
  # Load packages necessary
  cat(paste0("Loading custom functions onto cluster...", "\n"))
  clusterEvalQ(cl, {
    #  Load Libraries
    libs <- c("nimble", "coda", "tidyverse")
    lapply(libs,require, character.only = T)
    source(paste0(getwd(), "/Scripts/fn_cormack_jolly_seber_mod_nimble.R"))
    `%+%` <- function(a, b) paste0(a, b)
  })
  
  seeds <-  1:ncores
  
  cat(paste0("Running MCMC in Parallel (SLOW) ...", "\n"))
  out_list <- parLapply(cl      = cl,
                        X       = seeds, 
                        fun     = execute_cjs_nimble_pipeline, 
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
  samples <- as.mcmc.list(lapply(1:ncores, function(x) out_list[[x]]$samples))
  inits   <- lapply(1:ncores, function(x) out_list[[x]]$inits)
  
  return(list(samples = samples,
              inits   = inits,
              seed    = seeds))
  
}