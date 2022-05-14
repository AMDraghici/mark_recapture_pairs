# Nimble Functions for Running Jolly Seber Model ------------------------------------------------------------------------------------------------------------ 

# BUGS/JAGS Code
nimble_js_model <- nimbleCode({
  
  # Data Augmentation --------------------------------------------------------------------------------
  for(i in 1:n){
    z[i] ~ dbern(xi)
  }
  
  ## Compute population size
  N <- sum(z[1:n])
  
  # Recruit Likelihood -------------------------------------------------------------------------------
  for(i in 1:n){
    recruit[i, 1] ~  dbern(eps[1])
    for(t in 2:(k-1)){
      recruit[i,t] ~ dbern(recruit[i, t-1] + (1-recruit[i, t-1]) * eps[t])
    } 
  }
  
  # JS Likelihood -----------------------------------------------------------------------------------
  for(t in 2:k){
    for(i in 1:n){
      ## 1) Survival
      a[i,t] ~ dbern(phi[i] * a[i,t-1] * recruit[i,t]+ (1 - recruit[i,t]))
    }
  }
  
  for(t in 1:k){
    for(i in 1:n){
      x[i,t] ~ dbern(p[i] * a[i,t] * recruit[i,t] * z[i])
    }
  }
  
  
  # Priors--------------------------------------------------------------------------------------------
  
  # Split probs by sex 
  for(i in 1:n){
    p[i] <- female[i] * PF + (1-female[i]) * PM
    phi[i] <- female[i] * PhiF + (1-female[i])* PhiM
  }
  
  # # Data augmentation
  xi ~ dbeta(0.1,1.0)
  
  # Survival by sex
  PhiF ~ dbeta(1,1)
  PhiM ~ dbeta(1,1)
  
  # Recapture by sex
  PF ~ dbeta(1,1)
  PM ~ dbeta(1,1)
  
  # Recruitment 
  for(t in 1:k){
    eps[t] ~ dbeta(1,1)
  }
})


generate_init_js <- function(jags_data){
  
  #Unpack Variables -----------------------------------------------------------------
  
  # Known data and indices 
  k <- jags_data$k # Number of capture occasions
  n <- jags_data$n # Number of animals 
  female <- jags_data$female # Sex
  x <- jags_data$x # Recapture 
  z <- jags_data$z
  # CR data with missing components
  a <- jags_data$a  # Survival
  recruit <- jags_data$recruit # Recruitment
  
  # Recapture Prob and Survival Prob -------------------------------------------------
  PF <- rbeta(1,1,1)
  PM <- rbeta(1,1,1)
  PhiF <- rbeta(1,1,1)
  PhiM <- rbeta(1,1,1)
  eps <- rbeta(k, 1, 1)
  xi <- rbeta(1,1,0.1)
  
  # Sample augmentation
  for(i in 1:n){
    z[i] <- ifelse(is.na(z[i]), rbinom(1, 1, xi),z[i])
  }
  
  # Sample recruitment 
  # Recruitment 
  # Female Recruitment
  for(i in 1:n){
    recruit[i,1] <- ifelse(is.na(recruit[i,1]), rbinom(1, 1, eps[1]),recruit[i,1])
    for(t in 2:(k-1)){
      recruit[i,t] <- ifelse(is.na(recruit[i,t]), rbinom(1, 1, (recruit[i,t-1] + (1-recruit[i,t-1]) * eps[t])),recruit[i,t])
    } 
  }
  
  
  # Sample Survival
  for(i in 1:n){
    sexF <- female[i]
    phi <- PhiF * sexF + PhiM * (1-sexF) 
    for(t in 2:k){
      if(is.na(a[i, t])){
        a[i, t] <- rbinom(1, 1, phi * a[i, t-1] * recruit[i,t] + (1-recruit[i,t]))
      }
    }
  }
  
  # Add unknown status 
  build_NA_mat <- function(mat, jags_mat){
    mat_final <- matrix(NA,nrow = dim(mat)[1], ncol = dim(mat)[2])
    mat_final[is.na(jags_mat)] <- mat[is.na(jags_mat)]
    return(mat_final)
  }
  
  build_NA_vec <- function(vec, jags_vec){
    vec_final <- rep(NA, length(jags_vec))
    vec_final[is.na(jags_vec)] <- vec[is.na(jags_vec)]
    return(vec_final)
  }
  # Recruit init
  recruit <- build_NA_mat(recruit, jags_data$recruit)
  
  # Survival Init
  a <- build_NA_mat(a, jags_data$a)
  
  z <- build_NA_vec(z, jags_data$z)
  # Return Results ------------------------------------------------------------------
  
  # Store in object
  jags_inits <- list(
    PF      = PF,
    PM      = PM,
    PhiF    = PhiF,
    PhiM    = PhiM,
    eps     = eps,
    recruit = recruit,
    a       = a,
    z       = z,
    p       = PF * female + (1-female) * PM,
    phi     = PhiF * female + (1-female) * PhiM,
    xi      = xi
  )
  
  # Return Initial Values for a single chain
  return(jags_inits)
  
}

# Compile Model
compile_jolly_seber_nimble <- function(jags_data,
                                       params = NULL){
  
 
  # Generating Initial Values
  cat("Generating Initial Values...", "\n")
  nimble_inits <- generate_init_js(jags_data)
  
  # Construct Nimble Objects 
  cat("Organizing Data for Nimble...", "\n")
  
  nimble_js_constants <- list(
    n      = jags_data$n,
    k      = jags_data$k,
    female = jags_data$female
  )
  
  
  nimble_js_dat <- list(
    z       = jags_data$z,
    recruit = jags_data$recruit,
    a       = jags_data$a,
    x       = jags_data$x
  )
  
  if(!is.null(params)){
    cat("User-specified Params Detected...","\n")
    cat("Using params := ", "\n")
    cat(params, "\n")
    
    nimble_params <- params
  } else {
    cat("Params argument is NULL...","\n")
    nimble_params <- c("PF","PM","PhiF","PhiM", "eps","xi")
    cat("Using params := ", "\n")
    cat(nimble_params, "\n")
  }
  
  nimble_dims <- list(p       = c(nimble_js_constants$n),
                      phi     = c(nimble_js_constants$n))
  
  cat("Building Model Nodes in Nimble (SLOW)...", "\n")
  jsModel <- nimbleModel(code       = nimble_js_model, 
                         constants  = nimble_js_constants, 
                         inits      = nimble_inits,
                         data       = nimble_js_dat,
                         dimensions = nimble_dims)
  
  
  # jsModel$simulate()
  lp_init <- jsModel$calculate()
  print(paste0("LP from initial values is ", round(lp_init,3)))
  
  cat("Compiling Graphical Model in C++ (SLOW)...", "\n")
  compile_js <- compileNimble(jsModel, showCompilerOutput = F)
  
  node_names <- jsModel$getNodeNames()[jsModel$getNodeType(jsModel$getNodeNames())=="stoch"]
  
  # [Note] SafeDepare.... warnings are annoying so suppress messages 
  cat("Configuring Markov Chain Monte Carlo Process (SLOW)...", "\n")
  jsConf  <- suppressMessages(configureMCMC(jsModel,
                                            print = F,
                                            nodes = node_names,
                                            multivariateNodesAsScalars = T, 
                                            onlySlice = F))
  
  jsConf$removeSampler("x",print = F)
 
  print(jsConf)
  
  cat("Adding Monitors and Constructing MCMC...", "\n")
  jsConf$addMonitors(nimble_params)
  jsMCMC  <- buildMCMC(jsConf)
  
  cat("Compiling MCMC Samplers (SLOW)...", "\n")
  CjsMCMC <- compileNimble(jsMCMC, project = jsModel)
  
  cat("Project Defined, MCMC and Model are compiled...", "\n")
  
  cat("Returning Model Object...", "\n")
  
  return(list(CjsMCMC = CjsMCMC, 
              nimble_inits = nimble_inits))
  
  
}

# Get Samples from Model
run_nimble <- function(CmdlMCMC, 
                       niter,
                       nburnin,
                       thin,
                       seed = F){
  
  cat("MCMC Sampling from Model...","\n")
  samples <- runMCMC(CmdlMCMC,
                     niter             = niter,
                     nburnin           = nburnin, 
                     thin              = thin, 
                     setSeed           = seed,
                     samplesAsCodaMCMC = TRUE)
  
  cat("Converting to CODA samples...","\n")
  coda.samples <- as.mcmc(samples)
  cat("Returning Output...","\n")
  return(coda.samples)
}
