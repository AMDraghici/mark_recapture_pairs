# 1. Simulating Datasets for study---------------------------------------------------------------------------------
#    Code runs using embarrassing parallel and should work on both Linux/Windows/MAC

#Generate Data
sim_cr_dat <- function(parameter_list, iterations, ncores = detectCores() - 1){
  
  #Assign Core count (will not use more than the system has minus one)
  cat("Initializing cluster for data generation...\n")
  numcores <- ncores
  
  #Set up Cluster
  cl <- makeCluster(numcores)
  
  # Upload Functions, Variables, and Load Libraries
  cat("Pushing objects to children...\n")
  
  #Get Script Path
  path2scripts <- getwd() %+% "/Scripts"
  
  #Export Variables to clusters
  export <- list(
    "parameter_list", "iterations","path2scripts"
  )
  
  clusterExport(cl, export, envir = environment())
  clusterEvalQ(cl, source(paste0(path2scripts,"/00_fn_sim_pair_data.R")))
  
  # Set Random Seeds
  clusterSetRNGStream(cl)
  
  # Dummy Function
  f <- function(i,parameter_list) {
    return(sim_dat(parameter_list))
  }
  
  #Simulate Data
  cat("Generating data...\n")
  
  cr_dat_list <- parLapply(cl, 1:iterations, f, parameter_list)
  
  cat("Data generation complete...\n")
  
  # Close Cluster
  stopCluster(cl)
  
  #Return Results
  return(cr_dat_list)
}

# 2. Fitting JAGS Model--------------------------------------------------------------------------------------------
# Code runs using embarrassing parallel and should work on both Linux/Windows/MAC

# Build initial values for jags to use
generate_init <- function(jags_data, debug = F){
  
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
  psi_raw <- array(NA, dim = c(nf, nm, k))   # Update Exclusion matrix PSI (ignore repartner structure at time 1 since no repartnering has occured yet)
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
      
      # Draw Survival Event
      if(is.na(af[i,t+1])){
        af[i, t+1] <- rbinom(1,1, phi.totalF[i,t] * af[i,t] * recruit_f[i,t] + (1-recruit_f[j,t]))
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
    am = am
  )
  
  # Data objects to ensure pairs are monogamous and inform latent states
  # Only returned to help error checking initialization 
  jags_debug <- list(
    histories = histories,
    single_female = single_female,
    phi.totalF = phi.totalF,
    male_taken_jt = male_taken_jt, 
    prob_repartner = prob_repartner, 
    psi_raw = psi_raw,
    psi_cond = psi_cond,
    psi_cond2 = psi_cond2
  )
  
  if(debug){
    # Note the extra list structure if doing debugging
    jags_init_out <- list(jags_inits = jags_inits,
                          jags_debug = jags_debug)
  } else {
    jags_init_out <- jags_inits
  }
  
  # Return Initial Values for a single chain
  return(jags_init_out)
}

#Process Jags in Parallel (multiple chains same data)
run_jags_parallel <- function(jags_data, 
                              jags_model, 
                              jags_params, 
                              par_settings, 
                              out_dir,
                              save=TRUE,
                              outname = NULL){
  
  #Make sure you don't run the chain for hours and then lose the data
  if(missing(out_dir)){
    if(save == TRUE){
      stop("Specify an output directory using out_dir or set save == FALSE")
    }
  }
  
  #Create Initial Values for each chain using custom initial function
  jags_init <- list()
  
  #First chain uses flat transition initialization then remaining are random
  for(i in 1:par_settings$n.chains){
    jags_init[[i]] <- generate_init(jags_data)
  }
  
  #Start timer
  timer <- proc.time()
  
  #Assign clusters and pass through dclone
  n.cores <- detectCores() - 1
  workers <- pmin(n.cores, par_settings$n.chains)
  
  cat("Setting up " %+% workers %+% " parallel workers ... \n")
  cl <- makePSOCKcluster(workers)
  tmp <- clusterEvalQ(cl, library(dclone))
  
  cat("Initializing graph nodes and adapting chains with " %+% par_settings$n.adapt %+% " iterations ... \n")
  #Construct graph nodes and initialize in parallel
  parJagsModel(cl = cl, 
               name = 'pair_swap', 
               file = jags_model, 
               data = jags_data,
               n.chains = par_settings$n.chains, 
               n.adapt = par_settings$n.adapt,
               inits = jags_init)
  
  #Burn-in each chain in parallel
  cat("Burning-in chains with "%+% par_settings$n.burn %+% " iterations ... \n")
  parUpdate(cl = cl, 
            object = 'pair_swap', 
            n.iter = par_settings$n.burn)
  
  #Sample from distribution in parallel
  cat("Sampling distribution with " %+% par_settings$n.iter %+% " iterations and thinning of " %+% par_settings$n.thin %+%" ... \n ")
  jags_samples <- parCodaSamples(cl = cl, 
                                 model = 'pair_swap', 
                                 variable.names = jags_params, 
                                 n.iter = par_settings$n.iter, 
                                 thin = par_settings$n.thin)
  
  #Release clusters
  stopCluster(cl)
  
  #Stop timer and report time
  time.taken <- proc.time() - timer
  cat("MCMC has finished - time statistics are: ...\n")
  print(time.taken)
  
  ## Save output    
  if(save==TRUE){
    cat("Saving samples and initial values to output directory...\n")
    #Time stamp
    tstamp <- strftime(Sys.time(),format="%Y-%m-%d_%H:%M")
    
    #If no name is given just use the timestamp
    if(is.null(outname)){
      outfile_name <- "/jags_samples_" %+% tstamp
      outinit_name <- "/jags_inits_" %+% tstamp
    } else {
      outfile_name <- "/jags_samples_" %+% outname
      outinit_name <- "/jags_inits_" %+% outname
    }
    
    #Save JAGS Samples
    outfile_samples <- out_dir %+% outfile_name %+% ".rds"
    saveRDS(jags_samples,outfile_samples)
    
    #Save Jags Initial Values
    outfile_inits <- out_dir %+% outinit_name %+% ".rds"
    saveRDS(jags_init,outfile_inits)
  }
  
  #Return MCMC Results
  return(jags_samples)
}

# Run single chain of jags without parallel and return results + initial values
run_jags <- function(jags_data,
                     jags_model,
                     jags_params, 
                     par_settings,
                     debug = F){
  
  # Generate initial values
  if(debug){
    jags_init_out <- generate_init(jags_data, debug)
    jags_inits <- jags_init_out[["jags_inits"]]
    cat("Saving initial values for debugging....\n")
    saveRDS(jags_init_out, getwd() %+% "/jags_init_out_debug.rds")
    saveRDS(jags_data, getwd() %+% "/jags_data_out_debug.rds")
  } else {
    jags_inits <- generate_init(jags_data, debug = F)
  }
  
  # Construct Model Object
  jags_model_graph <- jags.model(file = jags_model, 
                                 data = jags_data,
                                 inits = jags_inits,
                                 n.chains = 1, 
                                 n.adapt = par_settings$n.adapt)
  
  # Burn - in process 
  update(jags_model_graph, par_settings$n.burn)
  
  # Sample from the model 
  jags_samples <- coda.samples(model = jags_model_graph, 
                               variable.names = jags_params, 
                               n.iter = par_settings$n.iter, 
                               thin = par_settings$n.thin)
  
  # Return inits and samples
  return(list(
    jags_inits = jags_inits,
    jags_samples = jags_samples
  ))
}

# Run jags in parallel with each chain representing a different dataset
run_jags_simulation_parallel <- function(jags_data_list,     # List of datasets to run model over 
                                         jags_model,         # Same model to be used for all runs
                                         jags_params,        # Same parameters to be saved for all runs
                                         par_settings,       # Same parameter settings for all runs 
                                         out_dir,            # Store in same directory if save = T 
                                         save=TRUE,          # Save results? 
                                         outname = NULL,     # Base name for files
                                         ncores  = detectCores() - 1 # number of cores
){ 
  
  #Assign Core count (will not use more than the system has minus one)
  cat("Initializing cluster for data modelling with JAGS...\n")
  numcores <- ncores
  
  #Set up Cluster
  cl <- makeCluster(numcores)
  
  # Upload Functions, Variables, and Load Libraries
  cat("Pushing objects to children...\n")
  
  #Get Script Path
  path2scripts <- getwd() %+% "/Scripts"
  
  #Export Variables to clusters
  export <- list(
    "jags_data_list", "jags_model","jags_params", "par_settings", "path2scripts"
  )
  
  clusterExport(cl, export, envir = environment())
  clusterEvalQ(cl, source(paste0(path2scripts,"/01_fn_model_code.R")))
  clusterEvalQ(cl, {
    library(boot)
    library(rjags)
  })
  
  # Set Random Seeds
  clusterSetRNGStream(cl)
  
  # Dummy Function
  f <- function(x,
                jags_data_list,
                jags_model, 
                jags_params,
                par_settings) {
    return(run_jags(jags_data_list[[x]],
                    jags_model,
                    jags_params, 
                    par_settings,
                    debug = F))
  }
  
  #Simulate Data
  cat("Running Models...\n")
  
  jags_results_list <- parLapply(cl, 1:length(jags_data_list), f, jags_data_list, jags_model, jags_params, par_settings)
  
  cat("Model Fitting Complete complete...\n")
  
  # Close Cluster
  stopCluster(cl)
  
  ## Save output    
  if(save==TRUE){
    cat("Saving samples and initial values to output directory...\n")
    #Time stamp
    tstamp <- strftime(Sys.time(),format="%Y-%m-%d_%H:%M")
    
    #If no name is given just use the timestamp
    if(is.null(outname)){
      outfile_name <- "/jags_samples_list" %+% tstamp
    } else {
      outfile_name <- "/jags_samples_list" %+% outname
    }
    
    #Save JAGS Samples and initial values
    outfile_samples <- out_dir %+% outfile_name %+% ".rds"
    saveRDS(jags_results_list,outfile_samples)
  }
  
  # Return list of outputs 
  return(jags_results_list)
}

# 3. Simulation Study Wrappers ---------------------------------------------------------------------------------------------------------------------

#Grid Generator (for multiple grid simulation)
generate_grid <- function(grid_base,to_vary){
  
  #Initialize Variables
  parameter_grid <- list()
  param_names <- names(to_vary)
  grid_entries <- bind_rows(to_vary[1])
  
  #Unpack Grid Entries 
  if(length(to_vary) > 1){
    for(i in 2:length(param_names)){
      param <- param_names[i]
      param_entries <- bind_rows(to_vary[param])
      grid_entries <- merge(grid_entries,param_entries,by=NULL)
    }
  }
  
  #Create Parameter Grid
  for(j in 1:nrow(grid_entries)){
    for(l in 1:length(param_names)){
      param <- param_names[l]
      dimension <- length(grid_base[param][[1]])
      grid_base[[param]] <- rep(as.vector(t(grid_entries[j,][param])),dimension)
    }
    parameter_grid[[j]] <- grid_base
  }
  #Return Results
  return(list(parameter_grid = parameter_grid,
              grid_entries = grid_entries))
}


#Add covariance part for grid simulation within mathematical bounds
add_cov_grid <- function(to_vary,parameter_list,param,by=0.5){
  
  #Covariance Bounds
  grid_entries <- generate_grid(parameter_list,to_vary)[[2]]
  #Add Parameters
  if(param == "gamma"){
    #Gamma Grid
    gam_stats <- compute_jbin_param_cjs(grid_entries$phi.f,grid_entries$phi.m)
    gam_bounds <- cbind(gam_stats[[5]],gam_stats[[6]])
    gam_range <- cbind(max(gam_bounds[,1]),min(gam_bounds[,2]))
    gam <- seq(gam_range[1],gam_range[2],by=by) %>% round(3)
    to_vary[["gamma"]] <- gam
  } else if(param== "rho"){
    #Rho Grid
    rho_stats <- compute_jbin_param_cjs(grid_entries$p.f,grid_entries$p.m)
    rho_bounds <- cbind(rho_stats[[5]],rho_stats[[6]])
    rho_range <- cbind(max(rho_bounds[,1]),min(rho_bounds[,2]))
    vrho <- seq(rho_range[1],rho_range[2],by=by) %>% round(3)
    to_vary[["rho"]] <- vrho
  }
  #Return Results
  return(to_vary)
}


# Run Simulation Study 
run_simulation <- function(grid_data, 
                           jags_model, 
                           jags_params,
                           par_settings,
                           out_dir,
                           iterations, 
                           ncores = detectCores() - 1){
  
  # Extract Grid Objects
  parameter_grid <- grid_data[["parameter_grid"]]
  grid_entries <- grid_data[["grid_entries"]]
  
  # Assign Batches based on each parameter setting
  batches <- grid_entries %>% nrow()
  cat("Initializing Simulation by assigning #" %+% batches %+% " batches" %+% "...\n")
  
  # Run simulation study in batches
  for(j in 1:batches){
    
    cat("Executing Batch #" %+% j %+% " of " %+% batches %+% "...\n")
    
    #Report List Slice
    cat("Generating data for Batch #" %+% j %+% "...\n")  
    
    # Generate Grid of Data
    jags_data_list <-  sim_cr_dat(parameter_list = parameter_grid[[j]], 
                                  iterations =  iterations,
                                  ncores = ncores)
    
    # Identify batches when post-processing
    outname <- "_batch_" %+% j %+% "_of_" %+% batches
    
    #Store Datasets
    cat("Saving data from Batch #" %+% j %+% "....\n")
    outdata_name <- out_dir %+%  "/jags_data_list" %+% outname %+% ".rds"
    saveRDS(jags_data_list, outdata_name)
    
    #Model Data
    cat("Running JAGS with " %+% ncores %+% " cores for Batch #" %+% j %+% "...\n")
    
    jags_results_list <- run_jags_simulation_parallel(jags_data_list = jags_data_list, 
                                                      jags_model = jags_model, 
                                                      jags_params = jags_params,
                                                      par_settings = par_settings,
                                                      out_dir = out_dir, 
                                                      save = T,
                                                      outname = outname,
                                                      ncores = ncores)
    
    #Report Completion 
    cat("Results stored for Batch #" %+% j %+% "....\n")
    
    # Create compact summary tables
    cat("Constructing summary table for Batch #" %+% j %+% "....\n")
    jags_summary_list <- lapply(1:iterations,
                                function(x) gather_posterior_summary(jags_results_list[[x]]$jags_samples))
    
    
    cat("Saving summary table for Batch #" %+% j %+% "....\n")
    outsumm_name <- out_dir %+%  "/jags_summary_list" %+% outname %+% ".rds"
    saveRDS(jags_summary_list, outsumm_name)
    
    
    cat("Modelling batch #" %+% j %+% " of " %+% batches %+% " complete...\n")
    cat("Freeing Memory....\n")
    
    # Drop data objects
    rm(jags_data_list)
    rm(jags_results_list)
    rm(jags_summary_list)
    
    # Silent Garbage Collection (free mem)
    invisible(gc())
  }
  
  cat("Simulation complete...\n")
}


# 4. Investigate Results for a Single Run (live data for example)-------------------------------------------------------------------

# Grab Summary statistics 
gather_posterior_summary <- function(fit){
  
  # Summarize Results
  summ <- summary(fit)
  
  # Put results into dataframe
  summ_stats <- as.data.frame(summ$statistics) %>% rownames_to_column(var = "Parameter")
  summ_quant <- as.data.frame(summ$quantiles) %>% rownames_to_column(var = "Parameter")
  post_stats <- inner_join(summ_stats, summ_quant, by = "Parameter") %>% 
    mutate(Parameter_Name = gsub("[[]","",gsub("[[0-9]]+","",inner_join(summ_stats, summ_quant, by = "Parameter")$Parameter))) # add unique par name
  
  # Return Results
  return(post_stats)
}

add_true_values <- function(post, param_list){
  
  true_vals <- c(param_list$p.f[1],
                 param_list$p.m[1],
                 param_list$phi.f[1],
                 param_list$phi.m[1],
                 param_list$delta[1],
                 # param_list$betas$beta0[1],
                 #param_list$betas$beta1[1],
                 c(1,rep(0,param_list$k-1)),
                 param_list$rho[1],
                 param_list$gam[1]
  )
  
  true_names <- c("PF","PM","PhiF","PhiM","delta",
                  "eps[" %+% 1:param_list$k %+% "]","rho","gamma")
  
  true_df <- data.frame("Parameter" = true_names, "true" = true_vals)
  
  # Join values together
  post2 <- post %>% inner_join(true_df, by = "Parameter")
  return(post2)
}

# Create Caterpillar Plot of Estimates
plot_caterpillar <- function(post_stats, 
                             params = c("PhiF","PhiM", "PF","PM"), 
                             yrange = NULL, 
                             title = "Caterpillar Plot of Posterior Estimates"
){
  p1 <- post_stats %>% 
    filter(Parameter_Name %in% params) %>% 
    mutate(Parameter = fct_reorder(Parameter,Mean)) %>% 
    ggplot() + 
    geom_linerange(aes(x = Parameter, ymax = `97.5%`, ymin = `2.5%`), alpha = 0.5, size = 1, color = "skyblue")  +
    geom_linerange(aes(x = Parameter, ymax = `25%`, ymin = `75%`), size = 2.0, alpha = 1, color = "lightblue") + 
    geom_point(aes(x = Parameter, y = Mean), size = 5, alpha = 0.75, color = "lightblue") +
    labs(x = "Parameter", y = "Posterior Values", title = title) +
    theme(plot.title = element_text(hjust = 0.5))
  
  if(!is.null(yrange)){
    p1 <- p1 +  coord_cartesian(ylim = yrange)
  }
  
  return(p1)
}

# Simulation Study Functions-----------------------------------------------------------------------------------------------------------------

# 1. Shuffle female order experiment 
shuffle_f_data <- function(jags_data){
  
  # Objects that need to be reordered/used in reordering
  recruit_f <- jags_data$recruit_f
  amating_f <- jags_data$amating_f
  apairs_f <- jags_data$apairs_f
  arepartner <- jags_data$arepartner
  psi <- jags_data$psi
  af <- jags_data$af
  recap_f <- jags_data$recap_f
  nf <- jags_data$nf 
  
  # Sample new order for females (final index is a dummy spot and is unchanged)
  new_index <- c(sample(1:nf, size = nf, replace = FALSE),nf+1)
  
  #  Reorder objects and Return shuffled data set
  jags_shuffled <- list(
    nf                =   jags_data$nf,
    nm                =   jags_data$nm,
    k                 =   jags_data$k, 
    recruit_m         =   jags_data$recruit_m,
    amating_m         =   jags_data$amating_m,
    am                =   jags_data$am, 
    recap_m           =   jags_data$recap_m,
    recruit_f         =   recruit_f[new_index,],
    amating_f         =   amating_f[new_index,],
    apairs_f          =   apairs_f[new_index[-(nf+1)],], 
    arepartner        =   arepartner[new_index[-(nf+1)],], 
    psi               =   psi[new_index[-(nf+1)],,],
    af                =   af[new_index,], 
    recap_f           =   recap_f[new_index,]
  )
  
  return(jags_shuffled)
}

replicate_shuffled_data <- function(jags_data, nsamples){
  jags_data_list <- lapply(1:nsamples, function(x) shuffle_data(jags_data))
  return(jags_data_list)
}

