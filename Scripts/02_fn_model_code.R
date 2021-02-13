# 1. Simulating Datasets for study---------------------------------------------------------------------------------
#    Code runs using embarrassing parallel and should work on both Linux/Windows/MAC

#Generate Data
sim_cr_dat <- function(parameter_list,iterations, ncores = detectCores() - 1){
  
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
#    Code runs using embarrassing parallel and should work on both Linux/Windows/MAC



#Process Jags in Parallel
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
  #jags_init <- list()
  
  #First chain uses flat transition initialization then remaining are random
  # for(i in 1:par_settings$n.chains){
  #   jags_init[[i]] <- generate_init(jags_data,psi = psi_setting)
  # }
  # 
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
               n.adapt = par_settings$n.adapt)#,
               #inits = jags_init)
  
  #Burn-in each chain in parallel
  cat("Burning-in chains with "%+% par_settings$n.burn %+% " iterations ... \n")
  parUpdate(cl = cl, 
            object = 'pair_swap', 
            n.iter = par_settings$n.burn)
  
  #Sample from distribution in parallel
  cat("Sampling distribution with " %+% par_settings$n.iter %+% " iterations with thinning of " %+% par_settings$n.thin %+%" ... \n ")
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
