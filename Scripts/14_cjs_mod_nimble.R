#### Standard CJS Model 

nimble_cjs <- nimbleCode({
  # CJS Likelihood -----------------------------------------------------------------------------------
  for(i in 1:n){
    for(t in (initial_entry[i]+1):k){
      
      x[i,t] ~ dbern(prob = p[i,t])
      a[i,t] ~ dbern(prob = phi[i,t])
     
      phi[i,t] <- (female[i]*phiF + (1 - female[i])*phiM)*a[i,t-1]
      p[i,t] <- (female[i]*pF + (1 - female[i])*pM)*a[i,t]
      
      
    }
  }
  
  # Priors--------------------------------------------------------------------------------------------
  
  # Survival by sex
  phiF ~ dbeta(1,1)
  phiM ~ dbeta(1,1)
  
  # Recapture by sex
  pF ~ dbeta(1,1)
  pM ~ dbeta(1,1)
})

nimble_cjs_inits <- function(){
  
  #unpack data
  n <- cjs_data$n
  k <- cjs_data$k
  female <- cjs_data$female
  initial_entry <- cjs_data$initial_entry
  x <- cjs_data$x
  a <- cjs_data$a
  
  # Generate Inits
  phiF = runif(1,0,1)
  phiM = runif(1,0,1)
  pF = runif(1,0,1)
  pM = runif(1,0,1)
  
  #Construct derived inits (for some reason)
  phi <- matrix(NA, nrow = n, ncol = k)
  p <- matrix(NA, nrow = n, ncol = k)
  
  # Build derived initial Values from samples
  for(i in 1:n){
    for(t in (initial_entry[i]+1):k){
      phi[i,t] <- female[i]*phiF + (1-female[i])*phiM 
      p[i,t] <- female[i]*pF + (1-female[i])*pM  
      
      #If we have a missing value then populate with a random sample
      if(is.na(a[i,t])){
        a[i,t] <- rbinom(1, 1, phi[i,t]*a[i,t-1])
      }
      
      # Apply weights
      phi[i,t] <- phi[i,t]**a[i,t-1]
      p[i,t] <- p[i,t]*a[i,t]
    }
  }

  phi[,1] <- phi[,2]
  p[,1] <- p[,2]
  

  
  
  nimble_inits <- list(a = a,
                       phiF = phiF,
                       phiM = phiM,
                       pF = pF,
                       pM = pM,
                       phi = phi, 
                       p = p)
  
  # Return initial values
  return(nimble_inits)
} 



cjs_constants <- list(
  n = cjs_data$n,
  k = cjs_data$k,
  female = cjs_data$female,
  initial_entry = cjs_data$initial_entry
)

cjs_dat <- list(
  a = cjs_data$a,
  x = cjs_data$x
)

nimble_params <- c("phiF", "phiM", "pF", "pM")

dims <- list(phi = c(cjs_data$n, cjs_data$k),
             p = c(cjs_data$n, cjs_data$k))


cjsModel <- nimbleModel(code = nimble_cjs, 
                        constants = cjs_constants, 
                        inits = nimble_cjs_inits(),
                        data = cjs_dat,
                        dimensions = dims)


compile_cjs <- compileNimble(cjsModel, showCompilerOutput = TRUE)


cjsConf <- configureMCMC(cjsModel, print = TRUE)
#cjsConf$addSampler()
cjsConf$addMonitors(c("phiF","phiM","pF","pM"))
cjsMCMC <- buildMCMC(cjsConf)
CcjsMCMC <- compileNimble(cjsMCMC, project = cjsModel)




samples <- runMCMC(CcjsMCMC, niter = 1e5, nburnin = 1e5/2, thin = 10, setSeed = TRUE,samplesAsCodaMCMC = TRUE)


coda.samples <- as.mcmc(samples)
