#### Standard CJS Model 

nimble_cjs <- nimbleCode({
  # Multiple-Partners Joint-Fates Extended CJS Model-----------------------------------------------------------------------------
  #1. Data Augmentation----------------------------------------------------------------------------------------------------------------------
  
  # Nothing yet
  
  
  # 2. Recruitment into population-----------------------------------------------------------------------------------------------------------
  
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
    
    # 3. Decision to Mate -------------------------------------------------------------------------------------------------------------------
    
    # Female Mating Choice at time t
    for(i in 1:nf){
      amating_f[i,t] ~ dbern(af[i,t] * recruit_f[i,t] * delta)
    }
    
    # Male Mating Choice at time t
    for(j in 1:nm){
      amating_m[j,t] ~ dbern(am[j,t] * recruit_m[j,t] * delta)
    }
    
    # Choose to re-form pairs
    for(i in 1:nf){
      prob_repartner[i,t] <- ilogit(beta1*histories[i, apairs_f[i,t] , t]) * psi[i, apairs_f[i,t], t]
      arepartner[i,t] ~ dbern(prob_repartner[i,t] * amating_f[i,t] * amating_m[apairs_f[i,t],t])
    }
    
    # 4. Mate Selection -------------------------------------------------------------------------------------------------------------------
    # Use Categorical Distribution to classify mates
    
    # Build Homogeneous Partnership probabilities
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
    
    # Initialize choice selection
    psi_cond2[1, 1:(nm+1), t] <- c(psi_cond[1,1:nm,t],equals(sum(psi_cond[1,1:nm,t]),0))
    apairs_f[1,t+1] ~ dcat(psi_cond2[1, 1:(nm+1), t])
    single_female[1,t] <- equals(apairs_f[1,t+1],nm+1)
    
    # Attempts at partnerships forming
    # Monogamous pairings only
    for(i in 2:nf){
      
      # Remove Formed Pairs
      for(j in 1:nm){
        psi_cond2[i,j,t] <- psi_cond[i,j,t]*(1-sum(equals(apairs_f[1:(i-1),t+1],j)))
      }
      
      #Add case in which no pairs are available
      psi_cond2[i,(nm+1),t] <- equals(sum(psi_cond2[i,1:nm,t]),0)
      
      # Find mate for i
      apairs_f[i,t+1] ~ dcat(psi_cond2[i, 1:(nm+1), t])
      
      # Designate as single if no males available or if choosing not to mate
      single_female[i,t] <- equals(apairs_f[i,t+1],nm+1)
    }
    
    # Update Total History for Next Time Step
    for(i in  1:nf){
      for(j in  1:(nm+1)){
        histories[i, j, t+1] <- histories[i, j, t] + equals(apairs_f[i,t+1],j)*(1-single_female[i,t])
      }
    }
    
    # 5. Joint Survival ---------------------------------------------------------------------------------------------------------------------
    
    
    # Marginal Survival Event for Males in the Population (P[Y^M_T])
    for(j in 1:nm){
      am[j,t+1] ~ dbern(PhiM * am[j,t] * recruit_m[j,t])
    }
    
    # Marginal Recapture Event for Males in the Population (P[X^M_T|X^F_T])
    for(i in 1:nf){
      
      # Probability of female surviving given partnership and partner recapture status
      phi.totalF[i, t] <- single_female[i,t] * PhiF + # female was single
        (1 - single_female[i,t]) * (am[apairs_f[i,t+1],t+1] * (Phifm/PhiM) + # Male mated and female surived
                                      (1 - am[apairs_f[i,t+1],t+1]) * (Phif0/(1-PhiM))) # Male mated and female perished
      
      # Draw Survival Event
      af[i, t+1] ~ dbern(phi.totalF[i,t] * af[i,t] * recruit_f[i,t])
    }
    
    # 6. Joint Recapture --------------------------------------------------------------------------------------------------------------------
    
    
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
  
  
  # 7. Prior Distributions-------------------------------------------------------------------------------------------------------------------
  
  # Recruitment
  for(t in 1:k){
    eps[t] ~ dbeta(1,1)
  }
  
  # Attempt to Mate
  delta ~ dbeta(3,2)
  
  # Pairs reforming
  ##beta0 ~ dnorm(0, 1/4)
  beta1 ~ dnorm(0, 1/4)
  
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
  gamma_raw ~ dbeta(3,3)
  
  ###Frechet-Hoeffding Bounds for Correlation
  
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
  PhiF ~ dbeta(3,3)
  PhiM ~ dbeta(3,3)
  
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
  rho_raw ~ dbeta(3,3)
  
  ###Frechet-Hoeffding Bounds for Correlation
  
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
  PF ~ dbeta(3,3)
  PM ~ dbeta(3,3)
})

# nimble_cjs_inits <- function(){
#   
#   #unpack data
#   n <- cjs_data$n
#   k <- cjs_data$k
#   female <- cjs_data$female
#   initial_entry <- cjs_data$initial_entry
#   x <- cjs_data$x
#   a <- cjs_data$a
#   
#   # Generate Inits
#   phiF = runif(1,0,1)
#   phiM = runif(1,0,1)
#   pF = runif(1,0,1)
#   pM = runif(1,0,1)
#   
#   #Construct derived inits (for some reason)
#   phi <- matrix(NA, nrow = n, ncol = k)
#   p <- matrix(NA, nrow = n, ncol = k)
#   
#   # Build derived initial Values from samples
#   for(i in 1:n){
#     for(t in (initial_entry[i]+1):k){
#       phi[i,t] <- female[i]*phiF + (1-female[i])*phiM 
#       p[i,t] <- female[i]*pF + (1-female[i])*pM  
#       
#       #If we have a missing value then populate with a random sample
#       if(is.na(a[i,t])){
#         a[i,t] <- rbinom(1, 1, phi[i,t]*a[i,t-1])
#       }
#       
#       # Apply weights
#       phi[i,t] <- phi[i,t]**a[i,t-1]
#       p[i,t] <- p[i,t]*a[i,t]
#     }
#   }
# 
#   phi[,1] <- phi[,2]
#   p[,1] <- p[,2]
#   
# 
#   
#   
#   nimble_inits <- list(a = a,
#                        phiF = phiF,
#                        phiM = phiM,
#                        pF = pF,
#                        pM = pM,
#                        phi = phi, 
#                        p = p)
#   
#   # Return initial values
#   return(nimble_inits)
# } 
# 


cjs_constants <- list(
  nf = jags_data$nf,
  nm = jags_data$nm,
  k = jags_data$k
)

cjs_dat <- list(
  recruit_f = jags_data$recruit_f,
  recruit_m = jags_data$recruit_m,
  amating_f = jags_data$amating_f,
  amating_m = jags_data$amating_m,
  psi = jags_data$psi, #
  af = jags_data$af,
  am = jags_data$am,
  apairs_f = jags_data$apairs_f,#
  arepartner = jags_data$arepartner, #
  recap_f = jags_data$recap_f,
  recap_m = jags_data$recap_m
)

nimble_params <- c("PF","PM","rho","PhiF","PhiM","gamma","delta","beta0","beta1", "eps")

dims <- list(histories = c(cjs_constants$nf, cjs_constants$nm+1, cjs_constants$k+1),
             prob_repartner = c(cjs_constants$nf, cjs_constants$k),
             psi_raw = c(cjs_constants$nf, cjs_constants$nm, cjs_constants$k),
             male_taken_jt = c(cjs_constants$nm, cjs_constants$k),
             psi_cond = c(cjs_constants$nf, cjs_constants$nm, cjs_constants$k),
             psi_cond2 = c(cjs_constants$nf, cjs_constants$nm+1, cjs_constants$k),
             single_female = c(cjs_constants$nf, cjs_constants$k),
             phi.totalF = c(cjs_constants$nf, cjs_constants$k),
             p.totalF = c(cjs_constants$nf, cjs_constants$k)
             )


cjsModel <- nimbleModel(code = nimble_cjs, 
                        constants = cjs_constants, 
                        #inits = nimble_cjs_inits(),
                        data = cjs_dat,
                        dimensions = dims)


compile_cjs <- compileNimble(cjsModel, showCompilerOutput = TRUE)


cjsConf <- configureMCMC(cjsModel, print = TRUE)
#cjsConf$addSampler()
cjsConf$addMonitors(nimble_params)
cjsMCMC <- buildMCMC(cjsConf)
CcjsMCMC <- compileNimble(cjsMCMC, project = cjsModel)


samples <- runMCMC(CcjsMCMC, niter = 1e3, nburnin = 1e3, thin = 1, setSeed = TRUE, samplesAsCodaMCMC = TRUE)


coda.samples <- as.mcmc(samples)
