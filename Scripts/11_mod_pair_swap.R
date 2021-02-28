#### Multiple-Partners Joint-Fates Extended CJS Model 

model{
  
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
  
  # Conditional Partnership/Survival/Recapture Steps ----------------------------------------------------------------------------------------
  
  # Time 1 events ---------------------------------------------------------------------------------------------------------------------------
  
  # Mating at Time 1 
  
  # Female Mating at time 1
  for(i in 1:nf){
    amating_f[i,1] ~ dbern(delta*recruit_f[i,1])
  }
  
  # Male mating at time 1
  for(j in 1:nm){
    amating_m[j,1] ~ dbern(delta*recruit_m[j,1])
  }
  
  # Pair - Formation at Time 1
  
  # Build Partnership probabilities (unconditioned) 
  
  # Link function for individual components 
  for(i in 1:nf){
    for(j in 1:nm){
      eta[i,j, 1] <- beta0 # Linear function
      psi_raw[i, j ,1] <- exp(eta[i, j, 1])
    }
    
    # Normalize Multivariate Logit Function (how does jags not have this)
    psi[i,1:nm, 1] <- psi_raw[i,1:nm, 1]/sum(psi_raw[i, 1:nm , 1])
  }
  
  # Compute conditional psi (on recruitment and mating status)
  for(i in 1:nf){
    for(j in 1:nm){
      psi_cond[i,j,1] <- psi[i,j,1] * amating_f[i,1] * amating_m[j,1]
    }
  }
  
  # First attempted partnership
  apairs[1,1,1] ~ dbern(psi_cond[1,1,1])
  
  # Remaining attempts for female 1 
  for(j in 2:nm){
    apairs[1,j,1] ~  dbern(psi_cond[1,j,1] * (1 - sum(apairs[1, 1:(j-1), 1])))
  }
  
  # Attempts for females 2-nf 
  for(i in 2:nf){
    # First try
    apairs[i, 1, 1] ~  dbern(psi_cond[i,1,1]* (1 - sum(apairs[1:(i-1), 1, 1])))
    for(j in 2:nm){
      apairs[i,j,1] ~ dbern(psi_cond[i,j,1] * (1 - sum(apairs[i, 1:(j-1), 1])) * (1 - sum(apairs[1:(i-1), j, 1])))
    }
  }
  
  # If no pairs assigned then assign single to boundaries 
  for(i in 1:nf){ # Females
    single_female[i, 1] <- 1 - sum(apairs[i,1:nm,1]) # if no pairs then classify as single
  }
  
  for(j in 1:nm){ # Males
    single_male[j, 1] <- 1 - sum(apairs[1:nf,j,1]) # if no pairs then classify as single
  }
  
  # Survival at time 1 is ensured - therefore value is coded as 1 in the data for all individuals
  
  # Recapture at time 1 
  
  
  # Recapture Probability for pairs
  for(i in 1:nf){
    for(j in 1:nm){
      # Marginal probability of recapture for females
      p.totalF[i,j,1] <- apairs[i,j,1]*PF + (1 - apairs[i,j,1])*recap_f[i,1]
      # P[X^F_1]
      rfmat[i,j,1] ~ dbern(p.totalF[i,j,1])
      # Conditional probability of recapture for males
      p.totalM_F[i,j,1] <- apairs[i,j,1]*(rfmat[i,j,1]*(Pfm/PF) + (1 - rfmat[i,j,1])*(Pm0/(1-PF))) + 
        (1 - apairs[i,j,1])*recap_m[j,1]
      #P[X^M_T|X^F_1]
      rmmat[i,j,1] ~ dbern(p.totalM_F[i,j,1])        
    }
  }
  
  # Recapture probability for single Males
  for(j in 1:nm){
    # Probability of recapture for males
    p.totalM_F[nf+1,j,1] <- single_male[j,1]*PM + (1 - single_male[j,1])*recap_m[j,1]
    
    #P[Y^M_1]
    rmmat[nf+1,j,1] ~ dbern(p.totalM_F[nf+1,j,1] * recruit_m[j,1]) 
  }
  
  
  # Recapture Probability for single Females
  for(i in 1:nf){
    # Marginal probability of recapture for females
    p.totalF[i, nm+1, 1] <- single_female[i,1] * PF + (1 - single_female[i,1])*recap_f[i,1]
    
    # P[X^F_1]
    rfmat[i, nm+1, 1] ~ dbern(p.totalF[i,nm+1,1] * recruit_f[i,1])
  }
  
  
  # Model Events from t=2 to k --------------------------------------------------------------------------------------------------------------
  for(t in 2:k){
    
    # 3. Decision to Mate -------------------------------------------------------------------------------------------------------------------
    
    # Female Mating Choice at time t
    for(i in 1:nf){
      amating_f[i,t] ~ dbern(af[i,t-1] * recruit_f[i,t] * delta])
    }
    
    # Male Mating Choice at time t
    for(j in 1:nm){
      amating_m[j,t] ~ dbern(am[j,t-1] * recruit_m[j,t] * delta)
    } 
    
    
    # 4. Mate Selection -------------------------------------------------------------------------------------------------------------------
    # JAGs does not play well with multinomial trials
    # Solution - series of conditional binomial trials with constraints for singles
    # Probabilities have to be unconditioned after sampling
    
    # Build Partnership probabilities 
    
    # Recover total history 
    for(i in 1:nf){
      for(j in 1:nm){
        histories[i, j, t] <- sum(apairs[i, j, 1:(t-1)])
      }
    }
    
    # Linear function 
    eta[1:nf, 1:nm, t] <- beta0 + beta1*histories[1:nf, 1:nm, t]
    
    # Link function for individual components 
    for(i in 1:nf){
      for(j in 1:nm){
        psi_raw[i, j, t] <- exp(eta[i, j, t])
      }
      
      # Normalize Multivariate Logit Function (how does jags not have this)
      psi[i,1:nm,t] <- psi_raw[i,1:nm,t]/sum(psi_raw[i, 1:nm ,t])
    }
    
    # Compute conditional psi (on recruitment and mating status)
    for(i in 1:nf){
      for(j in 1:nm){
        psi_cond[i,j,t] <- psi[i,j,t] * amating_f[i,t] * amating_m[j,t]
      }
    }
    
    # First attempted partnership
    apairs[1,1,t] ~ dbern(psi_cond[1,1,t])
    
    # Remaining attempts for female 1 
    for(j in 2:nm){
      apairs[1,j,t] ~  dbern(psi_cond[1,j,t] * (1 - sum(apairs[1, 1:(j-1), t])))
    }
    
    # Attempts for females 2-nf 
    for(i in 2:nf){
      # First try
      apairs[i, 1, t] ~  dbern(psi_cond[i,1,t]* (1 - sum(apairs[1:(i-1), 1, t])))
      for(j in 2:nm){
        apairs[i,j,t] ~ dbern(psi_cond[i,j,t] * (1 - sum(apairs[i, 1:(j-1), t])) * (1 - sum(apairs[1:(i-1), j, t])))
      }
    }
    
    # If no pairs assigned then assign single to boundaries 
    for(i in 1:nf){ # Females
      single_female[i, t] <- 1 - sum(apairs[i,1:nm,t]) # if no pairs then classify as single
    }
    
    for(j in 1:nm){ # Males
      single_male[j, t] <- 1 - sum(apairs[1:nf,j,t]) # if no pairs then classify as single
    }
    
    # 5. Joint Survival ---------------------------------------------------------------------------------------------------------------------
    
    # 5a. Survival Probability for pairs
    for(i in 1:nf){
      for(j in 1:nm){
        
        # Marginal probability of survival for females
        phi.totalF[i,j,t] <- apairs[i,j,t]*PhiF + (1 - apairs[i,j,t])
        
        # P[Y^F_T]
        afmat[i,j,t] ~ dbern(phi.totalF[i,j,t])
        
        # Conditional probability of survival for males given female survival 
        phi.totalM_F[i,j,t] <- apairs[i,j,t]*(afmat[i,j,t]*(Phifm/PhiF) + (1 - afmat[i,j,t])*(Phim0/(1-PhiF))) +(1 - apairs[i,j,t])
        
        #P[Y^M_T|Y^F_T]
        ammat[i,j,t] ~ dbern(phi.totalM_F[i,j,t])        
      }
    }
    
    
    # 5b. Survival probability for single Males
    for(j in 1:nm){
      # Probability of survival for males
      phi.totalM_F[nf+1,j,t] <- single_male[j,t]*recruit_m[j,t]*(PhiM - 1) + 1
      #P[Y^M_T]
      ammat[nf+1,j,t] ~ dbern(phi.totalM_F[nf+1,j,t]*am[j,t-1]) 
    }
    
    
    # 5c. Survival Probability for single Females
    for(i in 1:nf){
      # Probability of survival for females
      phi.totalF[i,nm+1,t] <- single_female[i,t]*recruit_f[i,t]*(PhiF - 1) +  1
      #P[Y^F_T]
      afmat[i,nm+1,t] ~ dbern(phi.totalF[i,nm+1,t]*af[i,t-1]) 
    }
    
    # Recover Marginal Survival Distribution for future conditioning
    # This is a small hack to get around the way JAGS defines data
    # Missing values are populated based on joint survival distribution
    
    #Females
    for(i in 1:nf){
      # Dummy parameter (its always 1 or zero)
      phif_dummy[i,t] <- sum(afmat[i,1:(nm+1),t]*c(apairs[i,1:nm,t],single_female[i,t]))
      # Missing values 
      af[i,t] ~ dbern(phif_dummy[i,t]) 
    }
    
    #Males
    for(j in 1:nm){
      # Dummy parameter (its always 1 or zero)
      phim_dummy[j,t] <- sum(ammat[1:(nf+1),j,t]*c(apairs[1:nf, j, t], single_male[j,t]))
      # Missing values 
      am[j,t] ~ dbern(phim_dummy[j,t]) 
    } 
    
    # 6. Joint Recapture --------------------------------------------------------------------------------------------------------------------
    
    # 6a. Recapture Probability for pairs
    for(i in 1:nf){
      for(j in 1:nm){
        # Marginal probability of recapture for females
        p.totalF[i,j,t] <- apairs[i,j,t]*PF + (1 - apairs[i,j,t])*recap_f[i,t]
        # P[X^F_T]
        rfmat[i,j,t] ~ dbern(p.totalF[i,j,t]*af[i,t])
        # Conditional probability of recapture for males
        p.totalM_F[i,j,t] <- apairs[i,j,t]*(rfmat[i,j,t]*(Pfm/PF) + (1 - rfmat[i,j,t])*(Pm0/(1-PF))) + (1 - apairs[i,j,t])*recap_m[j,t]
        #P[X^M_T|X^F_T]
        rmmat[i,j,t] ~ dbern(p.totalM_F[i,j,t]*am[j,t])        
      }
    }
    
    # 6b. Recapture probability for single Males
    for(j in 1:nm){
      # Probability of recapture for males
      p.totalM_F[nf+1,j,t] <- single_male[j,t]*PM + (1 - single_male[j,t])*recap_m[j,t]
      
      #P[Y^M_T]
      rmmat[nf+1,j,t] ~ dbern(p.totalM_F[nf+1,j,t]* am[j,t] * recruit_m[j,t]) 
    }
    
    
    # 6c. Recapture Probability for single Females
    for(i in 1:nf){
      # Marginal probability of recapture for females
      p.totalF[i, nm+1, t] <- single_female[i,t] * PF + (1 - single_female[i,t])*recap_f[i,t]
      
      # P[X^F_T]
      rfmat[i, nm+1, t] ~ dbern(p.totalF[i,nm+1,t]* af[i,t] * recruit_f[i,t])
    }
    
  }
  
  
  # 7. Prior Distributions-------------------------------------------------------------------------------------------------------------------
  
  #Prior for linear terms in pair selection
  beta0 ~ dnorm(0,1)
  beta1 ~ dnorm(0,1)
  
  # Recruitment 
  for(t in 1:k){
    eps[t] ~ dbeta(1,1)
  }
  
  # Attempt to Mate 
  delta ~ dbeta(1,1)
  
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
  gamma ~ dunif(gl, gu)
  
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
  
  ##Correlation (with FH bounds)
  rho ~ dunif(rl, ru)
  
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
  PF ~ dbeta(1,1)
  PM ~ dbeta(1,1)
}