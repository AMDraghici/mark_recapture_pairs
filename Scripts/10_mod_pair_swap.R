#### Multiple-Partners Joint-Fates Extended Arnason-Schwarz Model 

model{
  ############################
  ### 1. Data Augmentation ###
  ############################
  
  # Nothing yet 
  
  ##################################
  # 2. Recruitment into population #
  ##################################
  
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
  
  ########################################################
  ### Conditional Partnership/Survival/Recapture Steps ###
  ########################################################
  
  # Time 1 events
  
  # Female Mating at time 1
  for(i in 1:nf){
    amating_f[i,1] ~ dbern(delta[1]*recruit_f[i,1])
  }
  
  # Male mating at time 1
  for(j in 1:nm){
    amating_m[j,1] ~ dbern(delta[1]*recruit_m[j,1])
  }
  
  # Model Joint Partnership/Survival/Recapture outcomes 
  for(t in 2:k){
    
    #######################
    # 3. Decision to Mate #
    #######################
    
    # Female Mating Choice at time t
    for(i in 1:nf){
      amating_f[i,t] ~ dbern(af[i,t-1] * recruit_f[i,t] * delta[t])
    }
    
    # Male Mating Choice at time t
    for(j in 1:nm){
      amating_m[j,t] ~ dbern(am[j,t-1] * recruit_m[j,t] * delta[t])
    } 
    
    #####################
    # 4. Mate Selection #
    #####################
    # JAGs does not play well with multinomial trials
    # Solution - series of conditional binomial trials with constraints for singles
    # Probabilities have to be unconditioned after sampling
    
    
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
    
    #####################
    # 5. Joint Survival #
    #####################
    
    # 5a. Survival Probability for pairs
    for(i in 1:nf){
      for(j in 1:nm){
        
        # Marginal probability of survival for females #apairs[i,j,t]
        phi.totalF[i,j,t] <- apairs[i,j,t]*PhiF[t-1] + (1 - apairs[i,j,t])
        
        # P[Y^F_T]
        afmat[i,j,t] ~ dbern(phi.totalF[i,j,t]*af[i,t-1])
        
        # Conditional probability of survival for males given female survival 
        phi.totalM_F[i,j,t] <- apairs[i,j,t]*(afmat[i,j,t]*(Phifm[t-1]/PhiF[t-1]) + (1 - afmat[i,j,t])*(Phim0[t-1]/(1-PhiF[t-1]))) + (1 - apairs[i,j,t])
        
        #P[Y^M_T|Y^F_T]
        ammat[i,j,t] ~ dbern(phi.totalM_F[i,j,t]*am[j,t-1])        
      }
    }
    
    
    # 5b. Survival probability for single Males
    for(j in 1:nm){
      # Probability of survival for males
      phi.totalM_F[nf+1,j,t] <- single_male[j,t]*(PhiM[t-1]) + (1 - single_male[j,t])
      #P[Y^M_T]
      ammat[nf+1,j,t] ~ dbern(phi.totalM_F[nf+1,j,t]*am[j,t-1]) 
    }
    
    
    # 5c. Survival Probability for single Females
    for(i in 1:nf){
      # Probability of survival for females
      phi.totalF[i,nm+1,t] <- single_female[i,t]*(PhiF[t-1]) +  (1 - single_female[i,t])
      #P[Y^F_T]
      afmat[i,nm+1,t] ~ dbern(phi.totalF[i,nm+1,t]*af[i,t-1]) 
    }
    
    # Recover Marginal Survival Distribution for future conditioning
    
    #Females
    for(i in 1:nf){
      # Dummy parameter (its always 1 or zero)
      phif_dummy[i,t] <- sum(afmat[i,,t]*c(apairs[i,,t],single_female[i,t]))
      # Missing values are populated based on asurv (hack to get around data/constant problem)
      af[i,t] ~ dbern(phif_dummy[i,t]) 
    }
    
    #Males
    for(j in 1:nm){
      # Dummy parameter (its always 1 or zero)
      phim_dummy[j,t] <- sum(ammat[,j,t]*c(apairs[,j,t], single_male[j,t]))
      # Missing values are populated based on asurv (hack to get around data/constant problem)
      am[j,t] ~ dbern(phim_dummy[j,t]) 
    } 
    
    ######################
    # 6. Joint Recapture #
    ######################
    
    # 6a. Recapture Probability for pairs
    for(i in 1:nf){
      for(j in 1:nm){
        # Marginal probability of recapture for females
        p.totalF[i,j,t] <- apairs[i,j,t]*PF[t] + (1 - apairs[i,j,t])*recap_f[i,t]
        # P[X^F_T]
        rfmat[i,j,t] ~ dbern(p.totalF[i,j,t]*af[i,t])
        # Conditional probability of recapture for males
        p.totalM_F[i,j,t] <- apairs[i,j,t]*(rfmat[i,j,t]*(Pfm[t]/PF[t]) + (1 - rfmat[i,j,t])*(Pm0[t]/(1-PF[t]))) + 
          (1 - apairs[i,j,t])*recap_m[j,t]
        #P[X^M_T|X^F_T]
        rmmat[i,j,t] ~ dbern(p.totalM_F[i,j,t]*am[j,t])        
      }
    }
    
    # 6b. Recapture probability for single Males
    for(j in 1:nm){
      # Probability of recapture for males
      p.totalM_F[nf+1,j,t] <- single_male[j,t]*PM[t] + (1 - single_male[j,t])*recap_m[j,t]
      
      #P[Y^M_T]
      rmmat[nf+1,j,t] ~ dbern(p.totalM_F[nf+1,j,t]*am[j,t]) 
    }
    
    
    # 6c. Recapture Probability for single Females
    for(i in 1:nf){
      # Marginal probability of recapture for females
      p.totalF[i,nm+1,t] <- single_female[i,t]*PF[t] + (1 - single_female[i,t])*recap_f[i,t]
      
      # P[X^F_T]
      rfmat[i,nm+1,t] ~ dbern(p.totalF[i,nm+1,t]*af[i,t])
    }
    
  }
  ###########################
  ### Prior distributions ###
  ###########################
  
  # Build Partnership probabilities (unconditioned)
  # !!!! Need to modify histories to be latent for now just assume its some random covariate
  for(t in 1:k){
    # Linear function 
    eta[1:nf, 1:nm, t] <- beta0 + beta1*histories[1:nf, 1:nm, t]
    
    # Link function for individual components 
    for(i in 1:nf){
      for(j in 1:nm){
        psi_raw[i, j ,t] <- exp(eta[i, j, t])
      }
      
      # Normalize Multivariate Logit Function (how does jags not have this)
      psi[i,1:nm,t] <- psi_raw[i,1:nm,t]/sum(psi_raw[i, 1:nm ,t])
    }
  }
  
  #prior for mating
  beta0 ~ dnorm(0,1)
  beta1 ~ dnorm(0,1)
  
  for(t in 1:k){
    
    ### Derived Parameters ####
    
    #Joint Survival probabilities for paired individuals
    Phi00[t] <- 1 - PhiF[t] - PhiM[t] + Phifm[t]
    Phif0[t] <- PhiF[t] - Phifm[t]
    Phim0[t] <- PhiM[t] - Phifm[t]
    Phifm[t] <- gamma[t]*sig.PhiF[t]*sig.PhiM[t] + PhiF[t]*PhiM[t]
    
    #Joint Capture probabilities for paired individuals
    P00[t] <- 1 - PF[t] - PM[t] + Pfm[t]
    Pf0[t] <- PF[t] - Pfm[t]
    Pm0[t] <- PM[t] - Pfm[t]
    Pfm[t] <- rho[t]*sig.PF[t]*sig.PM[t] + PF[t]*PM[t]
    
    ###Binomial SD for survival and recapture 
    sig.PhiF[t] <- sqrt(PhiF[t]*(1-PhiF[t]))
    sig.PhiM[t] <- sqrt(PhiM[t]*(1-PhiM[t]))
    sig.PF[t] <- sqrt(PF[t]*(1-PF[t]))
    sig.PM[t] <- sqrt(PM[t]*(1-PM[t]))
    
    ###Frechet-Hoeffding Bounds for Correlation
    
    # Survival Rates (Gamma)
    gu[t] <-  min(sqrt(OR.Phi[t]), 1/sqrt(OR.Phi[t])) 
    gl[t] <- -min(sqrt(OP.Phi[t]), 1/sqrt(OP.Phi[t])) 
    
    # Recapture Rates (Rho)
    ru[t] <-  min(sqrt(OR.Phi[t]), 1/sqrt(OR.Phi[t])) 
    rl[t] <- -min(sqrt(OP.Phi[t]), 1/sqrt(OP.Phi[t])) 
    
    # Odds Ratio and Product of Survival and Recapture Rates
    OP.Phi[t] <- odds.PhiF[t]*odds.phiM[t]
    OR.Phi[t] <- odds.PhiF[t]/odds.phiM[t]
    OP.P[t] <- odds.PF[t]*odds.PM[t]
    OR.P[t] <- odds.PF[t]/odds.PM[t]
    
    ### Odds of Survival and Recapture Rates
    odds.phiM[t] <- PhiM[t]/(1 - PhiM[t])
    odds.PhiF[t] <- PhiF[t]/(1 - PhiF[t])
    odds.PM[t] <- PM[t]/(1 - PM[t])
    odds.PF[t] <- PF[t]/(1 - PF[t])
    
    ### Prior Parameters ####
    
    ##Correlation S/R (with FH bounds)
    gamma[t] ~ dunif(gl[t], gu[t])
    rho[t] ~ dunif(rl[t], ru[t])
    
    ##Survival Rates M/F
    PhiF[t] ~ dbeta(1,1)
    PhiM[t] ~ dbeta(1,1)
    
    # Recapture Rates M/F
    PF[t] ~ dbeta(1,1)
    PM[t] ~ dbeta(1,1)
    
    # Recruitment 
    eps[t] ~ dbeta(1,1)
    
    # Mating
    delta[t] ~ dbeta(1,1)
    
  } 
}