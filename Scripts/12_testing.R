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
    
    # Compute conditional psi (on recruitment and mating status)
    for(i in 1:nf){
      for(j in 1:nm){
        psi_cond[i,j,t] <- psi[i,j,t] * amating_f[i,t] * amating_m[j,t]
      }
    }
    
    # Quadrant 2 and 4
    for(i in (nf+1):n){
      for(j in 1:n){
        psi_cond[i,j,t] <- psi[i,j,t]
      }
    }
    
    #Quadrant 3
    for(i in 1:nf){
      for(j in (nm+1):n){
        psi_cond[i,j,t] <- psi[i,j,t]
      }
    }

    # Assign mate index
    for(i in 1:n){
      pf[i,t] ~ dcat(psi[i, 1:n, t])
    }


    
    #pf[1:n,t] ~ dsample(rep(1,n),n)
    #########################
    # Parameter Constraints #
    #########################
    
    #Partner choice cannot occur more than once at t (no polygamy)
   #is_chosen_sum[1:n,t] <- apairs[1:n,1:n,t] %*% t(rep(1,n))
   # for(i in 1:nm){
   #  is_chosen_sum[i, t] <- sum(pf[1:n,t] == i)
   #  zeroes_mat[i,t] ~ dinterval(is_chosen_sum[i,t], 1) # force each choice to occur between [0, 1] times.
   # }

    #####################
    # 5. Joint Survival #
    #####################
    
    
    # 5a. Survival Probability for apairs
    for(i in 1:nf){
      for(j in 1:nm){
        
        # Marginal probability of survival for females #apairs[i,j,t]
        phi.totalF[i,j,t] <- equals(pf[i,t],j)*PhiF[t-1] + 
          (1 - equals(pf[i,t],j))*af[i,t-1]
        
        # P[Y^F_T]
        afmat[i,j,t] ~ dbern(phi.totalF[i,j,t])
        
        # Conditional probability of survival for males
        phi.totalM_F[i,j,t] <- equals(pf[i,t],j)*(afmat[i,j,t]*(Phifm[t-1]/PhiF[t-1]) + 
                                         (1 - afmat[i,j,t])*(Phim0[t-1]/(1-PhiF[t-1]))) + 
          (1 - equals(pf[i,t],j))*am[j,t-1]
        
        #P[Y^M_T|Y^F_T]
        ammat[i,j,t] ~ dbern(phi.totalM_F[i,j,t])        
      }
    }
    
    
    # 5b. Survival probability for single Males
    for(i in (nf+1):n){
      for(j in 1:nm){
        # Probability of survival for males
        phi.totalM_F[i,j,t] <- equals(pf[i,t],j)*(PhiM[t-1]) + 
          (1 - equals(pf[i,t],j))
        
        #P[Y^M_T]
        ammat[i,j,t] ~ dbern(phi.totalM_F[i,j,t]*am[j,t-1]) 
      }
    }
    
    
    # 5c. Survival Probability for single Females
    for(i in 1:nf){
      for(j in (nm+1):n){
        # Probability of survival for females
        phi.totalF[i,j,t] <- equals(pf[i,t],j)*(PhiF[t-1]) + 
          (1 - equals(pf[i,t],j))
        
        #P[Y^F_T]
        afmat[i,j,t] ~ dbern(phi.totalF[i,j,t]*af[i,t-1]) 
      }
    }

    # Recover Marginal Survival Distribution for future conditioning
    
    #Females
    for(i in 1:nf){
      # Dummy parameter (its always 1 or zero)
      phif_dummy[i,t] <- max(afmat[i,,t])
      # Missing values are populated based on asurv (hack to get around data/constant problem)
      af[i,t] ~ dbern(phif_dummy[i,t]) 
    }
    
    #Males
    for(j in 1:nm){
      # Dummy parameter (its always 1 or zero)
      phim_dummy[j,t] <- max(ammat[,j,t])
      # Missing values are populated based on asurv (hack to get around data/constant problem)
      am[j,t] ~ dbern(phim_dummy[j,t]) 
    } 
    
    ######################
    # 6. Joint Recapture #
    ######################
    
    
    
  ###########################
  ### Prior distributions ###
  ###########################
  
  # Build Partnership probabilities (unconditioned)
  # !!!! Need to modify histories to be latent for now just assume its some random covariate!!!
  for(t in 1:k){
    # Linear function 
    eta[1:n, 1:n, t] <- beta0 + beta1*histories[1:n, 1:n, t]
    #Softmax by row 
    for(i in 1:n){
      for(j in 1:n){
        psi_raw[i, j ,t] <- exp(eta[i, j, t])
      }
      psi[i,1:n,t] <- psi_raw[i,1:n,t]/sum(psi_raw[i, 1:n ,t])
    }
  }
  
  
  #prior for probs
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