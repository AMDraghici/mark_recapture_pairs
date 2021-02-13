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
      amating_f[i,1] ~ dbern(delta[1]*recruit_f[i,t])
    }
    
    # Male mating at time 1
    for(j in 1:nm){
      amating_m[j,1] ~ dbern(delta[1]*recruit_m[j,t])
    }
  
   # !!!!!!!!!! ADD THE OTHER TIME 1 EVENTS (IF NEEDED)
  
    # Model Joint Partnership/Survival/Recapture outcomes 
    for(t in 2:k){
    
      #######################
      # 3. Decision to Mate #
      #######################
      
      # Female Mating Choice at time t
      for(i in 1:nf){
        amating_f[i,t] ~ dbern((1-af[i,t-1]) + (af[i,t-1]) * recruit_f[i,t] * delta[t])
      }
       
      # Male Mating Choice at time t
      for(j in 1:nm){
         amating_m[j,t] ~ dbern((1-am[j,t-1]) + (am[j,t-1]) * recruit_m[j,t] * delta[t])
      } 
       
      #####################
      # 4. Mate Selection #
      #####################
       
      # Compute conditional psi 
      
      # Partnership status given recruitment and mating
      psi_cond[1:nf,1:nm,t] <- psi[1:nf, 1:nm, t] * (recruit_f[1:nf,t] %*% t(recruit_m[1:nm,t])) * (amating_f[1:nf,t] %*% t(amating_m[1:nm,t]))
      # Fill out dummy quadrants of the matrix (2,3,4)
      psi_cond[(nf+1):n, (nm+1):n, t] <- psi[(nf+1):n, (nm+1):n, t]
      psi_cond[(nf+1):n, 1:nm, t] <- psi[(nf+1):n, 1:nm, t]
      psi_cond[1:nf, (nm+1):n, t] <- psi[1:nf, (nm+1):n, t]
      
      # Assign mate index
      for(i in 1:n){
        apairs[i,,t] ~ dsample(psi_cond[i, 1:n, t], 1)
      } 
       
      #########################
      # Parameter Constraints #
      #########################
      
      # Partner choice cannot occur more than once at t (no polygamy)
      is_chosen_sum[1:n,t] <- apairs[1:n,1:n,t] %*% rep(1,n)
      for(i in 1:n){
        zeroes_vec[i,t] ~ dinterval(is_chosen_sum[i,t], 1) # force each choice to occur between [0, 1] times. 
      }
      
      #####################
      # 5. Joint Survival #
      #####################
      
      for(i in 1:n){
        for(j in 1:n){
          # Building out the conditional survival states 
          #Survival Matrix (previous state indicators determine prob)
          phi.total1[i, j, t] <- equals(apairs[i,j,t], 1)*
            (equals(af[i,t-1], 1)*equals(am[j,t-1], 1)*Phi00[t-1] + 
               equals(af[i,t-1], 1)*equals(am[j,t-1], 0)*(1-PhiF[t-1]) + 
               equals(af[i,t-1], 0)*equals(am[j,t-1], 1)*(1-PhiM[t-1]) + 
               equals(af[i,t-1], 0)*equals(am[j,t-1], 0)) + equals(apairs[i,j,t], 0)
          
          phi.total2[i, j, t] <- equals(apairs[i,j,t], 1)*
            (equals(af[i,t-1], 1)*equals(am[j,t-1], 1)*Phif0[t-1] + 
               equals(af[i,t-1], 1)*equals(am[j,t-1], 0)*PhiF[t-1])
          
          phi.total3[i, j, t] <- equals(apairs[i,j,t], 1)*
            (equals(af[i,t-1], 1)*equals(am[j,t-1], 1)*Phim0[t-1] + 
               equals(af[i,t-1], 0)*equals(am[j,t-1], 1)*PhiM[t-1])
          
          phi.total4[i, j, t] <- equals(apairs[i,j,t], 1)*
            (equals(af[i,t-1], 1)*equals(am[j,t-1], 1)*Phifm[t-1])
          
          #Assign survival likelihoods 
          asurv[i, j, t] ~ dcat(c(phi.total1[i,j,t],phi.total2[i,j,t],phi.total3[i,j,t],phi.total4[i,j,t]))
        }
      }
      
      # Recover Marginal Survival Distribution 
      
      #Females
      for(i in 1:nf){
        # Dummy parameter (its always 1 or zero)
        phif_dummy[i,t] <- equals(max(asurv[i,,t]),4) + equals(max(asurv[i,,t]),2) 
        # Missing values are populated based on asurv (hack to get around data/constant problem)
        af[i,t] ~ dbern(phif_dummy[i,t]) 
      }
      
      #Males
      for(j in 1:nm){
        # Dummy parameter (its always 1 or zero)
        phim_dummy[j,t] <- equals(max(asurv[,j,t]),4) + equals(max(asurv[,j,t]),3) 
        # Missing values are populated based on asurv (hack to get around data/constant problem)
        af[j,t] ~ dbern(phim_dummy[j,t]) 
      } 
      
      ######################
      # 6. Joint Recapture #
      ######################
      
      # Building out the Joint recapture states 
      for(i in 1:n){
        for(j in 1:n){
          
          #Recapture Matrix (previous state indicators determine prob)
          p.total1[i, j, t] <- equals(asurv[i,j,t], 4) * P00[t] + 
            equals(asurv[i,j,t], 3) * (1-PF[t]) + 
            equals(asurv[i,j,t], 2) * (1-PM[t]) + 
            equals(asurv[i,j,t], 1)
          
          p.total2[i, j, t] <- equals(asurv[i,j,t], 4) * Pf0[t] + 
            equals(asurv[i,j,t], 3) * PF[t] 
          
          p.total3[i, j, t] <- equals(asurv[i,j,t], 4) * Pm0[t] + 
            equals(asurv[i,j,t], 2) * PM[t] 
          
          p.total4[i, j, t] <- equals(asurv[i,j,t], 4) * Pfm[t]

          # Assign recapture likelihoods 
          rpair[i, j, t] ~ dcat(c(p.total1[i, j, t],p.total2[i, j, t],p.total3[i, j, t],p.total4[i, j, t]))
        }
      }
    }
  
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
      psi[i, 1:n ,t] <- exp(eta[i,1:n,t]) # unormalized but JAGS handles this anyway
    }
  }
  
  #prior for probs
  beta0 ~ dnorm(0,1)
  beta1 ~ dnorm(0,1)
  
  for(t in 1:(K-1)){
    
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
    gu[t] <- min(sqrt(OR.Phi[t]), 1/sqrt(OR.Phi[t])) 
    gl[t] <- -min(sqrt(OP.Phi[t],1/sqrt(OP.Phi[t]))) 
    
    # Recapture Rates (Rho)
    ru[t] <- min(sqrt(OR.Phi[t]), 1/sqrt(OR.Phi[t])) 
    rl[t] <- -min(sqrt(OP.Phi[t],1/sqrt(OP.Phi[t]))) 
    
    # Odds Ratio and Product of Survival and Recapture Rates
    OP.Phi[t] <- odds.phiF[t]*odds.phiM[t]
    OR.Phi[t] <- odds.PF[t]/odds.PM[t]
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
    eps[t] ~ dunif(0,1)
    
  } 
}