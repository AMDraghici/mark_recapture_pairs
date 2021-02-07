#### Multiple-Partners Joint-Fates Extended Arnason-Schwarz Model 

model{
    ############################
    ### 1. Data Augmentation ###
    ############################
    
    # for(n in 1:N){
    #   
    # }
    
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
    
    for(i in 1:nf){
      amating_f[i,1] ~ dbern(delta[1]*recruit_f[i,t])
    }
    
    
    for(j in 1:nm){
      amating_m[j,1] ~ dbern(delta[1]*recruit_m[j,t])
    }
  
   # !!!!!!!!!! ADD THE OTHER TIME 1 EVENTS (IF NEEDED)
  
    for(t in 2:k){
    
      #######################
      # 3. Decision to Mate #
      #######################
      
      # Female Mating Choice
      for(i in 1:nf){
        amating_f[i,t] ~ dbern((1-af[i,t-1]) + (af[i,t-1]) * recruit_f[i,t] * delta[t])
      }
       
      # Male Mating Choice
      for(j in 1:nm){
       amating_m[j,t] ~ dbern((1-af[j,t-1]) + (af[j,t-1]) * recruit_m[j,t] * delta[t])
      } 
       
      #####################
      # 4. Mate Selection #
      #####################
       
      # Compute conditional psi 
      psi_cond[1:nf,1:nm,t] <- psi[1:nf, 1:nm, t] * (recruit_f[1:nf,t] %*% t(recruit_m[1:nm,t])) * (amating_f[1:nf,t] %*% t(amating_m[1:nm,t]))
      psi_cond[(nf+1):n, (nm+1):n, t] <- psi[(nf+1):n, (nm+1):n, t]
      
      # Assign mate index
      for(e in 1:n){
        apf[e,t] ~ dcat(psi_cond[e, 1:n ,t])
      } 
       
      #####################
      # 5. Joint Survival #
      #####################
      
      # Building out the conditional survival states 
      
      
      # Assign survival likelihoods 
      for(i in 1:n){
        for(j in 1:n){
          apair[i, j, t] ~ dcat(c(phi.total1[i,j,t],phi.total2[i,j,t],phi.total3[i,j,t],phi.total4[i,j,t]))
        }
      }
      
      # Recover Marginal Survival Distribution 
      
      
      
      ######################
      # 6. Joint Recapture #
      ######################
      
      # Building out the conditional survival states 
      
      # Assign recapture likelihoods 
      for(i in 1:n){
        for(j in 1:n){
          rpair[i, j, t] ~ dcat(c(p.total1[i,j,t],p.total2[i,j,t],p.total3[i,j,t],p.total4[i,j,t]))
        }
      }
    }
  
  #########################
  # Parameter Constraints #
  #########################
  
  # Partner choice cannot occur more than once at t (no polygamy)
  for(e in 1:n){
    is_chosen_sum[e] <- sum(apf[1:n]==j) # # of occurances for choice j 
    zeroes_vec[e] ~ dinterval(is_chosen_sum[e], 1) # force each choice to occur between [0, 1] times. 
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
      psi[i, 1:n ,t] <- exp(eta[i,1:n,t])/sum(eta[i,1:n,t])
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
    gamma[t] ~ dunif(gl[t],gu[t])
    rho[t] ~ dunif(rl[t],ru[t])
    
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


# 
# ##Iterate through each group (pairs and singles)
# for(i in 1:M){
#   #For group i cycle through each possible recapture
#   for(t in (first[i]+1):K){ 
#     
#     #########
#     #Divorce# 
#     #########
#     
#     #Temporary Divorce Transition Matrix - If Partner dies then we can only be unmated 
#     delta.total1[i,t] <-  equals(a[(i - 1) * K + t - 1,2],1) + 
#       equals(a[(i - 1) * K + t - 1,2],2) + 
#       equals(a[(i - 1) * K + t - 1,2],3) +
#       equals(a[(i - 1) * K + t - 1,2],4)*(1-delta[t-1])
#     
#     delta.total2[i,t] <- equals(a[(i - 1) * K + t - 1,2],4)*(delta[t-1])
#     
#     #Update mate state at time t
#     d[(i-1)*K+t,2] ~ dcat(c(delta.total1[i,t],delta.total2[i,t]))
#     
#     ##########
#     #Survival# 
#     ##########
#     
#     #Survival Matrix (previous state indicators determine prob)
#     phi.total1[i,t] <- equals(a[(i - 1) * K + t - 1,2],4)*equals(d[(i-1)*K+t,2],2)*Phi00[t-1] + 
#       equals(a[(i - 1) * K + t - 1,2],4)*equals(d[(i-1)*K+t,2],1)*(1-PhiF[t-1])*(1-PhiM[t-1]) + 
#       equals(a[(i - 1) * K + t - 1,2],2)*(1-PhiF[t-1]) + 
#       equals(a[(i - 1) * K + t - 1,2],3)*(1-PhiM[t-1]) + 
#       equals(a[(i - 1) * K + t - 1,2],1)
#     
#     phi.total2[i,t] <- equals(a[(i - 1) * K + t - 1,2],4)*equals(d[(i-1)*K+t,2],2)*Phif0[t-1] + 
#       equals(a[(i - 1) * K + t - 1,2],4)*equals(d[(i-1)*K+t,2],1)*PhiF[t-1]*(1-PhiM[t-1]) + 
#       equals(a[(i - 1) * K + t - 1,2],2)*PhiF[t-1]
#     
#     phi.total3[i,t] <- equals(a[(i - 1) * K + t - 1,2],4)*equals(d[(i-1)*K+t,2],2)*Phim0[t-1] + 
#       equals(a[(i - 1) * K + t - 1,2],4)*equals(d[(i-1)*K+t,2],1)*(1-PhiF[t-1])*PhiM[t-1] +
#       equals(a[(i - 1) * K + t - 1,2],3)*PhiM[t-1]
#     
#     phi.total4[i,t] <- equals(a[(i - 1) * K + t - 1,2],4)*equals(d[(i-1)*K+t,2],2)*Phifm[t-1] + 
#       equals(a[(i - 1) * K + t - 1,2],4)*equals(d[(i-1)*K+t,2],1)*PhiF[t-1]*PhiM[t-1]
#     
#     ###Update the current group state at time t using the categorical distribution
#     a[(i-1)*K+t,2] ~ dcat(c(phi.total1[i,t],phi.total2[i,t],phi.total3[i,t],phi.total4[i,t]))
#     
#     ###########
#     #Recapture# 
#     ###########
#     
#     #Recapture Matrix (current state indicators determine prob)
#     p.total1[i,t] <- equals(a[(i - 1) * K + t,2],4)*equals(d[(i-1)*K+t,2],2)*P00[t-1] + 
#       equals(a[(i - 1) * K + t,2],4)*equals(d[(i-1)*K+t,2],1)*(1-PF[t-1])*(1-PM[t-1]) + 
#       equals(a[(i - 1) * K + t,2],2)*(1-PF[t-1]) + 
#       equals(a[(i - 1) * K + t,2],3)*(1-PM[t-1]) + 
#       equals(a[(i - 1) * K + t,2],1)
#     
#     p.total2[i,t] <-equals(a[(i - 1) * K + t,2],4)*equals(d[(i-1)*K+t,2],2)*Pf0[t-1] + 
#       equals(a[(i - 1) * K + t,2],4)*equals(d[(i-1)*K+t,2],1)*PF[t-1]*(1-PM[t-1]) +  
#       equals(a[(i - 1) * K + t,2],2)*PF[t-1]
#     
#     p.total3[i,t] <-equals(a[(i - 1) * K + t,2],4)*equals(d[(i-1)*K+t,2],2)*Pm0[t-1] + 
#       equals(a[(i - 1) * K + t,2],4)*equals(d[(i-1)*K+t,2],1)*(1-PF[t-1])*PM[t-1] + 
#       equals(a[(i - 1) * K + t,2],3)*PM[t-1]
#     
#     p.total4[i,t] <-equals(a[(i - 1) * K + t,2],4)*equals(d[(i-1)*K+t,2],2)*Pfm[t-1] + 
#       equals(a[(i - 1) * K + t,2],4)*equals(d[(i-1)*K+t,2],1)*PF[t-1]*PM[t-1] 
#     
#     ##Update Recapture information with categorical distn
#     X[(i-1)*K+t,2] ~ dcat(c(p.total1[i,t],p.total2[i,t],p.total3[i,t],p.total4[i,t])) 
#     
#     
#   }
# } 
