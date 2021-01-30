#### Multiple-Partners Joint-Fates Extended Arnason-Schwarz Model 

model{
  ############################
  ### 1. Data Augmentation ###
  ############################
  
  for(n in 1:N){
    
  }
  
  ##################################
  # 2. Recruitment into population #
  ##################################
  
  for(t in 1:K){
    
  }
  
  ########################################################
  ### Conditional Partnership/Survival/Recapture Steps ###
  ########################################################
  
   for(t in 1:K){
    
    #######################
    # 3. Decision to Mate #
    #######################
     
    for(i in 1:N){
      
    }
     
    for(i in 1:Nf){ #Females (Mate Anchor Index)
      
      ############################
      # 3. Partnerships at time t#
      ############################
      
      
      for(j in 1:Nm){ #Males
        
        ###############################
        # 4. Survival states at time t#
        ############################### 
        
        
        ################################
        # 5. Recapture states at time t#
        ################################ 
        
      }
    }
     
  }
  
  #########################
  # Parameter Constraints #
  #########################
  
  # Male choice cannot occur more than once at t (no polygamy)
  # Exception "male 1" which represents being single at t (no likelihood contribution)
  for(j in 2:Nm){
    is_chosen_sum[j] <- sum(chosen_index[1:Nm]==j) # # of occurances for choice j 
    zeroes_vec[j] ~ dinterval(is_chosen_sum[j], 1) # force each choice to occur between [0, 1] times. 
  }
  
  
  ###########################
  ### Prior distributions ###
  ###########################
  
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
