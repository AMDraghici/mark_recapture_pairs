####Delta-Gamma-Rho-CJS Model - extended Rho-CJS from challenger (2010)###

model{
  ##Iterate through each group (pairs and singles)
  for(i in 1:M){
    #For group i cycle through each possible recapture
    for(j in (first[i]+1):K){ 
      
      #########
      #Divorce# 
      #########
      
      #Temporary Divorce Transition Matrix - If Partner dies then we can only be unmated 
      delta.total1[i,j] <-  equals(a[(i - 1) * K + j - 1,2],1) + 
        equals(a[(i - 1) * K + j - 1,2],2) + 
        equals(a[(i - 1) * K + j - 1,2],3) +
        equals(a[(i - 1) * K + j - 1,2],4)*(1-delta[j-1])
      
      delta.total2[i,j] <- equals(a[(i - 1) * K + j - 1,2],4)*(delta[j-1])
      
      #Update mate state at time j
      d[(i-1)*K+j,2] ~ dcat(c(delta.total1[i,j],delta.total2[i,j]))
      
      ##########
      #Survival# 
      ##########
      
      #Survival Matrix (previous state indicators determine prob)
      phi.total1[i,j] <- equals(a[(i - 1) * K + j - 1,2],4)*equals(d[(i-1)*K+j,2],2)*Phi00[j-1] + 
        equals(a[(i - 1) * K + j - 1,2],4)*equals(d[(i-1)*K+j,2],1)*(1-PhiF[j-1])*(1-PhiM[j-1]) + 
        equals(a[(i - 1) * K + j - 1,2],2)*(1-PhiF[j-1]) + 
        equals(a[(i - 1) * K + j - 1,2],3)*(1-PhiM[j-1]) + 
        equals(a[(i - 1) * K + j - 1,2],1)
      
      phi.total2[i,j] <- equals(a[(i - 1) * K + j - 1,2],4)*equals(d[(i-1)*K+j,2],2)*Phif0[j-1] + 
        equals(a[(i - 1) * K + j - 1,2],4)*equals(d[(i-1)*K+j,2],1)*PhiF[j-1]*(1-PhiM[j-1]) + 
        equals(a[(i - 1) * K + j - 1,2],2)*PhiF[j-1]
      
      phi.total3[i,j] <- equals(a[(i - 1) * K + j - 1,2],4)*equals(d[(i-1)*K+j,2],2)*Phim0[j-1] + 
        equals(a[(i - 1) * K + j - 1,2],4)*equals(d[(i-1)*K+j,2],1)*(1-PhiF[j-1])*PhiM[j-1] +
        equals(a[(i - 1) * K + j - 1,2],3)*PhiM[j-1]
      
      phi.total4[i,j] <- equals(a[(i - 1) * K + j - 1,2],4)*equals(d[(i-1)*K+j,2],2)*Phifm[j-1] + 
        equals(a[(i - 1) * K + j - 1,2],4)*equals(d[(i-1)*K+j,2],1)*PhiF[j-1]*PhiM[j-1]
      
      ###Update the current group state at time j using the categorical distribution
      a[(i-1)*K+j,2] ~ dcat(c(phi.total1[i,j],phi.total2[i,j],phi.total3[i,j],phi.total4[i,j]))
      
      ###########
      #Recapture# 
      ###########
      
      #Recapture Matrix (current state indicators determine prob)
      p.total1[i,j] <- equals(a[(i - 1) * K + j,2],4)*equals(d[(i-1)*K+j,2],2)*P00[j-1] + 
        equals(a[(i - 1) * K + j,2],4)*equals(d[(i-1)*K+j,2],1)*(1-PF[j-1])*(1-PM[j-1]) + 
        equals(a[(i - 1) * K + j,2],2)*(1-PF[j-1]) + 
        equals(a[(i - 1) * K + j,2],3)*(1-PM[j-1]) + 
        equals(a[(i - 1) * K + j,2],1)
      
      p.total2[i,j] <-equals(a[(i - 1) * K + j,2],4)*equals(d[(i-1)*K+j,2],2)*Pf0[j-1] + 
        equals(a[(i - 1) * K + j,2],4)*equals(d[(i-1)*K+j,2],1)*PF[j-1]*(1-PM[j-1]) +  
        equals(a[(i - 1) * K + j,2],2)*PF[j-1]
      
      p.total3[i,j] <-equals(a[(i - 1) * K + j,2],4)*equals(d[(i-1)*K+j,2],2)*Pm0[j-1] + 
        equals(a[(i - 1) * K + j,2],4)*equals(d[(i-1)*K+j,2],1)*(1-PF[j-1])*PM[j-1] + 
        equals(a[(i - 1) * K + j,2],3)*PM[j-1]
      
      p.total4[i,j] <-equals(a[(i - 1) * K + j,2],4)*equals(d[(i-1)*K+j,2],2)*Pfm[j-1] + 
        equals(a[(i - 1) * K + j,2],4)*equals(d[(i-1)*K+j,2],1)*PF[j-1]*PM[j-1] 
      
      ##Update Recapture information with categorical distn
      X[(i-1)*K+j,2] ~ dcat(c(p.total1[i,j],p.total2[i,j],p.total3[i,j],p.total4[i,j])) 
      
      
    }
  } 
  
  ### Prior distributions ###
  for(j in 1:(K-1)){
    
    #Joint Survival probabilities for paired individuals
    Phi00[j] <- 1 - PhiF[j] - PhiM[j] + Phifm[j]
    Phif0[j] <- PhiF[j] - Phifm[j]
    Phim0[j] <- PhiM[j] - Phifm[j]
    Phifm[j] <- gamma[j]*sig.PhiF[j]*sig.PhiM[j] + PhiF[j]*PhiM[j]
    
    #Joint Capture probabilities for paired individuals
    P00[j] <- 1 - PF[j] - PM[j] + Pfm[j]
    Pf0[j] <- PF[j] - Pfm[j]
    Pm0[j] <- PM[j] - Pfm[j]
    Pfm[j] <- rho[j]*sig.PF[j]*sig.PM[j] + PF[j]*PM[j]
    
    
    ###Upper and lower bounds for covariates (g - gamma(surv), r - rho(recap))
    gu[j] <- min(1,
                 (PhiF[j] - PhiF[j]*PhiM[j])/(sig.PhiF[j]*sig.PhiM[j]),
                 (PhiM[j] - PhiF[j]*PhiM[j])/(sig.PhiF[j]*sig.PhiM[j])
    ) 
    
    gl[j] <- max(-1,
                 (-PhiF[j]*PhiM[j])/(sig.PhiF[j]*sig.PhiM[j]),
                 (PhiF[j]+PhiM[j] - PhiF[j]*PhiM[j] - 1)/(sig.PhiF[j]*sig.PhiM[j])
    ) 
    
    ru[j] <- min(1,
                 (PF[j] - PF[j]*PM[j])/(sig.PF[j]*sig.PM[j]),
                 (PM[j] - PF[j]*PM[j])/(sig.PF[j]*sig.PM[j])
    ) 
    
    rl[j] <- max(-1,
                 (-PF[j]*PM[j])/(sig.PF[j]*sig.PM[j]),
                 (PF[j]+PM[j] - PF[j]*PM[j] - 1)/(sig.PF[j]*sig.PM[j])
    ) 
    
    ###Binomial standard deviation for survival and recapture 
    sig.PhiF[j] <- sqrt(PhiF[j]*(1-PhiF[j]))
    sig.PhiM[j] <- sqrt(PhiM[j]*(1-PhiM[j]))
    
    sig.PF[j] <- sqrt(PF[j]*(1-PF[j]))
    sig.PM[j] <- sqrt(PM[j]*(1-PM[j]))
    
    ##Uninformative priors for survival/recapture/divorce state
    delta[j] ~ dbeta(1,1)
    
    PhiF[j] ~ dbeta(1,1)
    PhiM[j] ~ dbeta(1,1)
    
    PF[j] ~ dbeta(1,1)
    PM[j] ~ dbeta(1,1)
    
    ##Pulling Coverates from unif distns on thier upper and lower bounds for j 
    gamma[j] ~ dunif(gl[j],gu[j])
    rho[j] ~ dunif(rl[j],ru[j])
    
  } 
}