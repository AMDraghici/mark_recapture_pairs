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
  
  # Initialize History Array (All Zero at time 1)
  for(i in 1:(nf+1)){
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
      prob_repartner[i,t] <- ilogit(beta0 + beta1*histories[i, apairs_f[i,t] , t])
      arepartner[i,t] ~ dbern(prob_repartner[i,t] * amating_f[i,t] * amating_m[apairs_f[i,t],t])
    }
    
    # 4. Mate Selection -------------------------------------------------------------------------------------------------------------------
    # JAGs does not play well with multinomial trials
    # Solution - series of conditional binomial trials with constraints for singles
    # Probabilities have to be unconditioned after sampling
    
    # Build Homogeneous Partnership probabilities 
    for(i in 1:nf){
      # Flat likelihood of mating conditional on decision to mate
      for(j in 1:nm){
        psi_raw[i, j, t] <- psi[i,j,t] * amating_f[i,t] * amating_m[j,t] * (1 - arepartner[i,t]) +
          arepartner[i,t]*equals(apairs_f[i,t],j)
      }
    }
    
    
    #  !!!!!!!!!!!!!!!!!!!!!!! NEED TO ADD FILTER FOR MALE SELECTION HERE  !!!!!!!!!!!!!!!!!!!!!!!
    
    # Attempts at partnerships forming
    # Monogamous pairings only 
    for(i in 1:nf){
     #  !!!!!!!!!!!!!!!!!!!!!!! NEED TO ADD FILTER FOR ALREADY CHOSEN PARTNERS HERE!!!!!!!!!!!!!!!!
      apairs_f[i,t+1] ~ dcat(c(psi_raw[i,1:nm,t],equals(sum(psi_raw[i,1:nm,t]),0)))
      single_female[i,t] = equals(apairs_f[i,t+1],nm+1)
    }
    
    # Update Total History for Next Time Step
    for(i in 1:nf){
      for(j in 1:nm){
        histories[i, j, t+1] <- equals(apairs_f[i,t+1],j)
      }
    }
    
    # 5. Joint Survival ---------------------------------------------------------------------------------------------------------------------
    
    
    # Marginal Survival Event for Males in the Population (P[Y^M_T])
    for(j in 1:nm){
      am[j,t+1] ~ dbern(PhiM * am[j,t])
    }
    
    # Marginal Recapture Event for Males in the Population (P[X^M_T|X^F_T])
    for(i in 1:nf){
      
      # Probability of female surviving given partnership and partner recapture status
      phi.totalF[i, t] <- single_female[i,t] * PhiF + # female was single
        (1 - single_female[i,t]) * (am[apairs_f[i,t+1],t+1] * (Phifm/PhiM) + # Male mated and female surived
                                    (1 - am[apairs_f[i,t+1],t+1]) * (Phif0/(1-PhiM))) # Male mated and female perished
      
      # Draw Survival Event 
      af[i, t+1] ~ dbern(phi.totalF[i,t] * af[i,t])    
    }
    
    # 6. Joint Recapture --------------------------------------------------------------------------------------------------------------------
    
    
    # Marginal Recapture Event for Males in the Population (P[X^M_T])
    for(j in 1:nm){
      recap_m[j,t] ~ dbern(PM * am[j,t+1])
    }
    
    
    # Marginal Recapture Event for females in the Population (P[X^F_T|X^M_T])
    for(i in 1:nf){
      
      # Probability of female being captured given partnership and partner recapture status
      p.totalF[i, t] <- single_female[i,t] * PM + # Female was single
        (1 - single_female[i,t]) * (recap_m[apairs_f[i,t+1],t] * (Pfm/PM) + # Male mated and female captured
                                    (1 - recap_m[apairs_f[i,t+1],t]) * (Pf0/(1-PM))) # Male mated and female not captured
      
      # Draw Recapture Probability 
      recap_f[i, t] ~ dbern(p.totalF[i,t]*af[i,t+1])    
    }
    
  }
  
  
  # 7. Prior Distributions-------------------------------------------------------------------------------------------------------------------
  
  # Recruitment 
  for(t in 1:k){
    eps[t] ~ dbeta(1,1)
  }
  
  # Attempt to Mate 
  delta ~ dbeta(1,1)
  
  # Pairs reforming
  ##beta0 ~ dnorm(0, 1/4)
  beta1 ~ dnorm(0, 1/4)
  beta0 ~ dnorm(0, 1/4)
  
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
  PhiF ~ dbeta(2,2)
  PhiM ~ dbeta(2,2)
  
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
  PF ~ dbeta(2,2)
  PM ~ dbeta(2,2)
}