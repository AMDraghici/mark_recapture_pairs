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
  for(i in 1:nf){
    for(j in 1:nm){
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
      logit(prob_repartner[i,t]) <- 0 + beta1*inprod(histories[i, 1:nm , t],  t(apairs[i + 1, 1:nm + 1, t])) 
      arepartner[i,t] ~ dbern(prob_repartner[i,t] * amating_f[i,t] * inprod(apairs[i + 1, 1:nm + 1, t], amating_m[1:nm,t]))
    }
    
    # 4. Mate Selection -------------------------------------------------------------------------------------------------------------------
    # JAGs does not play well with multinomial trials
    # Solution - series of conditional binomial trials with constraints for singles
    # Probabilities have to be unconditioned after sampling
    
    # Build Homogeneous Partnership probabilities 
    for(i in 1:nf){
      # Flat likelihood of mating conditional on decision to mate
      psi_female[i, 1:nm, t] <- amating_f[i,t] * amating_m[1:nm,t] * (1 - arepartner[i,t]) +#* (1 - apairs[i+1,1:nm+1,t]) + 
                             arepartner[i,t] * apairs[i+1,1:nm+1,t]
    }
      
    # Exclude impossible partnerships based on repartner status 
    for(j in 1:nm){
      psi_raw_male[j,1:nf,t] <- apairs[1:nf + 1,j+1, t] * arepartner[1:nf,t]
      psi_male[j,1:nf,t] <- (1-sum(psi_raw_male[j,1:nf,t]))*(1 - psi_raw_male[j,1:nf,t]) + psi_raw_male[j,1:nf,t] 
    }
      
    psi[1:nf, 1:nm, t] <- psi_female[1:nf, 1:nm, t] * t(psi_male[1:nm, 1:nf, t])
    
    # Attempts at partnerships forming
    # Monogamous pairings only 
    for(i in 1:nf){

      # Normalize with respect to males who found pairs
      for(j in 1:nm){
        psi_cond_raw[i,j,t] <- psi[i,j,t]*(1 - sum(apairs[1:i, j+1, t+1]))
      }
      psi_cond[i,1:nm,t] <-  psi_cond_raw[i,1:nm,t]/(sum(psi_cond_raw[i, 1:nm ,t]) + equals(sum(psi_cond_raw[i, 1:nm ,t]), 0))
      
      # Assign pairs using normalized probability density
      for(j in 1:nm){
        apairs[i+1,j+1,t+1] ~ dbern((psi_cond[i,j,t]/(sum(psi_cond[i,j:nm,t]) + equals(sum(psi_cond[i,j:nm,t]),0))) * (1 - sum(apairs[i+1, 1:j, t+1])))
      }
    }
    
    # If no pairs for male designate as single
    for(j in 1:nm){
      single_male[j, t] <- 1 - sum(apairs[1:nf+1,j+1,t+1]) # if no pairs then classify as single
    }
    
    
    # Update Total History for Next Time Step
    for(i in 1:nf){
      for(j in 1:nm){
        histories[i, j, t+1] <- sum(apairs[i+1, j+1, 1:(t+1)])
      }
    }
    
    # 5. Joint Survival ---------------------------------------------------------------------------------------------------------------------
    
    
    # Marginal Survival Event for Females in the Population (P[Y^F_T])
    for(i in 1:nf){
      af[i,t+1] ~ dbern(PhiF * af[i,t])
    }
    
    # Marginal Recapture Event for Males in the Population (P[X^M_T|X^F_T])
    for(j in 1:nm){
      
      # Probability of male surviving given partnership and partner recapture status
      phi.totalM[j, t] <- single_male[j,t] * PhiM + # Male was single
        (1 - single_male[j,t]) * (inprod(apairs[1:nf+1,j+1,t+1], af[1:nf,t+1]) * (Phifm/PhiF) + # Male mated and female surived
                                    (1 - inprod(apairs[1:nf+1,j+1,t+1], af[1:nf,t+1])) * (Phim0/(1-PhiF))) # Male mated and female perished
      
      # Draw Survival Event 
      am[j, t+1] ~ dbern(phi.totalM[j,t] * am[j,t])    
    }
    
    # 6. Joint Recapture --------------------------------------------------------------------------------------------------------------------
    
    
    # Marginal Recapture Event for Females in the Population (P[X^F_T])
    for(i in 1:nf){
      recap_f[i,t] ~ dbern(PF * af[i,t+1])
    }
    
    
    # Marginal Recapture Event for Males in the Population (P[X^M_T|X^F_T])
    for(j in 1:nm){
      
      # Probability of male being captured given partnership and partner recapture status
      p.totalM[j, t] <- single_male[j,t] * PM + # Male was single
        (1 - single_male[j,t]) * (inprod(apairs[1:nf+1,j+1,t+1], recap_f[1:nf,t]) * (Pfm/PF) + # Male mated and female captured
                                  (1 - inprod(apairs[1:nf+1,j+1,t+1], recap_f[1:nf,t])) * (Pm0/(1-PF))) # Male mated and female not captured
      
      # Draw Recapture Probability 
      recap_m[j, t] ~ dbern(p.totalM[j,t]*am[j,t+1])    
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