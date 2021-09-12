#### Multiple-Partners Joint-Fates Extended CJS Model 

# Data Dictionary --------------------------------------------------------------------------------------------------------------------------
# 
#  1. INDEXING
#  
#   - nf (integer): Number of females in the sample population
#   - nm (integer): Number of males in the sample population 
#   - k  (integer): Number of equidistant sampling occasions throughout the study 
#   
#  2. OBSERVED DATA
#  
#   - recap_f (binary - matrix: nrow = nf+1 , ncol = k): 
#         Recapture histories for females in the sample (1 = seen, 0 = unobserved)
#         Final row is a dummy entry and are used for indexing purposes in JAGS
#   - recap_m (binary - matrix: nrow = nm+1 , ncol = k):
#         Recapture histories for males in the sample (1 = seen, 0 = unobserved)
#         Final row is a dummy entry and are used for indexing purposes in JAGS
#
#
#  3. PARTIALLY OBSERVED DATA
#  
#  
#   - amating_f (binary - matrix nrow = nf +1 , ncol = k): 
#         Decision to mate for females in the sample from 1:k (1 = mating at k, 0 = not mating at k)
#         Final row is a dummy entry and are used for indexing purposes in JAGS
#   - amating_m (binary - matrix nrow = nm +1 , ncol = k): 
#         Decision to mate for males in the sample from 1:k (1 = mating at k, 0 = not mating at k)
#         Final row is a dummy entry and are used for indexing purposes in JAGS
#   - apairs_f (integer - matrix, nrow = nf, ncol = K + 1): 
#         Partnership records. Rows = Females, Cols = Time, Entries = ID of Male at time k. 
#         If female choses not to mate at time k then a dummy entry of (nm + 1) is assigned
#         First column is a dummy entry of nm + 1 for all females and is used for indexing purposes in JAGS
#
#   - af (binary - matrix nrow = nf+1, ncol = K + 1): 
#           Survival history females (1 for alive, 0 for dead)
#                First column and final row are dummy entries and are used for indexing purposes in JAGS
#   - am (binary - matrix nrow = nm+1, ncol = K + 1): 
#           Survival history males (1 for alive, 0 for dead)
#                First column and final row are dummy entries and are used for indexing purposes in JAGS
#
#   - arepartner (binary - matrix nrow = nf, ncol = k)
#           Repeat Partnership history. Rows = Females, Cols = Time
#           1 = Female at time k formed a pair with partner from time k - 1 
#           0 = otherwise 
#           Allows us to estimate fidelity for partners
#           Conditional on decision to mate and survival to avoid biasing the estimates
#
#   - recruit_f (binary - matrix, nrow = nf+1, ncol = k)
#           Recruitment into population status for females. 1 if currently in population 0 o/w
#           Value at k must be equal to 1 for all animals 
#         Final row is a dummy entry and are used for indexing purposes in JAGS
#   - recruit_m (binary - matrix, nrow = nm+1, ncol = k)
#           Recruitment into population status for males. 1 if currently in population 0 o/w
#           Value at k must be equal to 1 for all animals
#         Final row is a dummy entry and are used for indexing purposes in JAGS
#
#  4. MISCELLANEOUS DATA CONSTRUCTS from Pre-processing
#
#   - PSI (binary, array, nrow = nf, ncol = nm + 1, slices = k)
#       At each time point there is a matrix indicating possible pairings
#       A possible pairing is denoted by a 1 and an impossible pairing by a zero 
#       psi[10,10,2] = 1 indicates that male 10 and female 10 COULD form a pair at time 2
#       This is used to impose restrictions on the random partner sampling step to ensure..
#       ...that impossible pairs do not form as pairs are defined as one male and one female.
#       Final column is a dummy column for indexing purposes
#       
#--------------------------------------------------------------------------------------------------------------------------------------------

model{
  
  # 1. Recruitment into population-----------------------------------------------------------------------------------------------------------
  
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
    
    # 2. Decision to Mate -------------------------------------------------------------------------------------------------------------------
    
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
      prob_repartner[i,t] <- ilogit(beta0 + beta1*histories[i, apairs_f[i,t] , t]) * psi[i, apairs_f[i,t], t]
      arepartner[i,t] ~ dbern(prob_repartner[i,t] * amating_f[i,t] * amating_m[apairs_f[i,t],t])
    }
    
    # 3. Mate Selection -------------------------------------------------------------------------------------------------------------------
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
      
      # Is male j available at time t 
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
    apairs_f[1,t+1] ~ dcat(psi_cond2[1, 1:(nm+1), t] + 1e-36)
    single_female[1,t] <- psi_cond2[1, (nm+1), t]
    
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
      # Adding + 1e-15 is for numerical stability 
      apairs_f[i,t+1] ~ dcat(psi_cond2[i, 1:(nm+1), t] + 1e-36)
      
      # Designate as single if no males available or if choosing not to mate
      # Psi_cond2 == 1 means that must be single
      single_female[i,t] <- psi_cond2[i,(nm+1),t]
    }
    
    # Update Total History for Next Time Step
    for(i in  1:nf){
      for(j in  1:(nm+1)){
        histories[i, j, t+1] <- histories[i, j, t] + equals(apairs_f[i,t+1],j)*(1-single_female[i,t])
      }
    }
    
    # 4. Joint Survival ---------------------------------------------------------------------------------------------------------------------
    
    
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
    
    # 5. Joint Recapture --------------------------------------------------------------------------------------------------------------------
    
    
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
  
  
  # 6. Prior Distributions-------------------------------------------------------------------------------------------------------------------
  
  # Recruitment 
  for(t in 1:k){
    eps[t] ~ dbeta(1,1)
  }
  
  # Attempt to Mate 
  delta ~ dbeta(3,2)
  
  # Pairs reforming
  beta0 ~ dnorm(0, 2)
  beta1 ~ dnorm(0, 1)
  
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
}