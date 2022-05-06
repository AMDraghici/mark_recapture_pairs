#### Standard CJS Model 

model{
  
  # Data Augmentation --------------------------------------------------------------------------------
  for(i in 1:n){
    z[i] ~ dbern(xi)
  }

  ## Compute population size
  N <- sum(z[1:n])
  
  # Recruit Likelihood -------------------------------------------------------------------------------
  for(i in 1:n){
    recruit[i, 1] ~  dbern(eps[1])
    for(t in 2:(k-1)){
      recruit[i,t] ~ dbern(recruit[i, t-1] + (1-recruit[i, t-1]) * eps[t])
    } 
  }
  
  # JS Likelihood -----------------------------------------------------------------------------------
  for(i in 1:n){
    ## 1) Survival
    for(t in 2:k){
      a[i,t] ~ dbern(phi[i] * a[i,t-1] * recruit[i,t]+ (1 - recruit[i,t]))
    }
    
    for(t in 1:k){
      x[i,t] ~ dbern(p[i] * a[i,t] * recruit[i,t] * z[i])
    }
  }
  
  # Priors--------------------------------------------------------------------------------------------
  
  # Split probs by sex 
  for(i in 1:n){
    p[i] <- female[i] * pF + (1-female[i]) * pM
    phi[i] <- female[i] * phiF + (1-female[i])* phiM
  }
  
  # # Data augmentation
  xi ~ dbeta(1.0,1.0)

  # Survival by sex
  phiF ~ dbeta(1,1)
  phiM ~ dbeta(1,1)
  
  # Recapture by sex
  pF ~ dbeta(1,1)
  pM ~ dbeta(1,1)
  
  # Recruitment 
  for(t in 1:k){
    eps[t] ~ dbeta(1,1)
  }
  
}