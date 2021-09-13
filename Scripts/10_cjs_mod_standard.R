#### Standard CJS Model 

model{
  
  # CJS Likelihood -----------------------------------------------------------------------------------
  for(i in 1:n){
    for(t in (initial_entry[i]+1):k){
      x[i,t] ~ dbern(p[i] * a[i,t+1])
      a[i,t+1] ~ dbern(phi[i] * a[i,t])
    }
    
  }
  
  # Priors--------------------------------------------------------------------------------------------
  
  for(i in 1:n){
    p[i] <- female[i] * pF + (1-female[i]) * pM
    phi[i] <- female[i] * phiF + (1-female[i])* phiM
  }
  
  
  # Survival by sex
  phiF ~ dbeta(1,1)
  phiM ~ dbeta(1,1)
  
  # Recapture by sex
  pF ~ dbeta(1,1)
  pM ~ dbeta(1,1)
}