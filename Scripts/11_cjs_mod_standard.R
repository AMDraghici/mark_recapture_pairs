#### Standard CJS Model 

model{
  
  # CJS Likelihood -----------------------------------------------------------------------------------
  for(i in 1:n){
    for(t in (initial_entry[i]+1):k){
      x[i,t] ~ dbern(p[i] * a[i,t])
      a[i,t] ~ dbern(phi[i] * a[i,t-1])
    }
    
  }
  
  # Priors--------------------------------------------------------------------------------------------
  
  for(i in 1:n){
    p[i] <- equals(female[i], 1)* pF + equals(female[i], 0)* pM
    phi[i] <-  equals(female[i], 1) * phiF + equals(female[i], 0)* phiM
  }
  
  
  # Survival by sex
  phiF ~ dbeta(1,1)
  phiM ~ dbeta(1,1)
  
  # Recapture by sex
  pF ~ dbeta(1,1)
  pM ~ dbeta(1,1)
}