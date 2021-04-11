#### Standard CJS Model 

model{
  
  # CJS Likelihood -----------------------------------------------------------------------------------
  for(i in 1:n){
    for(t in (initial_entry[i]+1):k){
      
      phi[i,t] <- equals(female[i], 1) * phiF + equals(female[i], 0)*phiM
      p[i, t] <- equals(female[i], 1)*pF + equals(female[i], 0)*pM
      
      a[i,t] ~ dbern(phi[i,t] * a[i,t-1])
      x[i,t] ~ dbern(p[i, t] * a[i,t])
    }
  }
  
  # Priors--------------------------------------------------------------------------------------------
  
  # Survival by sex
  phiF ~ dbeta(1,1)
  phiM ~ dbeta(1,1)
  
  # Recapture by sex
  pF ~ dbeta(1,1)
  pM ~ dbeta(1,1)
}