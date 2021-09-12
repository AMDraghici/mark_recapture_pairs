model{
  
  for(t in 1:k){
    for(i in 2:n){
      for(j in 2:n){
        pairs[i,j,t] ~ dbern(psi[i,j,t]/(sum(psi[i,j:n,t]) + equals(sum(psi[i,j:n,t]),0))*(1 - sum(pairs[i, 1:(j-1), t])))
      }
    }
  }
  
  
}