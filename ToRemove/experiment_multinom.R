model{
  
  for(t in 1:k){
    for(i in 2:n){
      for(j in 2:n){
        pairs[i,j,t] ~ dmulti(psi[i,j,t],1) #*(1 - sum(pairs[i-1, 1:j, t]))
      }
    }
  }
  
  
}