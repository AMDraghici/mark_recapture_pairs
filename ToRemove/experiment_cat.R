model{
  
  # for(t in 1:k){
  #   for(i in 2:n){
  #       pairs2[i,t] ~ dcat(psi[i,1:n,t])
  #     }
  # }
  
  psi2[1:2] <- psi[1:2,1,1]
  psi2[1] <- psi[3,1,1]
  
}
