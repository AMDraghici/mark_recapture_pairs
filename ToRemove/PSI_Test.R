sort(lapply(2:29, function(x) which(!is.na(jags_data$apairs_f[,x]) & (jags_data$apairs_f[,x] < nm + 1)))[[t+1]]) #females
sort(lapply(2:29, function(x) jags_data$apairs_f[which(!is.na(jags_data$apairs_f[,x]) & (jags_data$apairs_f[,x] < nm + 1)),x])[[t+1]]) #males

list1 <- list()
list2 <- list()
for(i in 1:k){
  list1[[i]] <- which(colSums(jags_data$psi[1:nf,1:nm,i])==1)
  list2[[i]] <- which(rowSums(jags_data$psi[1:nf,1:nm,i])==1)
}

sort(list1[[t]]) # males

sort(list2[[t]]) #females

(jags_data$apairs_f[,t+1][!is.na(jags_data$apairs_f[,t+1])])

for(i in 1:k){
  if(!all.equal(list2[[i]],which(!is.na(jags_data$apairs_f[,i+1])))) stop(i)
}
