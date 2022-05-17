# UNIT TESTS

# Repartner check
for(t in 1:jags_data$k){
  index <- which(!is.na(jags_data$arepartner[,t]))
  print(all(jags_data$arepartner[index,t] == jags_data$repartner[index,t]))
} 

# apairs check
for(t in 1:jags_data$k){
  index <- which(!is.na(jags_data$apairs_f[,t]))
  print(all(jags_data$pairs_f[index,t] == jags_data$apairs_f[index,t]))
} 

for(t in 1:jags_data$k){
  index <- which(!is.na(jags_data$apairs_m[,t]))
  print(all(jags_data$pairs_m[index,t] == jags_data$apairs_m[index,t]))
} 

# mating check
for(t in 1:jags_data$k){
  index <- which(!is.na(jags_data$amating_f[1:jags_data$nf,t]))
  print(all(jags_data$amating_f[index,t] == jags_data$mating_f[index,t]))
} 

for(t in 1:jags_data$k){
  index <- which(!is.na(jags_data$amating_m[1:jags_data$nm,t]))
  print(all(jags_data$amating_m[index,t] == jags_data$mating_m[index,t]))
} 

t <- 3
index <- which(!is.na(jags_data$arepartner[,t]))
jags_data$arepartner[index,t]
jags_data$repartner[index,t]

for(t in 1:jags_data$k){
  rpair <- which((jags_data$rpair[1:nf,t]))
  print(all(jags_data$pairs_f[index,t-1] == jags_data$apairs_f[index,t]))
  print(all(jags_data$pairs_f[index,t-1] == jags_data$apairs_f[index,t]))
} 

i <- 100
data.frame(pairs_f = jags_data$pairs_f[i,],
           repartner = jags_data$repartner[i,],
           apairs_f = jags_data$apairs_f[i,2:31],
           arepartner = jags_data$arepartner[i,],
           recap_f = jags_data$recap_f[i,], 
           af = jags_data$af[i,2:31],
           sf = jags_data$sf[i,],
           amating_f = jags_data$amating_f[i,],
           mating_f = jags_data$mating_f[i,],
           recruit_f = jags_data$recruit_f[i,],
           time = 1:k)

\nf <-jags_data$nf
nm <-jags_data$nm
k <-jags_data$k


for(i in 1:nf){
  for(t in 1:k){
    if(is.na(jags_data$apf[i,t])) next
    if(jags_data$apf[i,t] == nm+1) next
    
    j <- jags_data$apf[i,t]
    
    if(j != which(jags_data$psi[i,,t]==1)) print(paste0("Bug at female", i, "time ", t))
    
    if(j != nim_inits$apairs_f[i,t]) print(paste0("Init Bug at female", i, "time ", t))
    
  }
}


for(i in 1:nf){
  for(t in 1:(k-1)){
    if(is.na(jags_data$arepartner[i,t])) next
    
    if(jags_data$arepartner[i,t] != nim_inits$arepartner[i,t]) print(paste0("Bug at female", i, "time ", t))
    
    
  }
}



