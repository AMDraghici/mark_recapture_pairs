# UNIT TESTS

# Repartner check
cat("TESTING arepartner","\n")
for(t in 1:(model_data$k-1)){
  x <- list()
  index <- which(!is.na(model_data$arepartner[,t]))
  x[[t]] <- all(model_data$arepartner[index,t] == model_data$repartner[index,t+1])
  if(!all(unlist(x))) cat(paste0("X-test Failure at ", which(unlist(x) == F), " ..."), "\n") 
} 

# apairs check
cat("TESTING apairs_f","\n")
for(t in 1:model_data$k){
  x <- list()
  index <- which(!is.na(model_data$apf[,t]))
  pairs_f_check <- model_data$pairs_f[index,t]
  pairs_f_check[pairs_f_check == model_data$true_pop_m + 1] <- model_data$sample_size_m + model_data$lm + 1
  x[[t]] <- all(pairs_f_check == model_data$apf[index,t])
  if(!all(unlist(x))) cat(paste0("X-test Failure at ", which(unlist(x) == F), " ..."), "\n") 
} 

cat("TESTING apairs_m","\n")
for(t in 1:model_data$k){
  x <- list()
  index <- which(!is.na(model_data$apairs_m[,t]))
  pairs_m_check <- model_data$apairs_m[index,t]
  pairs_m_check[pairs_m_check == model_data$true_pop_f + 1] <- model_data$sample_size_f + model_data$lf + 1
  x[[t]] <- all(pairs_m_check == model_data$apairs_m[index,t])
  if(!all(unlist(x))) cat(paste0("X-test Failure at ", which(unlist(x) == F), " ..."), "\n") 
} 

cat("TESTING amating_f","\n")
for(t in 1:model_data$k){
  x <- list()
  index <- which(!is.na(model_data$amating_f[1:model_data$nf,t]))
  x[[t]] <- all(model_data$amating_f[index,t] == model_data$mating_f[index,t])
  if(!all(unlist(x))) cat(paste0("X-test Failure at ", which(unlist(x) == F), " ..."), "\n") 
} 

cat("TESTING amating_m","\n")
for(t in 1:model_data$k){
  x <- list()
  index <- which(!is.na(model_data$amating_m[1:model_data$nm,t]))
  x[[t]] <- all(model_data$amating_m[index,t] == model_data$mating_m[index,t])
  if(!all(unlist(x))) cat(paste0("X-test Failure at ", which(unlist(x) == F), " ..."), "\n") 
} 


cat("TESTING PSI","\n")
nm <- model_data$nm
for(i in 1:model_data$nf){
  for(t in 1:model_data$k){
    if(is.na(model_data$apf[i,t])) next
    if(model_data$apf[i,t] == nm+1 & sum(model_data$psi[i,1:nm,t]) == 0){
      next
    } else if(model_data$apf[i,t] == nm+1 & sum(model_data$psi[i,1:nm,t]) > 0){
      stop(paste0("single female with partner options at time", t, "and index i:", i))
    }
    
    j <- model_data$apf[i,t]
    
    if(j != which(model_data$psi[i,,t]==1)|sum(model_data$psi[i,,t]) > 1) print(paste0("Female i:", i, "time ", t, " mismatch partner j:", j))
  }
}


cat("TESTING RECAPTURE FLOW MATING...","\n")
yf <- list()
ym <- list()

qf <- list()
qm <- list()

for(t in 1:model_data$k){
  
  index_f <- which(model_data$recap_f[1:nf,t] == 1)
  index_m <- which(model_data$recap_m[1:nm,t] == 1)
  
  index_q_f <- which(model_data$recap_f[1:nf,t] != 1)
  index_q_m <- which(model_data$recap_m[1:nm,t] != 1)
  
  yf[[t]] <- (all(model_data$amating_f[index_f,t] == model_data$mating_f[index_f,t]))
  ym[[t]] <- (all(model_data$amating_m[index_m,t] == model_data$mating_m[index_m,t]))  
  
  qf[[t]] <-  (all(is.na(model_data$amating_f[index_q_f,t])))
  qm[[t]] <-  (all(is.na(model_data$amating_m[index_q_m,t])))
} 

if(!all(unlist(yf))) cat(paste0("yf-test Failure at ", which(unlist(yf) == F), " ..."), "\n") 
if(!all(unlist(ym))) cat(paste0("ym-test Failure at ", which(unlist(ym) == F), " ..."), "\n") 
if(!all(unlist(qf))) cat(paste0("qf-test Failure at ", which(unlist(qf) == F), " ..."), "\n") 
if(!all(unlist(qm))) cat(paste0("qm-test Failure at ", which(unlist(qm) == F), " ..."), "\n") 


cat("TESTING RECAPTURE FLOW SURVIVAL","\n")
yf <- list()
ym <- list()

qf <- list()
qm <- list()

for(t in 1:model_data$k){
  
  nf <- model_data$nf
  nm <- model_data$nm
  
  index_f <- which(model_data$recap_f[1:nf,t] == 1)
  index_m <- which(model_data$recap_m[1:nm,t] == 1)
  
  index_q_f <- which(model_data$recap_f[1:nf,t] != 1)
  index_q_m <- which(model_data$recap_m[1:nm,t] != 1)
  
  yf[[t]] <- (all(model_data$sf[index_f,t] == model_data$af[index_f,t]))
  ym[[t]] <- (all(model_data$sm[index_m,t] == model_data$am[index_m,t]))  
  
} 

if(!all(unlist(yf))) cat(paste0("yf-test Failure at ", which(unlist(yf) == F), " ..."), "\n") 
if(!all(unlist(ym))) cat(paste0("ym-test Failure at ", which(unlist(ym) == F), " ..."), "\n") 
if(!all(unlist(qf))) cat(paste0("qf-test Failure at ", which(unlist(qf) == F), " ..."), "\n") 
if(!all(unlist(qm))) cat(paste0("qm-test Failure at ", which(unlist(qm) == F), " ..."), "\n") 
