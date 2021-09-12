library(bayess)
data(eurodip)
set.seed(1)
eurodip[eurodip > 0] <- 1
cjs_dat2 <- list(x = eurodip, n = nrow(eurodip), k = ncol(eurodip))
a <- matrix(NA, nrow = nrow(eurodip), ncol = ncol(eurodip))
initial_entry <- rep(NA,cjs_dat2$n)

for(i in 1:nrow(eurodip)){
  xi <- cjs_dat2$x[i,]
  last_seen <- max(which(xi==1))
  first_seen <- min(which(xi==1))
  
  a[i,first_seen:last_seen] <- 1
  initial_entry[i] <- first_seen 
}

cjs_dat2$a <- a
cjs_dat2$initial_entry <- initial_entry
cjs_dat2$female <- rbinom(cjs_dat2$n, 1, 0.5)
