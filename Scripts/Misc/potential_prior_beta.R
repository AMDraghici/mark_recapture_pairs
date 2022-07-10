

dbeta_corr <- function(x, par1, par2){
  tx <- (x + 1)/2
  beta_raw <- pbeta(tx,par1,par2) - pbeta(tx-0.001,par1,par2)
  return(2*beta_raw - 1)
}

dbeta_corr_trunc <- function(x, par1, par2, lb, ub){
  
  beta_corr_raw <- dbeta_corr(x, par1, par2)
  cdf_beta_corr_lower <- sum(dbeta_corr(seq(-1,lb, 0.001), par1, par2))
  cdf_beta_corr_upper <- sum(dbeta_corr(seq(-1,ub, 0.001), par1, par2))
  
  
  
  beta_corr_trunc_unnorm <- beta_corr_raw * (x >= lb & x <= ub)
  beta_corr_trunc_unnorm/(cdf_beta_corr_upper-cdf_beta_corr_lower)
  
}


generalized_pbeta <- function(x, par1, par2, b, a){
  tx <- (x - a)/(b-a) 
  return(pbeta(tx, par1, par2))
  
}

generalized_dbeta <- function(x, par1, par2, b, a){
   
  return(generalized_pbeta(x, par1, par2,b,a) - generalized_pbeta(x-0.001, par1, par2,b,a))
  
}

generalized_dtbeta <- function(x, par1, par2, b, a, lb, ub){
  fy <- generalized_dbeta(x, par1, par2, b, a) * (x>=lb) * (x<=ub)
  cdf_lb <-  generalized_pbeta(lb, par1, par2, b, a)
  cdf_ub <-  generalized_pbeta(ub, par1, par2, b, a)
  norm_fy <- cdf_ub - cdf_lb
  return(fy/norm_fy)
  
}

par1 <- 2
par2 <- 2
b <- 1
a <- -1
lb <- -0.09615 #-1/9 
ub <-0.90654  # 1.0
x <- seq(-1,1,0.001)


plot((generalized_dtbeta(x,par1,par2,b,a, lb, ub)), x = x, type = "l", col = "red", lty = 2)
lines((generalized_dtbeta(x,par1 = 10,par2 = 10,b,a, lb, ub)), x = x, type = "l", col = "blue", lty = 2)


lines((generalized_dbeta(x,par1,par2,b,a)), x = x, type = "l")

sum((generalized_dbeta(x,par1,par2,b,a)) * x)
sum((generalized_dtbeta(x,par1,par2,b,a, lb, ub)) * x)
