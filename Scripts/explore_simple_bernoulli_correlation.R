#Inline Paste0 function
`%+%` <- function(a, b) paste0(a, b)

#Odds function
odds <- function(p){
  out <- p/(1-p)
  return(out)
}

# odds ratio
or <- function(p,q){
  odds_p <- odds(p)
  odds_q <- odds(q)
  out <- odds_p/odds_q
  return(out)
}

# odds product
op <- function(p,q){
  odds_p <- odds(p)
  odds_q <- odds(q)
  out <- odds_p*odds_q
  return(out)
}

#Extract joint binomial parameters
compute_jbin_param <- function(prob.f,prob.m){
  #Compute Components
  prob.prod  <- prob.m * prob.f
  sig.prob.f <- sqrt(prob.f*(1-prob.f))
  sig.prob.m <- sqrt(prob.m*(1-prob.m))
  sig.prod   <- sig.prob.f*sig.prob.m
  ##Upper and Lower bounds for cov of joint binomial
  lub        <-  sqrt(pmin(or(prob.f,prob.m),1/or(prob.f,prob.m)))
  glb        <- -sqrt(pmin(op(prob.f,prob.m),1/op(prob.f,prob.m)))
  #Return values
  out <- list(prob.prod       = prob.prod,
              sig.prob.f      = sig.prob.f,
              sig.prob.m      = sig.prob.m,
              sig.prod        = sig.prod,
              cor_lower_bound = glb,
              cor_upper_bound = lub)
  return(out)
}

#Generate Derived Probabilities
compute_jbin <- function(prob.f,prob.m,corr){
  
  parameters      <- compute_jbin_param(prob.f,prob.m)
  cor_upper_bound <- parameters[["cor_upper_bound"]]
  cor_lower_bound <- parameters[["cor_lower_bound"]]
  sig.prob.f      <- parameters[["sig.prob.f"]]
  sig.prob.m      <- parameters[["sig.prob.m"]]
  
  # If we are on a boundary condition update correlation structure 
  boundary_f                <-  which(prob.f == 1| prob.f == 0)
  boundary_m                <-  which(prob.m == 1| prob.m == 0)
  boundary                  <-  sort(unique(c(boundary_f, boundary_m)))
  corr[boundary]            <- 0 #correlation is zero if 1 or 0 for both probs
  cor_upper_bound[boundary] <- 0
  cor_lower_bound[boundary] <- 0
  
  if(any(round(corr,3) > round(cor_upper_bound,3)) || any(round(corr,3) < round(cor_lower_bound,3))){
    stop("Correlation " %+% corr %+% " is outside of the bounds [" %+% 
           as.character(cor_lower_bound) %+% " , " %+%
           as.character(cor_upper_bound) %+% "]")
  }
  
  ###Joint Probability Distn 
  prob.mf <- corr * sig.prob.m * sig.prob.f + (prob.f*prob.m)
  prob.f0 <- prob.f - prob.mf
  prob.m0 <- prob.m - prob.mf
  prob.00 <- 1 - prob.f0 - prob.m0 - prob.mf
  
  #List of parameters
  out <- list(prob.mf = prob.mf,
              prob.f0 = prob.f0,
              prob.m0 = prob.m0,
              prob.00 = prob.00)
  
  #Return values
  return(out)
}

generate_correlated_binomial <- function(n,p1,p2,rho){
  x <- rbinom(n,1,p1)
  joint_dens <- compute_jbin(p1,p2,rho)
  p2_cond0 <- joint_dens$prob.mf/p1
  p2_cond1 <- joint_dens$prob.m0/(1-p1)
  p2_cond <- ifelse(x == 0, p2_cond1, p2_cond0)
  y <- rbinom(n,1,p2_cond)
  
  return(data.frame(x1=x,x2=y))
  
}

ll_fn <- function(pars,
                  x1,
                  x2,
                  p1,
                  p2){
  
  rho <- pars[1]
  
  joint_dens <- compute_jbin(p1,p2,rho)
  p2_cond0 <- joint_dens$prob.mf/p1
  p2_cond1 <- joint_dens$prob.m0/(1-p1)
  p2_cond <- ifelse(x1 == 0, p2_cond1, p2_cond0)
  
  ll1 <- dbinom(x = x1, size = 1, prob = p1, log = T)
  ll2 <- dbinom(x = x2, size = 1, prob = p2_cond, log = T)
  return(-sum(ll1+ll2))
}


n   <- 30
p1  <- 0.75
p2  <- 0.65
rho <- -0.4

df <- generate_correlated_binomial(n, p1,p2, rho)

lower <- compute_jbin_param(p1,p2)$cor_lower_bound
upper <- compute_jbin_param(p1,p2)$cor_upper_bound


nlminb(start          = runif(1,min = lower,max = upper),
       objective       = ll_fn,
       x1              = df$x1,
       x2              = df$x2,
       p1              = p1,
       p2              = p2,
       lower           = lower,
       upper           = upper)

rho_mesh <- seq(lower, upper, 0.001)

ll_mesh <- sapply(rho_mesh, function(r) ll_fn(r, x1 = df$x1, x2 = df$x2, p1 = p1, p2 = p2)) 

plot(y = -ll_mesh, x = rho_mesh, type = "l")


compute_survival_correlation <- function(x,
                                         phiF, 
                                         phiM,
                                         pMF){

  sigF <- sqrt(phiF * (1-phiF))
  sigM <- sqrt(phiM * (1-phiM))

  gamma <- ((x/pMF) - phiM * phiF)/(sigF * sigM)
  
  gl <- compute_jbin_param_cjs(phiF,phiM)$cor_lower_bound
  gu <- compute_jbin_param_cjs(phiF,phiM)$cor_upper_bound
  
  
  
  return(pmin(gu, pmax(gl, gamma)))
}


sum(1 * (ps_data$apairs_f != 211),na.rm = T)/prod(dim(ps_data$apairs_f))
ps_data <- sim_dat(param_list) # pair-swap data
N <- 0
n <- 0

nf <- ps_data$nf
nm <- ps_data$nm
first_capture_f <- ps_data$first_capture_f
k <- ps_data$k
apairs_f <- ps_data$apairs_f

for(i in 1:nf){
  for(t in (first_capture_f[i]+1):k){
    if(!is.na(apairs_f[i,t]) & !is.na(apairs_f[i,t-1])){
      if(apairs_f[i,t] != (nm+1) & (apairs_f[i,t] == apairs_f[i,t-1])){
        n <- n+1
        N <- N+1
      }
    } else {
      if(!is.na(apairs_f[i,t-1])){
        if(apairs_f[i,t-1] != (nm+1)){
          N <- N+1
        }
      }
    }
  }
}

compute_survival_correlation(n/N, 0.8,0.8,compute_jbin_cjs(0.75,0.75,0.5)$prob.mf)
