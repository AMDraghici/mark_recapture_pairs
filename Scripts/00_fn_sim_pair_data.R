## 1. Generic Functions -------------------------------------------------------------------------------------

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

#Logistic Function for probability
logit <- function(p){
  out <- log(odds(p))
  return(out)
}

# #Inverse Logistic Function
inv.logit <- function(x){
  out <- (1+exp(-x))^(-1)
  return(out)
}


# Numerically Stable Version
softmax <- function(par){
  n.par <- length(par)
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
}

## 2. S/R Distribution Functions -------------------------------------------------------------------------------------

#Extract joint binomial parameters
compute_jbin_param_cjs <- function(prob.f,prob.m){
  #Compute Components
  prob.prod <- prob.m * prob.f
  sig.prob.f <- sqrt(prob.f*(1-prob.f))
  sig.prob.m <- sqrt(prob.m*(1-prob.m))
  sig.prod <- sig.prob.f*sig.prob.m
  ##Upper and Lower bounds for cov of joint binomial
  lub <-  sqrt(pmin(or(prob.f,prob.m),1/or(prob.f,prob.m)))
  glb <- -sqrt(pmin(op(prob.f,prob.m),1/op(prob.f,prob.m)))
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
compute_jbin_cjs <- function(prob.f,prob.m,corr){
  #Extract Parameters
  parameters <- compute_jbin_param_cjs(prob.f,prob.m)
  cor_upper_bound <- round(parameters[["cor_upper_bound"]],3)
  cor_lower_bound <- round(parameters[["cor_lower_bound"]],3)
  corr <- round(corr, 3)
  sig.prob.f <- parameters[["sig.prob.f"]]
  sig.prob.m <- parameters[["sig.prob.m"]]
  
  # If we are on a boundary condition update correlation structure 
  boundary_f <-  which(prob.f == 1| prob.f == 0)
  boundary_m <-  which(prob.m == 1| prob.m == 0)
  boundary   <-  sort(unique(c(boundary_f, boundary_m)))
  corr[boundary] <- 0 #correlation technically wouldn't be zero but it makes the math work out this way
  cor_upper_bound[boundary] <- 0
  cor_lower_bound[boundary] <- 0
  
  if(any(corr > cor_upper_bound) || any(corr < cor_lower_bound)){
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


#Compute Survival from time t to t+1 for an entity  
survival <- function(previous_state, j, phi.m, phi.f, gam){
  
  # Compute Probability Distribution
  prob.distn <- compute_jbin_cjs(prob.m = phi.m, prob.f = phi.f, corr = gam)
  phi.mf <- prob.distn[["prob.mf"]]
  phi.f0 <- prob.distn[["prob.f0"]]
  phi.m0 <- prob.distn[["prob.m0"]]
  phi.00 <- prob.distn[["prob.00"]]
  
  # Construct Prob Matrix
  phi.matrix <- matrix(c(1,            0,         0,         0,
                         (1-phi.f[j]), phi.f[j],  0,         0,
                         (1-phi.m[j]), 0,         phi.m[j],  0,
                         phi.00[j],    phi.f0[j], phi.m0[j], phi.mf[j]), 
                       nrow = 4, ncol = 4, byrow = TRUE)
  
  #Compute Survival Outcome
  next_state <- which(rmultinom(1,1,phi.matrix[previous_state,]) == 1)
  return(next_state)
} 

#Compute Survival at time t for an entity  
recapture <- function(current_state, j, p.m, p.f, rho){
  
  # Recover Probability Distribution
  prob.distn <- compute_jbin_cjs(prob.f = p.f, prob.m = p.m, corr =  rho)
  p.mf <- prob.distn[["prob.mf"]]
  p.f0 <- prob.distn[["prob.f0"]]
  p.m0 <- prob.distn[["prob.m0"]]
  p.00 <- prob.distn[["prob.00"]]
  
  # Construct Prob Matrix
  p.matrix <- matrix(c(1,         0,       0,       0,
                       (1-p.f[j]),p.f[j],  0,       0,
                       (1-p.m[j]),0,       p.m[j],  0,
                       p.00[j],   p.f0[j], p.m0[j], p.mf[j]), 
                     nrow = 4, ncol = 4, byrow = TRUE)
  
  # Draw a random recapture outcome
  obs <- which(rmultinom(1,1,p.matrix[current_state,]) == 1)
  return(obs)
} 

# 3. Mating Prob ----------------------------------------------------------------------------------------------- 

mate <- function(previous_state, j, delta){
  
  # Chose to mate at time t? Previous state represents survival at t-1 and is either 0 or 1 
  # Vectorized and can represent many animals
  obs <- rbinom(length(previous_state), 1, delta[j]*previous_state)
  #Return Vector (or scalar) mating status 
  return(obs)
}

# 4. Construct Model Matrices------------------------------------------------------------------------------------

# Build out data structures and do easy initialization 

# p-q split between males and females
construct_sexes <- function(n, prop.female = 0.5, random = F){
  # Percent split 
  prop.male <- 1 - prop.female
  #Compute Genders by random sample or fixed count with prob_female = prop.female
  if(random == T){
    # Randomly sample gender with 
    sex <- sort(ifelse(rbinom(n,1,prop.female)==1,"F","M"))
  } else { 
    # Exactly 100 - prop.female*100% of population is male and prop.female*100% is female
    sex <- sort(c(rep("F",n * prop.female), rep("M", n * prop.male)))
  }
  
  # Return vector of sex designations
  return(sex)
}

construct_init_entry <- function(n, k, random = F, inits = NULL){
  
  # Select initial entry randomly?
  if(random == T & is.null(inits)){
    initial_entry <- as.integer(sample(c(1:round(k-1)),size=n,replace=T))
  } else {
    # If not do we have pre-specified values?
    if(!is.null(inits)){
      initial_entry <- inits 
      # If not just do fixed intervals across n and k 
    } else {
      initial_entry <- rep(1:(k-1),ceiling(n/(k-1)))
      initial_entry <- initial_entry[1:n]
      initial_entry <- sample(initial_entry, n, F) # then shuffle them
    }
  }

  # Return initial entries
  return(initial_entry)
  
}

# Build recruitment matrix
construct_recruit <- function(n, k, sex){ n
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  
  #Female's recruits
  recruit_f <- matrix(NA, nrow = (nf+1), ncol = k)
  recruit_f[nf+1,1:k] <- 1  # dummy spots for single pairs
  
  #Male's recruits
  recruit_m <- matrix(NA, nrow = (nm+1), ncol = k)
  recruit_m[nm+1,1:k] <- 1 # dummy spots for single pairs
  
  #Pair recruits
  recruit <- array(NA, dim = c(nf+1, nm+1, k))
  recruit[nf+1,nm+1,1:k] <- 1
  
  # Group as list
  recruit_list <- list(recruit_f = recruit_f, 
                       recruit_m = recruit_m,
                       recruit   = recruit)
  
  # Return recruit data 
  return(recruit_list)
  
}

#Skeleton Matrix for mate choice at k 
construct_mated <- function(n, k, sex){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  
  #Female's mate status
  mating_f <- matrix(NA, nrow = (nf+1), ncol = k)
  mating_f[nf+1,1:k] <- 1  # dummy spots for single pairs
  
  #Male's mate status
  mating_m <- matrix(NA, nrow = (nm+1), ncol = k)
  mating_m[nm+1,1:k] <- 1 # dummy spots for single pairs
  
  #Pair mate status
  mating <- array(NA, dim = c(nf+1, nm+1, k))
  mating[1:(nf+1),nm+1,1:k] <- 1
  mating[(nf+1),1:(nm+1),1:k] <- 1
  
  # Group as list
  mating_list <- list(mating_f = mating_f, 
                      mating_m = mating_m,
                      mating   = mating)
  
  # Return mating data
  return(mating_list)
  
}

#Skeleton array for pair assignment at k
construct_pairs <- function(sex,k){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  
  #Female's partners
  pairs_f <- matrix(NA, nrow = (nf), ncol = k)
  #Male's partners
  pairs_m <- matrix(NA, nrow = (nm), ncol = k)
  #Pair Identities
  pairs   <- array(NA, dim = c(nf+1, nm+1, k))
  # List of objects
  pairs_list <- list(pairs_f = pairs_f, pairs_m = pairs_m, pairs = pairs)
  return(pairs_list)
}

# Event that previous pair reforms 
construct_repartner <- function(sex, k){
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  
  repartner <- matrix(NA, nrow = nf, ncol = k)
  repartner[,1] <- 0 
  return(repartner)
}

# Construct pair coefficients
construct_coef <- function(sex, k){
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  #Intercept
  intercept <- array(1, dim = c(nf+1,nm+1,k))
  # Number of times a pair occurred
  histories <- array(0, dim = c(nf+1,nm+1,k))
  coef_list <- list(intercept =intercept, histories = histories)
  return(coef_list)
}

# Survival matrices for males, females, and 
construct_survival <- function(sex, k){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  #Females
  sf <- matrix(NA, nrow = nf+1, ncol = k)
  sf[nf+1,1:k] <- 1 # dummy spots for single pairs
  #Males
  sm <- matrix(NA, nrow = nm+1, ncol = k)
  sm[nm+1,1:k] <- 1 # dummy spots for single pairs
  #Pairs (and pair-single/single-pair combinations)
  spair <- array(NA, dim = c(nf+1, nm+1, k))
  spair[nf+1,nm+1,1:k] <- 4
  # Group into list and return data object
  surv_matrices <- list(sf    = sf,
                        sm    = sm,
                        spair = spair)
  return(surv_matrices)
}

# Recapture Array for pairs (don't need separate ones as recapture is not latent)
construct_recapture <- function(sex, k){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  
  #Pairs (and pair-single/single-pair combinations)
  rpair <- array(NA, dim = c(nf+1, nm+1, k))
  rpair[nf+1,nm+1,1:k] <- 4
  
  return(rpair)
}

# 5. Initialize CR States---------------------------------------------------------------------------------------

initialize_recruit <- function(sex, k, initial_entry, recruit_f, recruit_m, recruit){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  n  <- length(sex)
  
  # initial entry for males and females
  initial_entry_female <- initial_entry[1:nf]
  initial_entry_male   <- initial_entry[(nf+1):n]
  
  # Recruitment females
  for(i in 1:nf){
    recruit_f[i,initial_entry_female[i]:k] <- 1
    if(initial_entry_female[i] > 1) recruit_f[i,1:(initial_entry_female[i]-1)] <- 0
  }
  
  # Recruitment females
  for(i in 1:nm){
    recruit_m[i,initial_entry_male[i]:k] <- 1
    if(initial_entry_male[i] > 1) recruit_m[i,1:(initial_entry_male[i]-1)] <- 0
  }
  
  # Joint recruitment status
  for(i in 1:k){
    recruit[1:(nf+1),1:(nm+1),i] <- recruit_f[1:(nf+1),i] %*% t(recruit_m[1:(nm+1),i])
  }
  
  #Recruit List
  recruit_list <- list(recruit_f = recruit_f, 
                       recruit_m = recruit_m,
                       recruit   = recruit)
  
  # Return recruit data 
  return(recruit_list)
}


initialize_mating_choice <- function(n, k, initial_entry, delta, mating_f, mating_m, mating, sex){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  n <- length(sex)
  
  # initial entry for males and females
  initial_entry_female <- initial_entry[1:nf]
  initial_entry_male <- initial_entry[(nf+1):n]
  
  # Alive at initial entry
  previous_state <- rep(1, n)
  # Choose whether to mate or not
  initial_mating_choice <- mate(previous_state,initial_entry,delta)
  
  # Mating by females
  for(i in 1:nf){
    mating_f[i,initial_entry_female[i]] <- initial_mating_choice[i]
    if(initial_entry_female[i] > 1) mating_f[i,1:(initial_entry_female[i]-1)] <- 0
  }
  
  #Mating by males
  for(i in 1:nm){
    mating_m[i,initial_entry_male[i]] <- initial_mating_choice[i+nf]
    if(initial_entry_male[i] > 1) mating_m[i,1:(initial_entry_male[i]-1)] <- 0
  }
  
  #Joint mating distribution
  for(t in 1:k){
    mating[1:nf,1:nm,t] <- mating_f[1:nf,t] %*% t(mating_m[1:nm,t])
  }
  
  # Group results and join
  mating_list <- list(mating_f = mating_f,
                      mating_m = mating_m,
                      mating   = mating)
  
  return(mating_list)
}


initialize_partner_status <- function(n, coef_list, pairs, mating, recruit, initial_entry, sex){
  
  # When did the first animal get spotted
  time <- min(initial_entry) 
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  
  # Raw probability of mating (without conditions)
  psi_t <- matrix(1, nrow = nf+1, ncol = nm + 1) #1 + 0 * compute_partner_probs(sex, time, coef_list, c(1,0))

  # Flat probability of partnerships
  prob_mate <- psi_t * mating[,,time] * recruit[,,time]
  
  for(j in 1:(nf)){

    # Probability of partnerships forming for j -> 1:n
    probj <- prob_mate[j,]
    
    # Any pairs that have formed cannot form again
    if(j == 2){
      probj <- probj*c((1-pairs[1,1:nm,time]),1) #Remove pair j=1 from the prob (colsum only works on matrices)
    } else if(j > 2){
      probj <- probj*c((1-colSums(pairs[1:(j-1),1:nm,time])),1) #remove all pre-formed pairs
    }
    
    if(sum(probj[1:nm])==0){ 
      probj[nm+1] <- 1
    } else {
      probj[nm +1] <- 0
    }
    
    # Draw partnership
    pairs[j,,time] <- t(rmultinom(1,1,probj)) # assign a partner 
  }
  
  # Assign unmated males to singles row
  pairs[nf+1,,time] <- c(1-colSums(pairs[1:nf,1:nm,time]),1)
  
  # Return updated pairs matrix
  return(pairs)
}


initialize_survival_status <- function(sex, k, initial_entry, sf, sm){
  
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  n  <- length(sex)
  
  # Survival females
  # Assign survival status to initial entry
  for(i in 1:(k-1)){
    
    # Survival females
    enter_at_i_female <- which(initial_entry[1:nf] == i)
    sf[enter_at_i_female,i] <- 1
    
    #survival males 
    enter_at_i_male <- which(initial_entry[(nf+1):n] == i)
    sm[enter_at_i_male,i] <- 1
    
    # Turn the previous entries into 1 for likelihood calculation
    if(i > 1){
      sf[enter_at_i_female,1:(i-1)] <- 0
      sm[enter_at_i_male,  1:(i-1)] <- 0
    }
  }
  
  # Return Survival states
  init_surv_list <- list(sf = sf, sm = sm)
  return(init_surv_list)
  
}

# 6. Update Data structures from [t,t+1)------------------------------------------------------------------------

# Assign mate id and survival id to male/female individual level vectors
propogate_partner_state <- function(pairs, sex, pairs_f, pairs_m, time){

  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  n <- length(sex)
  
  # Recover index entries
  pairs_f[,time] <- pairs[1:nf,,time] %*% 1:(nm+1)
  pairs_m[,time] <- t((1:(nf+1)) %*% pairs[,1:nm,time])
  # Return updated information
  return(list(pairs_f = pairs_f,pairs_m = pairs_m))
}

# Convert joint survival to marginal states
propogate_surv_individual <- function(sf, sm, spair, time, sex){
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  
  surv_marginal_female <- c(0,1,0,1)
  surv_marginal_male <- c(0,0,1,1)
  
  # Only update females (not dummy states)
  for(i in 1:nf){
    female_state <- surv_marginal_female[max(spair[i,,time])]
    sf[i,time] <- female_state
  }
  
  # Only update males (not dummy states)
  for(i in 1:nm){
    male_state <- surv_marginal_male[max(spair[,i,time])]
    sm[i,time] <- male_state
  }
  
  # Group and return results 
  ind_surv_list <- list(sf = sf, sm = sm)
  return(ind_surv_list)
}


# Convert marginal recapture states (at the end)
propogate_recap_individual <- function(sex, k, rpair, recruit_f, recruit_m){
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  n <- length(sex)
  
  recap_marginal_female <- c(0,1,0,1)
  recap_marginal_male <- c(0,0,1,1)
  
  # Recapture individual level male and female 
  recap_f <- matrix(NA, nrow = nf+1, ncol = k)
  recap_f[(nf + 1),] <- 0
  recap_m <- matrix(NA, nrow = nm+1, ncol = k)
  recap_m[(nm + 1),] <- 0
  
  # Only update females (not dummy states)
  for(i in 1:nf){
    for(time in 1:k){
      female_observation <- recap_marginal_female[max(rpair[i,,time])]
      recap_f[i,time] <- female_observation
    }
    
    # HACKY ADJUSTMENT FOR WHEN AN INDIVIDUAL IS NEVER CAUGHT
    # If they're never seen force an observation at initial entry
    # Should introduce negligble bias 
    if(sum(recap_f[i,]) == 0){recap_f[i, min(which(recruit_f[i, ] == 1))] <- 1}
  }
  
  # Only update males (not dummy states)
  for(j in 1:nm){
    for(time in 1:k){
      male_observation <- recap_marginal_male[max(rpair[,j,time])]
      recap_m[j,time] <- male_observation
    }
    
    # HACKY ADJUSTMENT FOR WHEN AN INDIVIDUAL IS NEVER CAUGHT
    # If they're never seen force an observation at initial entry
    # Should introduce negligble bias 
    if(sum(recap_m[j,]) == 0){recap_m[j, min(which(recruit_m[j, ] == 1))] <- 1}
    
  }
  
 
  
  
  
  
  # Group and return results 
  ind_recap_list <- list(recap_f = recap_f, recap_m = recap_m)
  return(ind_recap_list)
}

# Convert marginal to joint survival states
propogate_surv_pairs <- function(sex, sf, sm, spair, time, pairs_f, pairs_m){
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  n <- length(sex)
  
  # Covers male/female pairs and female singles
  for(i in 1:nf){
    j <- pairs_f[i,time]
    
    state_vector <- c(
      (1-sf[i,time])*(1-sm[j,time]), #both dead
      sf[i,time]*(1-sm[j,time]), # female alive
      (1-sf[i,time])*sm[j,time], # male alive
      sf[i,time]*sm[j,time] # both alive 
    )
    
    status <- which(state_vector == 1)
    
    if(j == nm+1){
      state_value <- ifelse(status == 4, 2, ifelse(status == 3, 1, status))
    } else {
      state_value <- status
    }
    
   
    spair[i,j,time] <- state_value
  }
  
  for(j in 1:nm){
    i <- pairs_m[j,time]
    
    state_vector <- c(
      (1-sf[i,time])*(1-sm[j,time]), #both dead
      sf[i,time]*(1-sm[j,time]), # female alive
      (1-sf[i,time])*sm[j,time], # male alive
      sf[i,time]*sm[j,time] # both alive 
    )
    
    status <- which(state_vector == 1)
    
    if(i == (nf+1)){
      state_value <- ifelse(status == 4, 3, ifelse(status == 2, 1, status))
    } else {
      state_value <- status
    }
    
    spair[i,j,time] <- state_value
  }
  
  #Remaining entries are "dead"
  spair[,,time][is.na(spair[,,time])] <- 1
  
  # Return joint survival density
  return(spair)
  
}

compute_partner_lodds <- function(sex, time, coef_list, betas){
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  
  # Assign memory 
  eta_t <- matrix(0, nrow = nf+1, ncol = nm+1)
  
  beta_vec <- unlist(betas)
  
  # Add all the vectorized coefficients
  for(i in 1:length(beta_vec)){
    eta_t <- eta_t + coef_list[[i]][,,time]*beta_vec[i]
  }
  
  # Return log-odds of partnerships at time t
  return(eta_t)
}


compute_partner_probs <- function(sex, time, coef_list, betas){
  
  # Grab log-odds of pairs forming in order
  eta_t <- compute_partner_lodds(sex, time, coef_list, betas)
  
  psi_t <- matrix(NA, nrow = nrow(eta_t), ncol = ncol(eta_t))
  
  # convert to probability using softmax 
  for(x in 1:nrow(psi_t)){
    psi_t[x,] <- softmax(eta_t[x,])
  }
 
  
  #Return raw probability of mating (unconditional)
  return(psi_t)
}

# Update History
update_history <- function(coef_list, pairs, time, sex){
  
  # Number of females and males (ignore singles)
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  
  # Add partnerships at time t-1 to current partnership history at t
  coef_list[["histories"]][1:nf,1:nm,time] <- coef_list[["histories"]][1:nf,1:nm,time-1] +
    pairs[1:nf,1:nm,time-1]
  # Return coefficients 
  return(coef_list)
}

# compute mating decision at t
compute_mating <- function(sex, time, delta, recruit, recruit_f, recruit_m,
                           mating_f, mating_m, mating, sf, sm){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  
  # Mating by females
  for(i in 1:nf){
    
    # Initial Entry
    init_f_i <- min(which(recruit_f[i,]==1))
    
    # If just entered then initialize function has assigned status 
    if(init_f_i >= time) next
    
    # Otherwise update status 
    mating_f[i,time] <- mate(sf[i,time-1], time, delta)
    
  }
  
  # Mating by Males
  for(i in 1:nm){
    
    # Initial Entry
    init_m_i <- min(which(recruit_m[i,]==1))
    
    # If just entered then initialize function has assigned status 
    if(init_m_i >= time) next
    
    # Otherwise update status 
    mating_m[i,time] <- mate(sm[i,time-1], time, delta)
    
  }
  
  #Joint mating distribution
  mating[1:nf,1:nm,time] <- mating_f[1:nf,time] %*% t(mating_m[1:nm,time])
  
  # Group results and join
  mating_list <- list(mating_f = mating_f,
                      mating_m = mating_m,
                      mating = mating)
  
  return(mating_list)
}

# Does pair ij from t-1 reform with probability inv.logit(beta0 + beta1*history[i,j])
compute_re_partner <- function(repartner, sex, coef_list, betas, pairs, mating, time){
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  
  # Compute probability of reforming pairs
  for(i in 1:nf){
    previous_partner_i <- which(pairs[i,,time-1]==1) # Previous mate 
    history_ij <- coef_list[[2]][i,previous_partner_i,time] # number of times this pair formed
    prob <- inv.logit(c(1.0,history_ij) %*% unlist(betas)) # probability of reforming
    prob <- prob * mating[i,previous_partner_i,time] # both have to be willing to mate to repair
    prob <- prob * (1-pairs[i,nm+1,time-1]) # if previously single then set to zero (st prob of new mate isnt biased downward)
    repartner[i,time] <- rbinom(1,1,prob = prob) #  simulate re-pair status
  }
  return(repartner)
}


# compute partnerships at t

compute_partnerships <- function(sex, coef_list, pairs, mating, time, repartner){
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  
  # Raw probability of mating (without conditions)
  psi_t <- matrix(1, nrow = nf+1, ncol = nm + 1) # + 0 * compute_partner_probs(sex, time, coef_list, c(1,0))
  
  # Probability of partnerships with mating and recruitment included (surv is baked into mating)
  prob_mate <- psi_t * mating[,,time]

  male_taken <- rep(0,nm)
  
  for(j in 1:nm){
    male_taken[j] <- sum(pairs[1:nf,j,time-1] * repartner[1:nf, time])
    if(any(male_taken > 1|male_taken <0)) browser()
  }

  for(i in 1:nf){
    for(j in 1:nm){
      prob_mate[i,j] <- prob_mate[i,j] *(1-pairs[i,j,time-1])* (1 - repartner[i,time]) *(1-male_taken[j])  + 
        pairs[i,j,time-1] * repartner[i,time]
    }
  }
  
  for(j in 1:nf){
    # Probability of partnerships forming for j -> 1:n
    # If repartner = 1 then a pair reforms 
    probj <- prob_mate[j,]
    
    # Any pairs that have formed cannot form again
    if(j == 2){
      probj <- probj*c((1-pairs[1,1:nm,time]),1) #Remove pair j=1 from the prob (colsum only works on matrices)
    } else if(j > 2){
      probj <- probj*c((1-colSums(pairs[1:(j-1),1:nm,time])),1) #remove all pre-formed pairs
    }
    
    if(sum(probj[1:nm])==0){ 
      probj[nm+1] <- 1
    } else {
      probj[nm+1] <- 0
    }
    
    if(any(probj < 0)) browser()
    
    # If repartnering just set to previous partner
    if(repartner[j,time] ==1){
      pairs[j,,time] <- 0
      pairs[j,which(pairs[j,,time-1]==1),time] <- 1
    } else {
      # Draw partnership
      pairs[j,,time] <- t(rmultinom(1,1,probj)) # assign a partner 
    }
    
  }
  
  # Assign unmated males to singles row
  pairs[nf+1,,time] <- c(1-colSums(pairs[1:nf,1:nm,time]),1)
  
  # Return updated pairs matrix
  return(pairs)
}

# compute_partnerships <- function(sex, coef_list, betas, pairs, mating, recruit, time){
#   
#   nf <- length(sex[sex == "F"]) # Number of females
#   nm <- length(sex[sex == "M"]) # Number of males
#   
#   # Raw probability of mating (without conditions)
#   psi_t <- compute_partner_probs(sex, time, coef_list, betas)
#   
#   # Probability of partnerships with mating and recruitment included (surv is baked into mating)
#   prob_mate <- psi_t * mating[,,time] * recruit[,,time]
#   
#   for(j in 1:nf){
#     # Probability of partnerships forming for j -> 1:n
#     probj <- prob_mate[j,]
#     
#     # Any pairs that have formed cannot form again
#     if(j == 2){
#       probj <- probj*c((1-pairs[1,1:nm,time]),1) #Remove pair j=1 from the prob (colsum only works on matrices)
#     } else if(j > 2){
#       probj <- probj*c((1-colSums(pairs[1:(j-1),1:nm,time])),1) #remove all preformed pairs
#     }
#     
#     # Draw partnership
#     pairs[j,,time] <- t(rmultinom(1,1,probj)) # assign a partner 
#   }
#   
#   # Assign unmated males to singles row
#   pairs[nf+1,,time] <- c(1-colSums(pairs[1:nf,1:nm,time]),1)
#   
#   # Return updated pairs matrix
#   return(pairs)
# }

# Compute survival state at t
compute_survival <- function(spair, sf, sm, pairs_f, pairs_m, recruit_f, recruit_m, time, phi.m, phi.f, gam){
  
  # Survival for pairs and single females
  for(i in 1:nrow(pairs_f)){
    
    # Identify i's partner
    partner_i <- pairs_f[i, time]
    single_check <- (partner_i == (nrow(pairs_m) + 1))
    
    # Is initial entry now or later?
    init_f_i <- min(which(recruit_f[i,]==1))
    init_m_i <- min(which(recruit_m[partner_i,]==1))
    
    if(init_f_i > time & init_m_i > time) next
    
    # Previous survival states
    surv_i <- sf[i, time-1]*(init_f_i < time) + 1*(init_f_i == time)
    surv_m <- (sm[partner_i, time-1]*(init_m_i < time) + 1*(init_m_i == time))*(1-single_check)
    
    # Recover joint survival state for this new couple at t-1
    previous_state_i <- 1 + surv_i + 2*surv_m
    
    if(previous_state_i == 1) next
    
    # If you've entered now then you're alive at t so probability becomes 1 
    phiF <- phi.f
    phiM <- phi.m 
    
    phiF[time-1] <- phi.f[time-1]*(init_f_i < time) + 1*(init_f_i == time)   
    phiM[time-1] <- (phi.m[time-1]*(init_m_i < time) + 1*(init_m_i == time)) * (1-single_check)
    
    # Compute joint survival at time t
    spair[i,partner_i, time] <- survival(previous_state = previous_state_i,
                                         j = time-1,
                                         phi.m = phiM, 
                                         phi.f = phiF, 
                                         gam = gam)
  }
  
  # Survival for single males
  single_males <- which(pairs_m[,time]==(nrow(pairs_f)+1))
  
  for(j in single_males){
    
    # Is initial entry now or later?
    init_m_j <- min(which(recruit_m[j,]==1))
    
    if(init_m_i > time) next
    
    # Previous survival states
    surv_j <- sm[j, time-1]*(init_m_j < time) + 1*(init_m_j == time)
    
    # Recover joint survival state for this new couple at t-1
    previous_state_j <- 1 + 2*surv_j
    
    # If you've entered now then you're alive t so probability becomes 1 
    phiM <- phi.m
    
    phiM[time-1] <- phi.m[time-1]*(init_m_j < time) + 1*(init_m_j == time)   
    
    # Compute joint survival at time t
    spair[nrow(pairs_f)+1, j, time] <- survival(previous_state = 
                                                previous_state_j,
                                                j = time-1,
                                                phi.m = phiM,
                                                phi.f = rep(0,length(phi.f)),
                                                gam = gam)
  }
  
  #Remaining entries are "dead"
  spair[,,time][is.na(spair[,,time])] <- 1
  
  # Return joint survival density
  return(spair)
}


# Compute recapture state at t
compute_recapture <- function(sex, rpair, spair, pairs_f, pairs_m, time, p.m, p.f, rho, recruit_f, recruit_m){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  
  # Compute capture status
  for(i in 1:nf){
    j <- pairs_f[i, time]
    single_check <- (j == (nrow(pairs_m) + 1))
    current_state <- spair[i,j, time]
    i_first <- min(which(recruit_f[i,] == 1))
    j_first <- min(which(recruit_m[j,] == 1))
    # if(i_first == time) pF <- rep(1, length(p.f)) else if(i_first > time) pF <- rep(0, length(p.f)) else pF <- p.f
    # if(j_first == time) pM <- rep(1, length(p.m)) else if(j_first > time) pM <- rep(0, length(p.m)) else pM <- p.m
    if(i_first > time) pF <- rep(0, length(p.f)) else pF <- p.f
    if(j_first > time) pM <- rep(0, length(p.m)) else pM <- p.m
    rpair[i,j, time] <- recapture(current_state = current_state, 
                                  j = time, 
                                  p.m = pM * (1-single_check), 
                                  p.f = pF, 
                                  rho = rho)
  }
  
  single_males <- which(pairs_m[,time]==(nf+1)) #which(spair[nf+1,1:nm,time] != 1)
  
  for(j in single_males){
    current_state <- spair[nf+1,j, time]
    j_first <- min(which(recruit_m[j,] == 1))
    #if(j_first == time) pM <- rep(1, length(p.m)) else if(j_first > time) pM <- rep(0, length(p.m)) else pM <- p.m
    if(j_first > time) pM <- rep(0, length(p.m)) else pM <- p.m
    
    
    rpair[nf+1, j, time] <- recapture(current_state = current_state,
                                      j = time,
                                      p.m = pM,
                                      p.f = rep(0, length(pM)), #p.f,
                                      rho = rho)
  }
  
  # Turn rest to unobserved 
  rpair[,,time][is.na(rpair[,,time])] <- 1
  
  # Return recapture-pairs matrix
  return(rpair)
}

compute_hidden_recruitment <- function(recruit_f, recruit_m, recap_f, recap_m){
  
  # First captures
  first_female <- sapply(1:(nrow(recap_f)-1), function(x) min(which(recap_f[x,] == 1)))
  first_male   <- sapply(1:(nrow(recap_m)-1), function(x) min(which(recap_m[x,] == 1)))
  
  
  for(i in 1:length(first_female)){
    
    if(first_female[i] == 1| first_female[i] == Inf) next
    
    recruit_f[i, 1:(first_female[i]-1)] <- NA
    
  }
  
  
  for(i in 1:length(first_male)){
    
    if(first_male[i] == 1| first_male[i] == Inf) next
    
    recruit_m[i, 1:(first_male[i]-1)] <- NA
    
  }
  
  return(list(recruit_f = recruit_f,
              recruit_m = recruit_m))

  
}


compute_hidden_survival <- function(pairs_f, pairs_m, rpair, spair, sf, sm, k, sex){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"])  
  n <- length(sex)
  
  # Produce inferred survival states
  af <- matrix(NA, nrow = nrow(sf), ncol = ncol(sf))
  am <- matrix(NA, nrow = nrow(sm), ncol = ncol(sm))
  
  # States based on observations
  female_state <- c(NA, 1,  NA, 1) 
  male_state   <- c(NA, NA, 1,  1)
  pair_state <- c(NA, NA, NA, 4)
  
  # Add observations directly from joint recapture matrix (pairs + single females)
  for(i in 1:nf){
    for(time in 1:k){
      
      # Joint Recapture of pair i and pf[i,t] at time t
      x_ijt <- rpair[i, pairs_f[i,time] , time]
      
      # Uncouple survival states
      af[i, time] <- female_state[x_ijt]
      am[pairs_f[i,time], time] <- male_state[x_ijt]
    }
  }
  
  # Add observations directly from joint recapture matrix (single males)
  for(j in 1:nm){
    for(time in 1:k){
      if(is.na(pairs_m[j, time])) browser()
      # skip if its not a single male
      if(pairs_m[j,time] != (nf+1)) next
      
      # Joint Recapture of pair i and pf[i,t] at time t
      x_ijt <- rpair[(nf+1), j , time]
      
      # Uncouple survival states
      am[j, time] <- male_state[x_ijt]
    }
  }
  
  
  # Go back and add inferred states based on time
  
  # Females
  for(i in 1:nf){
    
    # If all unknown skip
    if(all(is.na(af[i,]))) next
    
    # Last time seen alive
    last_alive <- max(which(af[i,]==1))
    
    # If the last time they were seen alive was time 1 then skip 
    if(last_alive == 1) next 
    
    # Change all other states prior to last alive to alive
    af[i, 1:last_alive] <- 1    
    
    rm(last_alive)
  }
  
  #Males
  for(i in 1:nm){
    
    # If all unknown skip
    if(all(is.na(am[i,]))) next
    
    # Last time seen alive
    last_alive <- max(which(am[i,]==1))
    
    # If the last time they were seen alive was time 1 then skip 
    if(last_alive == 1) next 
    
    # Change all other states prior to last alive to alive
    am[i, 1:last_alive] <- 1    
    
    rm(last_alive)
  }
  
  # Must be alive at time 1 (if not recruited considered alive (or going to be alive) but not in pop)
  # This line is necessary because of the fact that some animals will be unobserved throughout the whole study//
  # and we need to force their first capture to 1
  af[,1] <- 1
  am[,1] <- 1
  
  # Set dummy states to known ones 
  af[(nf+1),] <- 1
  am[(nm+1),] <- 1
  
  # Store results in a list
  # add dummy row 1 for the repartner mechanism (in model index t is t+1)
  state_list <- list(af = af, #cbind(rep(1,nrow(af)),af), 
                     am = am) #cbind(rep(1,nrow(am)),am))
  
  # Return List
  return(state_list)
}

compute_hidden_pairs <- function(pairs_f, pairs_m, recap_f, recap_m, k, sex, repartner, mating_f, mating_m, show_unmated = F){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"])  
  n <- length(sex)
  
  # Produce inferred survival states
  apairs_f <- matrix(NA, nrow = nrow(pairs_f), ncol = ncol(pairs_f))
  apairs_m <- matrix(NA, nrow = nrow(pairs_m), ncol = ncol(pairs_m))
  apairs <- array(NA, dim =c(nf+1,nm+1,k+1))
  amating_f <- matrix(NA, nrow = nf+1, ncol = k)
  amating_m <- matrix(NA, nrow = nm+1, ncol = k) 
  arepartner <- matrix(NA, nrow = nf, ncol = k)

  # States based on observations
  pair_state <- c(NA, NA, NA, 1)
  
  # Add observations directly from joint recapture matrix
  
  # Females 
  for(i in 1:nf){
    for(time in 1:k){
      
      if(time >=2){
        x_ijt_minus_1 <- x_ijt 
      } else{
        x_ijt_minus_1 <- 4
      }
      
      # Joint Recapture of pair i and pf[i,t] at time t
      # x_ijt <- rpair[i, pairs_f[i,time] , time]
      x_ijt <- 1 + recap_f[i, time] + recap_m[pairs_f[i,time], time] * 2
      
      # Uncouple Pair Index
      apairs_f[i, time] <- pair_state[x_ijt] * pairs_f[i, time]
      amating_f[i, time] <- pair_state[x_ijt]
      
      # Observed Repartnership
      arepartner[i,time] <- pair_state[x_ijt] * pair_state[x_ijt_minus_1] * repartner[i, time]
      
      if(!is.na(pair_state[x_ijt]))  amating_m[ apairs_f[i, time], time] <- pair_state[x_ijt]
      
      
      # If we change the setting to show mating status of observed individuals 
      # We must update the case in which only females are observed
      # Then if she is not mating, or if her pair status is set to no partner
      # We update the data to reflect this
      if(show_unmated & x_ijt == 2){
        amating_f[i,time] <- mating_f[i, time]
        if(mating_f[i,time] == 0|pairs_f[i,time]==nm+1){
          arepartner[i,time] <- 0
          apairs_f[i, time] <- nm + 1
          # amating_f[i,time] <- 0
          # mating_f[i, time] <- 0
        }
      }
      
      
      
    }
  }
  
  # Males
  for(j in 1:nm){
    for(time in 1:k){
      
      # Joint Recapture of pair i and pf[i,t] at time t
      # x_ijt <- rpair[pairs_m[j, time], j , time]
      x_ijt <- 1 + recap_f[pairs_m[j, time], time] + recap_m[j, time] * 2
      
      # Uncouple Pair Index
      apairs_m[j, time] <- pair_state[x_ijt] * pairs_m[j, time]
      amating_m[j, time] <- pair_state[x_ijt]
      if(!is.na(pair_state[x_ijt]))  amating_f[ apairs_m[j, time], time] <- pair_state[x_ijt]
      
      # If we change the setting to show mating status of observed individuals 
      # We must update the case in which only males are observed
      # Then if he is not mating, or if his pair status is set to no partner
      # We update the data to reflect this
      if(show_unmated & x_ijt == 3){
        amating_m[j,time] <- mating_m[j, time]
        if(mating_m[j,time] == 0|pairs_m[j,time]==nf+1){
          apairs_m[j, time] <- nf + 1
          # amating_m[j,time] <- 0
          # mating_m[j, time] <- 0
        }
      }
      
      
    }
  }
  
  # Add details to joint density
  
  # All pairs (nf+1:n are dummy states)
  for(i in 1:nf){
    
    # If all unknown skip
    if(all(is.na(apairs_f[i,]))) next
    
    # We can only update known information
    time_index <- which(!is.na(apairs_f[i,]))
    
    # Add known pairings and designate known non-pairs
    for(time in time_index){
      # Assign states 
      apairs[i+1,,time+1] <- 0 # make all pairs unmated
      if(apairs_f[i,time] == nm + 1){next} # if single next 
      apairs[,apairs_f[i, time]+1, time+1] <- 0 # make all male partner pairs unavailable
      apairs[i+1,apairs_f[i, time]+1, time+1] <- 1 # assign pairing to the correct male/female combo at time t
      
    }
  }
  
  
  for(time in 2:k){
    
    # If all previous or current partners unknown for females then next, arepartner must be zero for whole slate 
    if(all(is.na(apairs_f[1:nf,time-1]))) next
    if(all(is.na(apairs_f[1:nf,time]))) next
    
    for(i in 1:nf){
      
      mask1 <- is.na(apairs_f[i,time-1]) # is previous partner unknown
      mask2 <- is.na(apairs_f[i,time])   # is current partner unknown
      
      if(mask1 & mask2){
        next  # if both are unknown repartner is unknown
      } else if(mask1 & !mask2){
        # Current partner is known, previous is unknown
        # Need to check if current partner had a different mate last time, if so then set to zero
        
        # Grab all known pairings from before
        previous_partners <- apairs_f[1:nf,time-1]  # Grab data 
        previous_partners <- previous_partners[!is.na(previous_partners)] # remove NA
        previous_partners <- previous_partners[previous_partners != (nm + 1)] # drop single designation 
        
        # Now check if current partner had a different mate at t-1
        current_partner                 <- apairs_f[i,time]
        current_single                  <- current_partner == (nm+1)
        current_partner_previous_status <-  any(previous_partners==current_partner) 
        
        # if match is found or single then this cant be a repartnership
        if(current_partner_previous_status|current_single) arepartner[i,time] <- 0
        
        
      } else if(!mask1 & mask2){
        # Current partner is unknown, previous is known
        # Need to check if previous has a new partner, if so then arepartner = 0
        
        # Grab known pairings at 
        current_partners <- apairs_f[1:nf,time]  # Grab data 
        current_partners <- current_partners[!is.na(current_partners)] # remove NA
        current_partners <- current_partners[current_partners != (nm + 1)] # drop single designation 
        
        # Check if previous partner has a different mate at t
        previous_partner <- apairs_f[i,time-1]
        previous_single <-  previous_partner == (nm+1)
        previous_partner_current_status <- any(current_partners==previous_partner) 
        
        # if match is found or single then this cant be a repartnership
        if(previous_partner_current_status|previous_partner_current_status) arepartner[i,time] <- 0
        
      } else if(!mask1 & !mask2){
        # Current and previous partner are known and not single then if both same arepartner = 1, otherwise arepartner = 0
        current_partner  <-  apairs_f[i,time]
        previous_partner <-  apairs_f[i,time-1]
        current_single   <-  current_partner == (nm+1)
        previous_single  <-  previous_partner == (nm+1)
        partner_match    <-  current_partner ==  previous_partner
        
        # Scenarios in which arepartner=0
        if(current_single|previous_single|!partner_match) arepartner[i,time] <- 0
        
        # Scenarios in which arepartner=1
        if(partner_match & !current_single & !previous_single) arepartner[i,time] <- 1
        
        # Otherwise just leave as NA
      }
      
    }
  }
  
  
  # Dummy states always available
  amating_f[(nf+1),1:k] <- 1
  amating_m[(nm+1),1:k] <- 1
  apairs[1:(nf+1),1,] <- 0
  apairs[1,1:(nm+1),] <- 0
  apairs[,,1] <- 0 
  arepartner[,1] <- 0 
    
  # # add dummy index (used to make JAGS/NIMBLE code cleaner)
  # apairs_f <- cbind(rep((nm+1),nrow(apairs_f)),apairs_f)
  # apairs_m <- cbind(rep((nf+1),nrow(apairs_m)),apairs_m)
  
  
  # Build index of possible pairings 
  # Used for homogenous pair assignment mechanism
  psi <- apairs
  psi[is.na(psi)] <- 1
  psi <- psi[1:nf+1,1:nm+1,]
  psi <- psi[,,1:k+1]
  
  psi_array <- array(NA,dim = c(nf,nm+1,k))
  
  for(t in 1:k){
    psi_array[,,t] <- cbind(psi[,,t],rep(0,nf))
  }
  
  # Store results in a list
  pairs_list <- list(apairs_m = apairs_m,
                     apairs_f = apairs_f,
                     apairs = apairs,
                     amating_f = amating_f,
                     amating_m = amating_m,
                     arepartner = arepartner,
                     psi = psi_array)
  
  # Return List
  return(pairs_list)
}

# 7. Simulate Data ---------------------------------------------------------------------------------------------

simulate_cr_data <- function(n,
                             k, 
                             prop.female,
                             delta,
                             phi.f, 
                             phi.m, 
                             gam, 
                             p.f, 
                             p.m, 
                             rho, 
                             betas, 
                             rand_init = T,
                             init = NULL,
                             show_unmated){
  
  # Make sure the number of individuals simulated is even
  if(!n %% 2 == 0){
    n <- n + 1
    if(!is.null(init)){
      init[n] <- sample(k-1,1)
    }
  }
  
  # Generate SKeleton Data Structures
  sex <- construct_sexes(n = n, prop.female = prop.female)
  initial_entry <- construct_init_entry(n = n, k = k, random = rand_init, init = init) 
  recruit_list <- construct_recruit(n = n, k = k, sex = sex)
  recruit_f <- recruit_list[["recruit_f"]]
  recruit_m  <- recruit_list[["recruit_m"]]
  recruit <- recruit_list[["recruit"]]
  mating_list <- construct_mated(n=n, k=k, sex = sex)
  mating_f <- mating_list[["mating_f"]] 
  mating_m <- mating_list[["mating_m"]]  
  mating <- mating_list[["mating"]]
  repartner <- construct_repartner(sex = sex, k = k)
  pairs_list <- construct_pairs(sex = sex, k = k)
  pairs_f <- pairs_list[["pairs_f"]]
  pairs_m <- pairs_list[["pairs_m"]]
  pairs <- pairs_list[["pairs"]]
  coef_list <- construct_coef(sex = sex, k = k)
  survival_list <- construct_survival(sex = sex, k = k)
  sf <- survival_list[["sf"]]
  sm <- survival_list[["sm"]]
  spair <- survival_list[["spair"]]
  rpair <- construct_recapture(sex = sex, k = k)
  
  # Initialize Data 
  initial_time <- min(initial_entry)
  recruit_list <- initialize_recruit(sex = sex, k = k, initial_entry = initial_entry, recruit_f = recruit_f, recruit_m = recruit_m, recruit = recruit)
  recruit_f <- recruit_list[["recruit_f"]]
  recruit_m <- recruit_list[["recruit_m"]]
  recruit <- recruit_list[["recruit"]]
  mating_list_init <- initialize_mating_choice(n = n, k = k, initial_entry = initial_entry,delta = delta, mating_f = mating_f, mating_m = mating_m,  mating = mating, sex = sex)
  mating_f <- mating_list_init[["mating_f"]] 
  mating_m <- mating_list_init[["mating_m"]]  
  mating <- mating_list_init[["mating"]]
  pairs <- initialize_partner_status(n = n, coef_list = coef_list, pairs = pairs, mating = mating, recruit = recruit, initial_entry = initial_entry, sex = sex)
  coef_list <- update_history(coef_list = coef_list, pairs = pairs, time = initial_time + 1, sex = sex)
  pairs_ind_list <- propogate_partner_state(pairs = pairs, sex =  sex, pairs_f = pairs_f, pairs_m = pairs_m, time = initial_time)
  pairs_f <- pairs_ind_list[["pairs_f"]]
  pairs_m <- pairs_ind_list[["pairs_m"]]
  init_surv_list <-  initialize_survival_status(sex = sex, k =  k, initial_entry = initial_entry, sf= sf, sm = sm)
  sf <- init_surv_list[["sf"]]
  sm <- init_surv_list[["sm"]]
  
  # Survival must happen at least on first capture, debug if not
  if(any(rowSums(sf, na.rm = T) < 1)) browser()
  if(any(rowSums(sm, na.rm = T) < 1)) browser()
  
  spair <- propogate_surv_pairs(sex = sex, sf =  sf, sm =  sm, spair = spair, time = initial_time, pairs_f = pairs_f, pairs_m = pairs_m)
  rpair <- compute_recapture(sex = sex, 
                             rpair = rpair,
                             spair = spair, 
                             pairs_f = pairs_f,
                             pairs_m = pairs_m, 
                             time = initial_time,
                             p.m = p.m, 
                             p.f = p.f, 
                             rho = rho,
                             recruit_f = recruit_f,
                             recruit_m = recruit_m)
  
  # Simulate Fates
  for(time in (initial_time+1):k){
    # Compute mating status at t 
    mating_list_t <- compute_mating(sex = sex, time = time, delta = delta, 
                                    recruit = recruit, recruit_f = recruit_f, recruit_m = recruit_m, 
                                    mating_f = mating_f, mating_m = mating_m, mating = mating, 
                                    sf = sf, sm = sm)
    mating_f <- mating_list_t[["mating_f"]] 
    mating_m <- mating_list_t[["mating_m"]]  
    mating <- mating_list_t[["mating"]]
    
    # Will previous pairs from time t-1 reform? 
    repartner <- compute_re_partner(repartner = repartner, 
                                    sex =  sex, 
                                    coef_list = coef_list,
                                    betas = betas, 
                                    pairs = pairs, 
                                    mating = mating, 
                                    time = time)
    
    # Compute partnership probability at t based on surv at t-1
    pairs <- compute_partnerships(sex = sex,
                                  coef_list = coef_list, 
                                  pairs = pairs,
                                  mating = mating,
                                  time = time,
                                  repartner = repartner)
    # Update partner histories going into next survival check (at time k we don't need this)
    if(time < k){
      coef_list <- update_history(coef_list = coef_list,
                                  pairs =  pairs,
                                  time =  time + 1,
                                  sex = sex)
    }

    pairs_ind_list <- propogate_partner_state(pairs = pairs, 
                                              sex =  sex, 
                                              pairs_f = pairs_f, 
                                              pairs_m = pairs_m, 
                                              time = time)
    pairs_f <- pairs_ind_list[["pairs_f"]]
    pairs_m <- pairs_ind_list[["pairs_m"]]
    
    # Compute survival probability at t based on partners at t 
    spair <- compute_survival(spair = spair, 
                              sf = sf, 
                              sm =  sm, 
                              pairs_f = pairs_f,
                              pairs_m = pairs_m, 
                              recruit_f = recruit_f, 
                              recruit_m = recruit_m,
                              time = time,
                              phi.m = phi.m, 
                              phi.f = phi.f,
                              gam = gam)
    
    sind_list_t <- propogate_surv_individual(sf = sf,
                                             sm = sm, 
                                             spair = spair, 
                                             time = time, 
                                             sex = sex)
    sf <- sind_list_t[["sf"]]
    sm <- sind_list_t[["sm"]]
    
    # Survival must happen at least on first capture, debug if not
    if(any(rowSums(sf, na.rm = T) < 1)) browser()
    if(any(rowSums(sm, na.rm = T) < 1)) browser()
    
    # Compute recapture probability at t based on survival at t
    rpair <- compute_recapture(sex = sex,
                               rpair =  rpair,
                               spair = spair, 
                               pairs_f = pairs_f,
                               pairs_m = pairs_m, 
                               time =  time, 
                               p.m = p.m,
                               p.f = p.f, 
                               rho = rho, 
                               recruit_f = recruit_f,
                               recruit_m = recruit_m)
  }
  
  # Grab individual recaptures
  recap_ind_list <- propogate_recap_individual(sex = sex, k = k, rpair = rpair, recruit_f = recruit_f, recruit_m = recruit_m)
  recap_f <- recap_ind_list[["recap_f"]] 
  recap_m <- recap_ind_list[["recap_m"]] 

  # Build partially observed/latent data variables
  
  # Hidden recruitment (Mask unknown states for recruitment variables)
  recruit_list <- compute_hidden_recruitment(recruit_f = recruit_f, 
                                             recruit_m = recruit_m,
                                             recap_f = recap_f,
                                             recap_m = recap_m)
  # Now they have can NAs
  recruit_f <- recruit_list[["recruit_f"]]
  recruit_m <- recruit_list[["recruit_m"]]
  
  # Hidden Survival
  asurv_list <- compute_hidden_survival(pairs_f = pairs_f, 
                                        pairs_m = pairs_m,
                                        rpair = rpair, 
                                        spair = spair,
                                        sf = sf, 
                                        sm = sm, 
                                        k = k, 
                                        sex =  sex)
  af <- asurv_list[["af"]]
  am <- asurv_list[["am"]]
  
  # Hidden Partnerships and mate choice
  apairs_list <- compute_hidden_pairs(pairs_f  = pairs_f,
                                      pairs_m = pairs_m, 
                                      recap_f = recap_f,
                                      recap_m = recap_m,
                                      k = k,
                                      sex = sex,
                                      repartner = repartner,
                                      mating_m = mating_m,
                                      mating_f = mating_f,
                                      show_unmated = show_unmated)
  
  apairs  <- apairs_list[["apairs"]]
  amating_f <- apairs_list[["amating_f"]]
  amating_m <- apairs_list[["amating_m"]]
  arepartner <- apairs_list[["arepartner"]]
  apairs_f <- apairs_list[["apairs_f"]]
  apairs_m <- apairs_list[["apairs_m"]]
  psi <- apairs_list[["psi"]]

  # Return JAGS/NIBMLE (and true) Data
  model_data <- list(
    
    # Known data 
    n = n, # Number of animals sampled
    k = k, # Number of occasions
    nf = length(sex[sex == "F"]), # Number of females
    nm = length(sex[sex == "M"]), # Number of males
    sex = sex, # Sex of sampled individuals
    initial_entry = initial_entry, # When did they enter the population
    recruit_f = recruit_f, 
    recruit_m = recruit_m,
    
    # Latent States (true values - hidden in real data)
    mating_f = mating_f, # Mate status of females at t (+ dummy)
    mating_m = mating_m, # Mate status of males at t (+ dummy)
    pairs_f = pairs_f, # partners of females 
    pairs_m = pairs_m, # partners of males
    pairs = pairs, # pair histories
    known_histories = coef_list[["histories"]], # number of occasions a pair occurred
    sf = sf, # true survival of females
    sm = sm, # true survival of males
    spair = spair, # true survival of partnerships
    repartner = repartner, # whether partner from time t-1 was repicked (female x time)
    psi = psi, # Pairs that may exist (not excluded due to already formed pairs)
    rpair = rpair, # 2D Recapture Matrix
    
    # Observed /Inferred states (Missing Values are possible)
    af = af,  # Female Survival with missing values
    am = am,  # Male Survival with missing values
    apairs  = apairs, # Joint Pairs Matrices (array across time)
    apairs_f = apairs_f,
    apairs_m = apairs_m,
    arepartner = arepartner, # repartner with inferred states 
    amating_f = amating_f, # Mating Status Females at T
    amating_m = amating_m,  # Mating Status Males at T
    recap_f = recap_f, # Observed Recapture of Females
    recap_m = recap_m # Observed Recapture of Males
  )
  
  #Return Model object
  return(model_data)
}

# Wrapper for do-call (helps parallel)
sim_dat <- function(parameter_list){
  model_data <- do.call(simulate_cr_data, parameter_list)
}

# 8. Format Data for Alternate Models -------------------------------------------------------------------------------


format_to_cjs <- function(model_data){
  
  x <- rbind(model_data$recap_f[1:model_data$nf,],model_data$recap_m[1:model_data$nm,])
  a <- rbind(model_data$af[1:model_data$nf,],model_data$am[1:model_data$nm,])
  initial_entry <- c()
  
  # Add female initial capture
  for(i in 1:model_data$nf){
    initial_entry[i] <- min(which(model_data$recruit_f[i,] == 1))
  }
  
  # Add male initial capture
  for(j in (model_data$nf+1):(model_data$nf + model_data$nm)){
    initial_entry[j] <- min(which(model_data$recruit_m[j-model_data$nf,] == 1))
  }
  
  
  # Store results in list
  results <- list(n = model_data$nf + model_data$nm,
                     k = model_data$k,
                     female = c(rep(1, model_data$nf),rep(0, model_data$nm)), 
                     initial_entry = initial_entry,
                     x = x,
                     a = a)
  
  # Return Standard CJS Data
  return(results)
}



format_to_js <- function(model_data){
  
  x <- rbind(model_data$recap_f[1:model_data$nf,],model_data$recap_m[1:model_data$nm,])
  a <- rbind(model_data$af[1:model_data$nf,],model_data$am[1:model_data$nm,])
  recruit <-  rbind(model_data$recruit_f[1:model_data$nf,],model_data$recruit_m[1:model_data$nm,])

  # Store results in list
  model_data <- list(n = model_data$nf + model_data$nm,
                     k = model_data$k,
                     female = c(rep(1, model_data$nf),rep(0, model_data$nm)), 
                     recruit = recruit,
                     x = x,
                     a = a)
  
  # Return Standard CJS Data
  return(model_data)
}

#----------------------------------------------------------------------------------------------------------------------