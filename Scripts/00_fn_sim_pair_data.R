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

# odds ratio
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


#Softmax Function
# softmax <- function(vector){
#   out <- exp(vector)/sum(exp(vector))
#   return(out)
# }

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
  lub <- sqrt(pmin(or(prob.f,prob.m),1/or(prob.f,prob.m)))
  glb <- -sqrt(pmin(op(prob.f,prob.m),1/op(prob.f,prob.m)))
  #Return values
  out <- list(prob.prod = prob.prod,
              sig.prob.f = sig.prob.f,
              sig.prob.m = sig.prob.m,
              sig.prod = sig.prod,
              cor_lower_bound = glb,
              cor_upper_bound = lub)
  return(out)
}

#Generate Derived Probabilities
compute_jbin_cjs <- function(prob.f,prob.m,corr){
  #Extract Parameters
  parameters <- compute_jbin_param_cjs(prob.f,prob.m)
  cor_upper_bound <- round(parameters[[6]],3)
  cor_lower_bound <- round(parameters[[5]],3)
  corr <- round(corr, 3)
  sig.prob.f <- parameters[[2]]
  sig.prob.m <- parameters[[3]]
  
  # If we are on a boundary condition update correlation structure 
  boundary_f <-  which(prob.f == 1| prob.f == 0)
  boundary_m <- which(prob.m == 1| prob.m == 0)
  boundary <- sort(unique(c(boundary_f, boundary_m)))
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
  prob.distn <- compute_jbin_cjs(phi.m, phi.f, gam)
  phi.mf <- prob.distn[[1]]
  phi.f0 <- prob.distn[[2]]
  phi.m0 <- prob.distn[[3]]
  phi.00 <- prob.distn[[4]]
  
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
  prob.distn <- compute_jbin_cjs(p.m, p.f, rho)
  p.mf <- prob.distn[[1]]
  p.f0 <- prob.distn[[2]]
  p.m0 <- prob.distn[[3]]
  p.00 <- prob.distn[[4]]
  
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

construct_init_entry <- function(n, k, random = F, inits = NULL, sort_init = F){
  
  # Select initial entry randomly?
  if(random == T & is.null(inits)){
    initial_entry <- as.integer(sample(c(1:round(k/2)),size=n,replace=T))
  } else {
    # If not do we have pre-specified values?
    if(!is.null(inits)){
      initial_entry <- inits 
      # If not just do fixed intervals across n and k 
    } else {
      initial_entry <- rep(1:(k-1),ceiling(n/(k-1)))
      initial_entry <- initial_entry[1:n]
    }
  }
  # Arrange by who enters first if desired 
  # No impact on simulation but makes data look cleaner
  if(sort_init == T) initial_entry <- sort(initial_entry)
  
  # Return initial entries
  return(initial_entry)
  
}

# Build recruitment matrix
construct_recruit <- function(n, k, sex){
  
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
  recruit[nm+1,nf+1,1:k] <- 1
  
  # Group as list
  recruit_list <- list(recruit_f = recruit_f, 
                       recruit_m = recruit_m,
                       recruit = recruit)
  
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
  mating[1:(nm+1),nf+1,1:k] <- 1
  mating[(nm+1),1:(nf+1),1:k] <- 1
  
  # Group as list
  mating_list <- list(mating_f = mating_f, 
                      mating_m = mating_m,
                      mating = mating)
  
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
  pairs <- array(NA, dim = c(nf+1, nm+1, k))
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
  surv_matrices <- list(sf =sf,
                        sm = sm,
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
  n <- length(sex)
  
  # initial entry for males and females
  initial_entry_female <- initial_entry[1:nf]
  initial_entry_male <- initial_entry[(nf+1):n]
  
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
                       recruit = recruit)
  
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
  for(i in 1:k){
    mating[1:nf,1:nm,i] <- mating_f[1:nf,i] %*% t(mating_m[1:nm,i])
  }
  
  # Group results and join
  mating_list <- list(mating_f = mating_f,
                      mating_m = mating_m,
                      mating = mating)
  
  return(mating_list)
}


initialize_partner_status <- function(n, coef_list, pairs, mating, recruit, initial_entry, sex){
  
  # When did the first animal get spotted
  i <- min(initial_entry) 
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  
  psi_t <- compute_partner_probs(sex, i, coef_list, c(1,1))
  
  # Bias selection towards mates at the start
  # Small hack but should work out to setting prob of single = 0 unless no mates are available
  psi_t[1:nf,1:nm] <- psi_t[1:nf,1:nm] + 1e6 
  
  # Flat probability of partnerships
  prob_mate <- psi_t * mating[,,i] * recruit[,,i]
  
  for(j in 1:(nf)){
    # Probability of partnerships forming for j -> 1:n
    probj <- prob_mate[j,]
    
    # Any pairs that have formed cannot form again
    if(j == 2){
      probj <- probj*c((1-pairs[1,1:nm,i]),1) #Remove pair j=1 from the prob (colsum only works on matrices)
    } else if(j > 2){
      probj <- probj*c((1-colSums(pairs[1:(j-1),1:nm,i])),1) #remove all pre-formed pairs
    }
    
    # Draw partnership
    pairs[j,,i] <- t(rmultinom(1,1,probj)) # assign a partner 
  }
  
  # Assign unmated males to singles row
  pairs[nf+1,,i] <- c(1-colSums(pairs[1:nf,1:nm,i]),1)
  
  # Return updated pairs matrix
  return(pairs)
}


initialize_survival_status <- function(sex, k, initial_entry, sf, sm){
  
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  n <- length(sex)
  
  # Survival females
  # Assign survival status to initial entry
  for(i in 1:(k-1)){
    
    # Survival females
    enter_at_i_female <- which(initial_entry[1:nf] == i)
    sf[enter_at_i_female,i] <- 1
    
    #survival males 
    enter_at_i_male <- which(initial_entry[(nf+1):n] == i)
    sm[enter_at_i_male,i] <- 1
    
    # Turn the previous entries into NA for likelihood calculation
    if(i > 1){
      sf[enter_at_i_female,1:(i-1)] <- 0
      sm[enter_at_i_male,1:(i-1)] <- 0
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
  pairs_f[,time] <- pairs[1:nf,,time] %*% 1:(nf+1)
  pairs_m[,time] <- t(t(1:(nm+1)) %*% pairs[,1:nm,time])
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
propogate_recap_individual <- function(sex, k, rpair){
  
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
  }
  
  # Only update males (not dummy states)
  for(i in 1:nm){
    for(time in 1:k){
      male_observation <- recap_marginal_male[max(rpair[,i,time])]
      recap_m[i,time] <- male_observation
    }
    
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
  
  # convert to probability using softmax 
  psi_t <- sapply(1:nrow(eta_t),function(x) softmax(eta_t[x,]))
  
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
    prob <- prob * (1-pairs[i,,time-1][nm+1]) # if previously single then set to zero (st prob of new mate isnt biased downward)
    repartner[i,time] <- rbinom(1,1,prob = prob) #  simulate re-pair status
  }
  return(repartner)
}


# compute partnerships at t

compute_partnerships <- function(sex, coef_list, pairs, mating, time, repartner){
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  
  # Raw probability of mating (without conditions)
  psi_t <- compute_partner_probs(sex, time, coef_list, c(1,0))
  
  # Probability of partnerships with mating and recruitment included (surv is baked into mating)
  prob_mate <- psi_t * mating[,,time]
  
  # Make chance of not forming a pair if you're choosing to mate ~ 0
  prob_mate[nf+1,] <- .Machine$double.eps
  prob_mate[,nm+1] <- .Machine$double.eps
  
  # Incorporate repartnership probs
  for(j in 1:nf){
    prob_mate[j, ] <- prob_mate[j,] * (1 - repartner[j,time]) + pairs[j,,time-1] * repartner[j,time]
  }
  
  for(i in 1:nm){
    if(any(prob_mate[,i]==1)){
      impossible_pair_index <- !(1:(nm+1)) %in% which(prob_mate[,i]==1)
      prob_mate[impossible_pair_index,i] <- 0
    }
  }
  
  for(j in 1:nf){
    # Probability of partnerships forming for j -> 1:n
    # If repartner = 1 then a pair reforms 
    probj <- prob_mate[j,]/sum(prob_mate[j,]) 
    
    # Any pairs that have formed cannot form again
    if(j == 2){
      probj <- probj*c((1-pairs[1,1:nm,time]),1) #Remove pair j=1 from the prob (colsum only works on matrices)
    } else if(j > 2){
      probj <- probj*c((1-colSums(pairs[1:(j-1),1:nm,time])),1) #remove all preformed pairs
    }
    
    # Draw partnership
    pairs[j,,time] <- t(rmultinom(1,1,probj)) # assign a partner 
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
    single_check <- (partner_i == nrow(pairs_f) + 1)
    
    # Is initial entry now or later?
    init_f_i <- min(which(recruit_f[i,]==1))
    init_m_i <- min(which(recruit_m[partner_i,]==1))
    
    # Previous survival states
    surv_i <- sf[i, time-1]*(init_f_i != time) + 1*(init_f_i == time)
    surv_m <- (sm[partner_i, time-1]*(init_m_i != time) + 1*(init_m_i == time))*(1-single_check)
    
    # Recover joint survival state for this new couple at t-1
    previous_state_i <- 1 + surv_i + 2*surv_m
    
    # If you've entered now then you're alive t so probability becomes 1 
    phiF <- phi.f
    phiM <- phi.m 
    
    phiF[time-1] <- phi.f[time-1]*(init_f_i != time) + 1*(init_f_i == time)   
    phiM[time-1] <- phi.m[time-1]*(init_m_i != time) + 1*(init_m_i == time)   
    
    # Compute joint survival at time t
    spair[i,partner_i, time] <- survival(previous_state_i,time-1,phiM, phiF, gam)
  }
  
  # Survival for single males
  single_males <- which(pairs_m[,time]==(nrow(pairs_f)+1))
  
  for(j in single_males){
    
    # Is initial entry now or later?
    init_m_j <- min(which(recruit_m[j,]==1))
    
    # Previous survival states
    surv_j <- sm[j, time-1]*(init_m_j != time) + 1*(init_m_j == time)
    
    # Recover joint survival state for this new couple at t-1
    previous_state_j <- 1 + 2*surv_j
    
    # If you've entered now then you're alive t so probability becomes 1 
    phiM <- phi.m
    
    phiM[time-1] <- phi.m[time-1]*(init_m_j != time) + 1*(init_m_j == time)   
    
    # Compute joint survival at time t
    spair[nrow(pairs_f)+1, j, time] <- survival(previous_state_j,time-1,phiM, rep(0,length(phi.f)), gam)
  }
  
  #Remaining entries are "dead"
  spair[,,time][is.na(spair[,,time])] <- 1
  
  # Return joint survival density
  return(spair)
}


# Compute recapture state at t
compute_recapture <- function(sex, rpair, spair, pairs_f, pairs_m, time, p.m, p.f, rho){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  
  # Compute capture status
  for(i in 1:nf){
    j <- pairs_f[i, time]
    current_state <- spair[i,j, time]
    rpair[i,j, time] <- recapture(current_state, time, p.m, p.f, rho)
  }
  
  single_males <- which(spair[nf+1,1:nm,time] != 1)
  
  for(j in single_males){
    current_state <- spair[nf+1,j, time]
    rpair[nf+1, j, time] <- recapture(current_state, time, p.m, p.f, rho)
  }
  
  # Turn rest to unobserved 
  rpair[,,time][is.na(rpair[,,time])] <- 1
  
  # Return recapture-pairs matrix
  return(rpair)
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
  
  # Set dummy states to known ones 
  af[(nf+1),] <- 1
  am[(nm+1),] <- 1
  
  # Store results in a list
  state_list <- list(af = cbind(rep(1,nrow(af)),af),
                     am = cbind(rep(1,nrow(am)),am))
  
  # Return List
  return(state_list)
}



compute_hidden_pairs <- function(pairs_f, pairs_m, rpair, k, sex, repartner){
  
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
  for(i in 1:nf){
    for(time in 1:k){
      
      if(time >=2){
        x_ijt_minus_1 <- x_ijt 
      } else{
        x_ijt_minus_1 <- 4
      }
      
      # Joint Recapture of pair i and pf[i,t] at time t
      x_ijt <- rpair[i, pairs_f[i,time] , time]
      
      # Uncouple Pair Index
      apairs_f[i, time] <- pair_state[x_ijt] * pairs_f[i, time]
      amating_f[i, time] <- pair_state[x_ijt]
      
      # Observed Repartnership
      arepartner[i,time] <- pair_state[x_ijt] * pair_state[x_ijt_minus_1] * repartner[i, time]
      
      if(!is.na(pair_state[x_ijt]))  amating_m[ apairs_f[i, time], time] <- pair_state[x_ijt]
    }
  }
  
  for(j in 1:nm){
    for(time in 1:k){
      
      # Joint Recapture of pair i and pf[i,t] at time t
      x_ijt <- rpair[pairs_m[j, time], j , time]
      
      # Uncouple Pair Index
      apairs_m[j, time] <- pair_state[x_ijt] * pairs_m[j, time]
      amating_m[j, time] <- pair_state[x_ijt]
      if(!is.na(pair_state[x_ijt]))  amating_f[ apairs_m[j, time], time] <- pair_state[x_ijt]
    }
  }
  
  # Add details to joint density
  
  # All pairs (nf+1:n are dummy states)
  for(i in 1:nf){
    
    # If all unknown skip
    if(all(is.na(apairs_f[i,]))) next
    
    time_index <- which(!is.na(apairs_f[i,]))
    
    for(time in time_index){
      # Assign states 
      apairs[i+1,,time+1] <- 0
      apairs[,apairs_f[i, time]+1, time+1] <- 0
      apairs[i+1,apairs_f[i, time]+1, time+1] <- 1
      
    }
  }
  
  # Dummy states always available
  amating_f[(nf+1),1:k] <- 1
  amating_m[(nm+1),1:k] <- 1
  
  apairs[1:(nf+1),1,] <- 0
  apairs[1,1:(nm+1),] <- 0
  apairs[,,1] <- 0 
  arepartner[,1] <- 0 
    
  # Store results in a list
  pairs_list <- list(apairs_m = apairs_m,
                     apairs_f = apairs_f,
                     apairs = apairs,
                     amating_f = amating_f,
                     amating_m = amating_m,
                     arepartner = arepartner)
  
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
                             rand_sex = F,
                             rand_init = T,
                             init = NULL){
  
  # Make sure the number of individuals simulated is even
  if(!n %% 2 == 0){
    n <- n + 1
    if(!is.null(init)){
      init[n] <- 1
    }
  }
  
  # Generate SKeleton Data Structures
  #browser()
  sex <- construct_sexes(n, prop.female, rand_sex)
  initial_entry <- construct_init_entry(n, k, rand_init,init) 
  recruit_list <- construct_recruit(n, k, sex)
  recruit_f <- recruit_list[["recruit_f"]]
  recruit_m  <- recruit_list[["recruit_m"]]
  recruit <- recruit_list[["recruit"]]
  mating_list <- construct_mated(n, k, sex)
  mating_f <- mating_list[["mating_f"]] 
  mating_m <- mating_list[["mating_m"]]  
  mating <- mating_list[["mating"]]
  repartner <- construct_repartner(sex, k)
  pairs_list <- construct_pairs(sex, k)
  pairs_f <- pairs_list[["pairs_f"]]
  pairs_m <- pairs_list[["pairs_m"]]
  pairs <- pairs_list[["pairs"]]
  coef_list <- construct_coef(sex, k)
  survival_list <- construct_survival(sex, k)
  sf <- survival_list[["sf"]]
  sm <- survival_list[["sm"]]
  spair <- survival_list[["spair"]]
  rpair <- construct_recapture(sex,k)
  
  # Initialize Data 
  initial_time <- min(initial_entry)
  recruit_list <- initialize_recruit(sex, k, initial_entry, recruit_f, recruit_m, recruit)
  recruit_f <- recruit_list[["recruit_f"]]
  recruit_m <- recruit_list[["recruit_m"]]
  recruit <- recruit_list[["recruit"]]
  mating_list_init <- initialize_mating_choice(n, k, initial_entry,delta, mating_f, mating_m, mating, sex)
  mating_f <- mating_list_init[["mating_f"]] 
  mating_m <- mating_list_init[["mating_m"]]  
  mating <- mating_list_init[["mating"]]
  pairs <- initialize_partner_status(n, coef_list, pairs, mating, recruit, initial_entry, sex)
  coef_list <- update_history(coef_list, pairs, initial_time + 1, sex)
  pairs_ind_list <- propogate_partner_state(pairs, sex, pairs_f, pairs_m, initial_time)
  pairs_f <- pairs_ind_list[["pairs_f"]]
  pairs_m <- pairs_ind_list[["pairs_m"]]
  init_surv_list <-  initialize_survival_status(sex, k, initial_entry, sf, sm)
  sf <- init_surv_list[["sf"]]
  sm <- init_surv_list[["sm"]]
  spair <- propogate_surv_pairs(sex, sf, sm, spair, initial_time, pairs_f, pairs_m)
  rpair <- compute_recapture(sex, rpair, spair, pairs_f, pairs_m, initial_time, p.m, p.f, rho)
  
  # Simulate Fates
  for(time in (initial_time+1):k){
    
    # Compute mating status at t 
    mating_list_t <- compute_mating(sex, time, delta, 
                                    recruit, recruit_f, recruit_m, 
                                    mating_f, mating_m, mating, 
                                    sf, sm)
    mating_f <- mating_list_t[["mating_f"]] 
    mating_m <- mating_list_t[["mating_m"]]  
    mating <- mating_list_t[["mating"]]
    
    # Will previous pairs from time t-1 reform? 
    repartner <- compute_re_partner(repartner, sex, coef_list, betas, pairs, mating, time)
    
    # Compute partnership probability at t based on surv at t-1
    pairs <- compute_partnerships(sex, coef_list, pairs, mating, time, repartner)
    # Update partner histories going into next survival check (at time k we don't need this)
    if(time < k){
      coef_list <- update_history(coef_list, pairs, time + 1, sex)
    }

    pairs_ind_list <- propogate_partner_state(pairs, sex, pairs_f, pairs_m, time)
    pairs_f <- pairs_ind_list[["pairs_f"]]
    pairs_m <- pairs_ind_list[["pairs_m"]]
    
    # Compute survival probability at t based on partners at t 
    spair <- compute_survival(spair, sf, sm, pairs_f, pairs_m, recruit_f, recruit_m, time, phi.m, phi.f, gam)
    sind_list_t <- propogate_surv_individual(sf, sm, spair, time, sex)
    sf <- sind_list_t[["sf"]]
    sm <- sind_list_t[["sm"]]
    
    # Compute recapture probability at t based on survival at t
    rpair <- compute_recapture(sex, rpair, spair, pairs_f, pairs_m, time, p.m, p.f, rho)
  }
  
  # Grab individual recaptures
  recap_ind_list <- propogate_recap_individual(sex, k, rpair)
  recap_f <- recap_ind_list[["recap_f"]] 
  recap_m <- recap_ind_list[["recap_m"]] 

  # Build partially observed/latent data variables
  
  # Hidden Survival
  asurv_list <- compute_hidden_survival(pairs_f, pairs_m, rpair, spair, sf, sm, k, sex)
  af <- asurv_list[["af"]]
  am <- asurv_list[["am"]]
  
  # Hidden Partnerships and mate choice
  apairs_list <- compute_hidden_pairs(pairs_f, pairs_m, rpair, k, sex, repartner)
  apairs  <- apairs_list[["apairs"]]
  amating_f <- apairs_list[["amating_f"]]
  amating_m <- apairs_list[["amating_m"]]
  arepartner <- apairs_list[["arepartner"]]
  
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
    
    # Observed /Inferred states (Missing Values are possible)
    af = af,  # Female Survival with missing values
    am = am,  # Male Survival with missing values
    apairs  = apairs, # Joint Pairs Matrices (array across time)
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

# 7. Format Data for Different Purposes -------------------------------------------------------------------------------


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
  model_data <- list(n = model_data$nf + model_data$nm,
                     k = model_data$k,
                     female = c(rep(1, model_data$nf),rep(0, model_data$nm)), 
                     initial_entry = initial_entry,
                     x = x,
                     a = a)
  
  # Return Standard CJS Data
  return(model_data)
}

#----------------------------------------------------------------------------------------------------------------------