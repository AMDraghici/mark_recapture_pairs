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

#Inverse Logistic Function
inv.logit <- function(x){
  out <- (1+exp(-x))^(-1)
  return(out)
}


#Softmax Function
softmax <- function(vector){
  out <- exp(vector)/sum(exp(vector))
  return(out)
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
  
  # Chose to mate at time t? Previous state is 0/1 here
  obs <- rbinom(length(previous_state), 1, delta[j]*previous_state)

  return(obs)
}

# 4. Construct Model Matrices------------------------------------------------------------------------------------

# Build out data structures and do easy initialization 

# p-q split between males and females
construct_sexes <- function(n, prop.female = 0.5, random = F){
  # Percent split 
  prop.male <- 1 - prop.female
  #Compute Genders by random sample or fixed count
  if(random == T){
    sex <- sort(ifelse(rbinom(n,1,prop.female)==1,"F","M"))
  } else {
    sex <- sort(c(rep("F",n * prop.female), rep("M", n * prop.male)))
  }
  
  # Return vector of sex designations
  return(sex)
}

construct_init_entry <- function(n, k, random = F, inits = NULL, sort_init = F){
  
  # Select initial entry randomly?
  if(random == T){
    initial_entry = as.integer(sample(c(1:round(k/2)),size=n,replace=T))
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
  recruit_f <- matrix(NA, nrow = n, ncol = k)
  recruit_f[(nf+1):n,1:k] <- 1  # dummy spots for single pairs
  
  #Male's recruits
  recruit_m <- matrix(NA, nrow = n, ncol = k)
  recruit_m[(nm+1):n,1:k] <- 1 # dummy spots for single pairs
  
  #Pair recruits
  recruit <- array(NA, dim = c(n, n, k))
  recruit[(nm+1):n,,1:k] <- 1
  recruit[,(nf+1):n,1:k] <- 1
  
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
  mating_f <- matrix(NA, nrow = n, ncol = k)
  mating_f[(nf+1):n,1:k] <- 1  # dummy spots for single pairs
  
  #Male's mate status
  mating_m <- matrix(NA, nrow = n, ncol = k)
  mating_m[(nm+1):n,1:k] <- 1 # dummy spots for single pairs
  
  #Pair mate status
  mating <- array(NA, dim = c(n, n, k))
  mating[(nm+1):n,,1:k] <- 1
  mating[,(nf+1):n,1:k] <- 1
  
  # Group as list
  mating_list <- list(mating_f = mating_f, 
                       mating_m = mating_m,
                       mating = mating)
  
  # Return mating data
  return(mating_list)
  
}

#Skeleton array for pair assignment at k
construct_pairs <- function(n,k){
  #Female's partners
  pf <- matrix(NA, nrow = n, ncol = k)
  #Male's partners
  pm <- matrix(NA, nrow = n, ncol = k)
  #Pair Identities
  pairs <- array(NA, dim = c(n, n, k))
  # List of objects
  pairs_list <- list(pf = pf, pm = pm, pairs = pairs)
  return(pairs_list)
}

# Construct pair coefficients
construct_coef <- function(sex, k){
  #Intercept
  intercept <- array(1, dim = c(n,n,k))
  # Number of times a pair occurred
  histories <- array(0, dim = c(n,n,k))
  coef_list <- list(intercept =intercept, histories = histories)
  return(coef_list)
}

# Survival matrices for males, females, and 
construct_survival <- function(sex, k){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  #Females
  sf <- matrix(NA, nrow = n, ncol = k)
  sf[(nf+1):n,1:k] <- 0 # dummy spots for single pairs
  #Males
  sm <- matrix(NA, nrow = n, ncol = k)
  sm[(nm+1):n,1:k] <- 0 # dummy spots for single pairs
  #Pairs (and pair-single/single-pair combinations)
  spair <- array(NA, dim = c(length(sex), length(sex), k))
  # Group into list and return data object
  surv_matrices <- list(sf =sf,
                        sm = sm,
                        spair = spair)
  return(surv_matrices)
}

# Recapture Array for pairs (don't need separate ones as recapture is not latent)
construct_recapture <- function(n, k){
  #Pairs (and pair-single/single-pair combinations)
  rpair <- array(NA, dim = c(n, n, k))
  return(rpair)
}

# 5. Initialize CR States---------------------------------------------------------------------------------------

initialize_recruit <- function(n, k, sex, initial_entry, recruit_f, recruit_m, recruit){
 
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
   
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
    recruit[1:nf,1:nm,i] <- recruit_f[1:nf,i] %*% t(recruit_m[1:nm,i])
  }
  
  #Recruit List
  recruit_list <- list(recruit_f = recruit_f, 
                       recruit_m = recruit_m,
                       recruit = recruit)
  
  # Return recruit data 
  return(recruit_list)
}


initialize_mating_choice <- function(n, initial_entry, delta, mating_f, mating_m, mating, sex){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  
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


initialize_partner_status <- function(n, coef_list, betas, pairs, mating, recruit, initial_entry){

  # When did the first animal get spotted
  i <- min(initial_entry) 
  
  psi_t <- compute_partner_probs(n, i, coef_list, betas)
  
  # Flat probability of partnerships
  prob_mate <- psi_t * mating[,,i] * recruit[,,i]
  
  for(j in 1:n){
    # Probability of partnerships forming for j -> 1:n
    probj <- prob_mate[j,]
    
    # Any pairs that have formed cannot form again
    if(j == 2){
      probj <- probj*(1-pairs[1,,i]) #Remove pair j=1 from the prob (colsum only works on matrices)
    } else if(j > 2){
      probj <- probj*(1-colSums(pairs[1:(j-1),,i])) #remove all preformed pairs
    }
    
    # Draw partnership
    pairs[j,,i] <- t(rmultinom(1,1,probj)) # assign a partner 
  }
  
  # Return updated pairs matrix
  return(pairs)
}


initialize_survival_status <- function(n, initial_entry, sex, sf, sm){
  
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  
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
propogate_partner_state <- function(pairs, n, pf, time){
  # Recover index entries
  pf[,time] <- pairs[,,time] %*% (1:n)
  # Return updated information
  return(pf)
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

# Convert marginal to joint survival states
propogate_surv_pairs <- function(sf, sm, spair, time, pf,n){
  
  for(i in 1:n){
    j <- pf[i,time]
    
    state_vector <- c(
                      (1-sf[i,time])*(1-sm[j,time]), #both dead
                      sf[i,time]*(1-sm[j,time]), # female alive
                      (1-sf[i,time])*sm[j,time], # male alive
                      sf[i,time]*sm[j,time] # both alive 
                      )
    state_value <- which(state_vector == 1)
    spair[i,j,time] <- state_value
  }
  
  #Remaining entries are "dead"
  spair[,,time][is.na(spair[,,time])] <- 1
  
  # Return joint survival density
  return(spair)
  
}

compute_partner_lodds <- function(n, time, coef_list, betas){
  
  # Assign memory 
  eta_t <- matrix(0, nrow = n, ncol = n)
  
  beta_vec <- unlist(betas)
  
  # Add all the vectorized coefficients
  for(i in 1:length(beta_vec)){
    eta_t <- eta_t + coef_list[[i]][,,time]*beta_vec[i]
  }
  
  # Return log-odds of partnerships at time t
  return(eta_t)
}


compute_partner_probs <- function(n, time, coef_list, betas){
  
  # Grab log-odds of pairs forming in order
  eta_t <- compute_partner_lodds(n, time, coef_list, betas)
  
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
compute_mating <- function(n, time, delta, recruit, recruit_f, recruit_m,
                           mating_f, mating_m, mating, sf, sm, sex){
  
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
  
  # Mating by females
  for(i in 1:nf){
    
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

# compute partnerships at t
compute_partnerships <- function(n, coef_list, betas, pairs, mating, recruit, time){
  
  # Raw probability of mating (without conditions)
  psi_t <- compute_partner_probs(n, time, coef_list, betas)
  
  # Probability of partnerships with mating and recruitment included (surv is baked into mating)
  prob_mate <- psi_t * mating[,,time] * recruit[,,time]
  
  for(j in 1:n){
    # Probability of partnerships forming for j -> 1:n
    probj <- prob_mate[j,]
    
    # Any pairs that have formed cannot form again
    if(j == 2){
      probj <- probj*(1-pairs[1,,time]) #Remove pair j=1 from the prob (colsum only works on matrices)
    } else if(j > 2){
      probj <- probj*(1-colSums(pairs[1:(j-1),,time])) #remove all preformed pairs
    }
    
    # Draw partnership
    pairs[j,,time] <- t(rmultinom(1,1,probj)) # assign a partner 
  }
  
  # Return updated pairs matrix
  return(pairs)
}
  
# Compute survival state at t
compute_survival <- function(spair, sf, sm, pf, time, phi.m, phi.f, gam, recruit_f, recruit_m){
 
  for(i in 1:nrow(pf)){
    
    # Identify i's partner
    partner_i <- pf[i, time]
    
    # Is initial entry now or later?
    init_f_i <- min(which(recruit_f[i,]==1))
    init_m_i <- min(which(recruit_m[partner_i,]==1))

    # Previous survival states
    surv_i <- sf[i, time-1]*(init_f_i != time) + 1*(init_f_i == time)
    surv_m <- sm[partner_i, time-1]*(init_m_i != time) + 1*(init_m_i == time)
    
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
  
  
  #Remaining entries are "dead"
  spair[,,time][is.na(spair[,,time])] <- 1
  
  # Return joint survival density
  return(spair)
}


# Compute recapture state at t
compute_recapture <- function(rpair, spair, pf, time, p.m, p.f, rho){
  
  # Compute capture status
  for(i in 1:n){
    j <- pf[i, time]
    current_state <- spair[i,j, time]
    rpair[i,j, time] <- recapture(current_state, time, p.m, p.f, rho)
  }
 
  # Turn rest to unobserved 
  rpair[,,time][is.na(rpair[,,time])] <- 1
  
  # Return recapture-pairs matrix
  return(rpair)
}

# 7. Simulate Data ---------------------------------------------------------------------------------------------

#############
# FOR TESTING (delete later)
#############

n <- 1000
k <- 10
prop.female <- 0.5
phi.f <- rep(0.5, k)
phi.m <- phi.f
gam <- rep(0.5, k)
p.f <- phi.f
p.m <- phi.f
rho <- gam
betas <- list(beta0 = 1, beta1 = 0.5)
delta <- rep(0.5, k)
rand_sex <- F
rand_init <- T
init <- NULL

#options(warn = 2)   
#options(warn = 0)   

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
  
  # Generate SKeleton Data Structures
  
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
  pairs_list <- construct_pairs(n, k)
  pf <- pairs_list[["pf"]]
  #pm <- pairs_list[["pm"]]
  pairs <- pairs_list[["pairs"]]
  coef_list <- construct_coef(sex, k)
  survival_list <- construct_survival(sex, k)
  sf <- survival_list[["sf"]]
  sm <- survival_list[["sm"]]
  spair <- survival_list[["spair"]]
  rpair <- construct_recapture(n,k)
  
  
  # Initialize Data 
  initial_time <- min(initial_entry)
  recruit_list <- initialize_recruit(n, k , sex, initial_entry, recruit_f, recruit_m, recruit)
  recruit_f <- recruit_list[["recruit_f"]]
  recruit_m <- recruit_list[["recruit_m"]]
  recruit <- recruit_list[["recruit"]]
  mating_list_init <- initialize_mating_choice(n, initial_entry,delta, mating_f,mating_m, mating, sex)
  mating_f <- mating_list_init[["mating_f"]] 
  mating_m <- mating_list_init[["mating_m"]]  
  mating <- mating_list_init[["mating"]]
  pairs <- initialize_partner_status(n, coef_list, betas, pairs, mating, recruit, initial_entry)
  coef_list <- update_history(coef_list, pairs, initial_time + 1, sex)
  pf <- propogate_partner_state(pairs, n, pf, time = initial_time)
  init_surv_list <-  initialize_survival_status(n, initial_entry, sex, sf, sm)
  sf <- init_surv_list[["sf"]]
  sm <- init_surv_list[["sm"]]
  spair <- propogate_surv_pairs(sf, sm, spair,initial_time, pf, n)
  rpair <- compute_recapture(rpair, spair, pf, initial_time, p.m, p.f, rho)
  
  # Simulate Fates
  for(time in (initial_time+1):k){
    # Compute recruitment at t
    # Ignore for now as I've randomized this for the time being 
    # When coming back to this think about initialization specifically wrt survival 
    
    # Compute mating status at t 
    
    mating_list_t <- compute_mating(n, time, delta, 
                                    recruit, recruit_f, recruit_m, 
                                    mating_f, mating_m, mating, 
                                    sf, sm, sex)
    mating_f <- mating_list_t[["mating_f"]] 
    mating_m <- mating_list_t[["mating_m"]]  
    mating <- mating_list_t[["mating"]]
    
    # Compute partnership probability at t based on surv at t-1
    pairs <- compute_partnerships(n, coef_list, betas, pairs, mating, recruit, time)
    
    # Update partner histories going into next survival check (at time k we don't need this)
    if(time < k){
      coef_list <- update_history(coef_list, pairs, time + 1, sex)
    }
    
    pf <- propogate_partner_state(pairs, n, pf, time)
    
    # Compute survival probability at t based on partners at t 
    spair <- compute_survival(spair, sf, sm, pf, time, phi.m, phi.f, gam, recruit_f, recruit_m)
    sind_list_t <- propogate_surv_individual(sf, sm, spair, time, sex)
    sf <- sind_list_t[["sf"]]
    sm <- sind_list_t[["sm"]]
    
    # Compute recapture probability at t based on survival at t
    rpair <- compute_recapture(rpair, spair, pf, time, p.m, p.f, rho)
    
  }
  
  # Build partially observed/latent data variables
  
  
  # Return JAGS/NIBMLE (and true) Data
  model_data <- list(
    n  = n, # Number of animals sampled
    k = k, # Number of occasions
    nf = length(sex[sex == "F"]), # Number of females
    nm = length(sex[sex == "M"]), # Number of males
    sex = sex, # Sex of sampled individuals
    initial_entry = initial_entry, # When did they enter the population
    mating = mating, # Mate status at time t
    pf = pf, # partners of females (including slots for singles)
    pm = pm, # partners of males (including slots for singles)
    pairs = pairs, # pair histories
    histories = histories, # number of occasions a pair occurred
    sf = sf, # true survival of females
    sm = sm, # true survival of males
    spair = spair, # true survival of partnerships
    rpair = rpair # recapture of pairs 
    )
  
  #Return Model object
  return(model_data)
}



# 7. Format Data for Different Purposes ------------------------------------------------------------------------

format_to_std_cjs <- function(model_data){
  
  
  
  # Store results in list
  model_data <- list()
  
  # Return Standard CJS Data
  return(model_data)
}