## Generic Functions -------------------------------------------------------------------------------------

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
  par1  <- sort(par, decreasing = TRUE)
  Lk    <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
}

## S/R Distribution Functions -------------------------------------------------------------------------------------

#Extract joint binomial parameters
compute_jbin_param_cjs <- function(prob.f,prob.m){
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
compute_jbin_cjs <- function(prob.f,prob.m,corr){
  
  parameters      <- compute_jbin_param_cjs(prob.f,prob.m)
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

# Construct Model Matrices------------------------------------------------------------------------------------

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
    sex <- c(rep("F",ceiling(n * prop.female)), rep("M", floor(n * prop.male)))
  }
  
  # Return vector of sex designations
  return(sort(sex))
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
  
  
  # Group as list
  recruit_list <- list(recruit_f = recruit_f, 
                       recruit_m = recruit_m)
  
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
  
  # Group as list
  mating_list <- list(mating_f = mating_f, 
                      mating_m = mating_m)
  
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
  pairs_list <- list(pairs_f = pairs_f, 
                     pairs_m = pairs_m, 
                     pairs   = pairs)
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
construct_histories <- function(sex, k){
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  # Number of times a pair occurred
  histories <- array(0, dim = c(nf+1,nm+1,k))
  return(histories)
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
  sm[nm+1,1:k] <- 1 # dummy spots for single pair
  
  # Group into list and return data object
  surv_matrices <- list(sf  = sf,
                        sm  = sm)
  return(surv_matrices)
}

# Recapture Array for pairs (don't need separate ones as recapture is not latent)
construct_recapture <- function(sex, k){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  
  #Females
  recap_f <- matrix(NA, nrow = nf+1, ncol = k)
  recap_f[nf+1,1:k] <- 0 # dummy spots for single pairs
  
  #Males
  recap_m <- matrix(NA, nrow = nm+1, ncol = k)
  recap_m[nm+1,1:k] <- 0 # dummy spots for single pairs
  
  # Group into list and return data object
  recap_matrices <- list(recap_f    = recap_f,
                         recap_m    = recap_m)
  return(recap_matrices)
  
  return(recap_matrices)
}

# 5. Initialize CR States---------------------------------------------------------------------------------------

initialize_recruit <- function(sex, 
                               k, 
                               initial_entry, 
                               recruit_f, 
                               recruit_m){
  
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
  
  #Recruit List
  recruit_list <- list(recruit_f = recruit_f, 
                       recruit_m = recruit_m)
  
  # Return recruit data 
  return(recruit_list)
}


initialize_mating_choice <- function(n, k, initial_entry, delta, mating_f, mating_m, sex){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  n <- length(sex)
  
  # initial entry for males and females
  initial_entry_female <- initial_entry[1:nf]
  initial_entry_male <- initial_entry[(nf+1):n]
  
  # Choose whether to mate or not
  initial_mating_choice <- rbinom(n, 1, delta[initial_entry])
  
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
  
  
  # Group results and join
  mating_list <- list(mating_f = mating_f,
                      mating_m = mating_m)
  
  return(mating_list)
}


initialize_partner_status <- function(n, pairs, mating_f, mating_m, initial_entry, sex){
  
  # When did the first animal get spotted
  time <- min(initial_entry) 
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  
  prob_mate <- matrix(0, nrow = nf+1, ncol = nm+1)
  # Flat probability of partnerships
  for(i in 1:nf){
    for(j in 1:nm){
      prob_mate[i,j] <- mating_f[i,time] * mating_m[j,time]
    }
  }
  
  females_mating <- sum(mating_f[1:nf,time])
  males_mating <- sum(mating_m[1:nm,time])
  osr <- females_mating/max(1, males_mating) 
  
  if(osr <= 1){ # Pick by female
    for(i in 1:(nf)){
      
      # Probability of partnerships forming for j -> 1:n
      probi <- prob_mate[i,]
      
      # Any pairs that have formed cannot form again
      if(i == 2){
        probi <- probi*c((1-pairs[1,1:nm,time]),1) #Remove pair j=1 from the prob (colsum only works on matrices)
      } else if(i > 2){
        probi <- probi*c((1-colSums(pairs[1:(i-1),1:nm,time])),1) #remove all pre-formed pairs
      }
      
      if(sum(probi[1:nm])==0){ 
        probi[nm+1] <- 1
      } else {
        probi[nm+1] <- 0
      }
      
      # Draw partnership
      pairs[i,,time] <- t(rmultinom(1,1,probi)) # assign a partner 
    }
    
    # Assign unmated males to singles row
    pairs[nf+1,,time] <- c(1-colSums(pairs[1:nf,1:nm,time]),1)
    
  } else { # Pick by male
    
    for(j in 1:(nm)){
      
      # Probability of partnerships forming for j -> 1:n
      probj <- prob_mate[,j]
      
      # Any pairs that have formed cannot form again
      if(j == 2){
        probj <- probj*c((1-pairs[1:nf,1,time]),1) #Remove pair j=1 from the prob (colsum only works on matrices)
      } else if(j > 2){
        probj <- probj*c((1-rowSums(pairs[1:nf, 1:(j-1), time])),1) #remove all pre-formed pairs
      }
      
      if(sum(probj[1:nf])==0){ 
        probj[nf+1] <- 1
      } else {
        probj[nf+1] <- 0
      }
      
      # Draw partnership
      pairs[,j,time] <- t(rmultinom(1,1,probj)) # assign a partner 
    }
    
    # Assign unmated females to singles column
    pairs[,nm+1,time] <- c(1-rowSums(pairs[1:nf,1:nm,time]),1)
  }
  
  
  # Return updated pairs matrix
  return(pairs)
}


initialize_survival_status <- function(sex, k, initial_entry, sf, sm){
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  n  <- length(sex)
  
  # Survival females
  # Assign survival status to initial entry
  for(t in 1:(k-1)){
    
    # Survival females
    enter_at_t_female <- which(initial_entry[1:nf] == t)
    sf[enter_at_t_female,t] <- 1
    
    #survival males 
    enter_at_t_male <- which(initial_entry[(nf+1):n] == t)
    sm[enter_at_t_male,t] <- 1
  }
  
  # Add survival to times prior (alive until recruited + perished)
  for(i in 1:nf){
    init_i <- min(which(sf[i,]==1))
    if(init_i == 1) next
    sf[i,1:(init_i-1)] <- 1
  }
  
  for(j in 1:nm){
    init_j <- min(which(sm[j,]==1))
    if(init_j == 1) next
    sm[j,1:(init_j-1)] <- 1
  }
  
  # Everyone is alive at start of study otherwise not in data
  sf[,1] <- 1
  sm[,1] <- 1
  
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


# Update History
update_history <- function(histories, pairs, time, sex){
  
  # Number of females and males (ignore singles)
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  
  # Add partnerships at time t-1 to current partnership history at t
  histories[1:nf,1:nm,time] <- histories[1:nf,1:nm,time-1] + pairs[1:nf,1:nm,time-1]
  # Return coefficients 
  return(histories)
}

# compute mating decision at t
compute_mating <- function(sex, 
                           time, 
                           delta, 
                           recruit_f,
                           recruit_m,
                           mating_f, 
                           mating_m, 
                           sf,
                           sm){
  
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
    mating_f[i,time] <- rbinom(1, 1, sf[i,time] * delta[time] * recruit_f[i,time])
    
  }
  
  # Mating by Males
  for(j in 1:nm){
    
    # Initial Entry
    init_m_j <- min(which(recruit_m[j,]==1))
    
    # If just entered then initialize function has assigned status 
    if(init_m_j >= time) next
    
    # Otherwise update status 
    mating_m[j,time] <- rbinom(1, 1, sm[j,time] * delta[time]* recruit_m[j,time])
    
  }
  
  # Group results and join
  mating_list <- list(mating_f = mating_f,
                      mating_m = mating_m)
  
  return(mating_list)
}

# Does pair ij from t-1 reform with probability inv.logit(beta0 + beta1*history[i,j])
compute_re_partner <- function(repartner, sex, histories, betas, pairs, mating_f, mating_m, time){
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  
  # Compute probability of reforming pairs
  for(i in 1:nf){
    previous_partner_i <- which(pairs[i,,time-1]==1) # Previous mate 
    history_ij <- histories[i,previous_partner_i,time] # number of times this pair formed
    prob <- inv.logit(sum(c(1.0,history_ij) * unlist(betas))) # probability of reforming
    prob <- prob * mating_f[i,time] * mating_m[previous_partner_i,time] # both have to be willing to mate to repair
    prob <- prob * (1-pairs[i,nm+1,time-1]) # if previously single then set to zero (st prob of new mate isnt biased downward)
    repartner[i,time] <- rbinom(1,1,prob = prob) #  simulate re-pair status
  }
  return(repartner)
}


# compute partnerships at t

compute_partnerships <- function(sex, pairs, mating_f, mating_m, time, repartner){
  
  nf <- length(sex[sex == "F"]) # Number of females
  nm <- length(sex[sex == "M"]) # Number of males
  
  prob_mate <- matrix(0, nrow = nf+1, ncol = nm+1)
  
  male_taken <- rep(0,nm)
  
  for(j in 1:nm){
    male_taken[j] <- sum(pairs[1:nf,j,time-1] * repartner[1:nf, time])
    if(male_taken[j] > 1|male_taken[j] <0) browser()
  }
  
  for(i in 1:nf){
    for(j in 1:nm){
      prob_mate[i,j] <-  mating_f[i,time] * mating_m[j,time] *(1-pairs[i,j,time-1])* (1 - repartner[i,time]) *(1-male_taken[j])  + 
        pairs[i,j,time-1] * male_taken[j] * repartner[i,time]
    }
  }
  
  
  females_mating <- sum(mating_f[1:nf,time])
  males_mating <- sum(mating_m[1:nm,time])
  osr <- females_mating/max(1, males_mating) 
  
  if(osr <= 1){ # Pick by female
    
    for(i in 1:nf){
      # Probability of partnerships forming for j -> 1:n
      # If repartner = 1 then a pair reforms 
      probi <- prob_mate[i,]
      
      # Any pairs that have formed cannot form again
      if(i == 2){
        probi <- probi*c((1-pairs[1,1:nm,time]),1) #Remove pair j=1 from the prob (colsum only works on matrices)
      } else if(i > 2){
        probi <- probi*c((1-colSums(pairs[1:(i-1),1:nm,time])),1) #remove all pre-formed pairs
      }
      
      if(sum(probi[1:nm])==0){ 
        probi[nm+1] <- 1
      } else {
        probi[nm+1] <- 0
      }
      
      if(any(probi < 0)) browser()
      
      # If repartnering just set to previous partner
      if(repartner[i,time] ==1){
        pairs[i,,time] <- 0
        pairs[i,which(pairs[i,,time-1]==1),time] <- 1
      } else {
        # Draw partnership
        pairs[i,,time] <- t(rmultinom(1,1,probi)) # assign a partner 
      }
      
    }
    
    # Assign unmated males to singles row
    pairs[nf+1,,time] <- c(1-colSums(pairs[1:nf,1:nm,time]),1)
    
  } else {
    for(j in 1:nm){
      # Probability of partnerships forming for j -> 1:n
      # If repartner = 1 then a pair reforms 
      probj <- prob_mate[,j]
      
      # Any pairs that have formed cannot form again
      if(j == 2){
        probj <- probj*c((1-pairs[1:nf,1,time]),1) #Remove pair j=1 from the prob (colsum only works on matrices)
      } else if(j > 2){
        probj <- probj*c((1-rowSums(pairs[1:nf,1:(j-1),time])),1) #remove all pre-formed pairs
      }
      
      if(sum(probj[1:nf])==0){ 
        probj[nf+1] <- 1
      } else {
        probj[nf+1] <- 0
      }
      
      if(any(probj < 0)) browser()
      
      previous_partner <- which(pairs[,j,time-1] == 1)
      
      if(previous_partner != (nf+1)){
        if(repartner[previous_partner,time] == 1){
          pairs[,j,time] <- 0
          pairs[previous_partner,j,time] <- 1
        } else {
          # Draw partnership
          pairs[,j,time] <- t(rmultinom(1,1,probj)) # assign a partner 
        }
      } else {
        # Draw partnership
        pairs[,j,time] <- t(rmultinom(1,1,probj)) # assign a partner 
      }
    }
    
    # Assign unmated males to singles row
    pairs[,nm+1,time] <- c(1-rowSums(pairs[1:nf,1:nm,time]),1)
  }
  
  # Return updated pairs matrix
  return(pairs)
}

# Compute survival state at t
compute_survival <- function(sf, sm, pairs_f, pairs_m, recruit_f, recruit_m, time, phi.m, phi.f, gam){
  
  # Grab number of females/males
  nf <- nrow(pairs_f)
  nm <- nrow(pairs_m)
  
  # Joint density 
  joint_surv_pmf <- compute_jbin_cjs(phi.f[time],phi.m[time],gam[time])
  
  # Sample males first
  for(j in 1:nm){
    # Always alive at time 1 or during first entry
    first_j <- min(which(recruit_m[j,]==1))
    if(time == 1|first_j == time){
      sm[j,time] <- 1
      next
    }
    sm[j,time] <- rbinom(1,1, prob = phi.m[time] * sm[j,time-1] * recruit_m[j,time-1] + (1-recruit_m[j,time-1]))
  }
  
  # Sample conditional females next 
  for(i in 1:nf){
    # Always alive at time 1
    # First entry needs to survive 
    first_i <- min(which(recruit_f[i,]==1))
    if(time == 1|first_i == time){
      sf[i,time] <- 1
      next
    }
    # Find if i has a partner
    j <- pairs_f[i, time-1]
    single_check <- 1*(j == (nm + 1))
    
    if(single_check != 1){
      male_previous_alive <- recruit_m[j,time-1] * sm[j,time-1]
    } else {
      male_previous_alive <- 0
    }
    
    
    # Conditional Probability of survival given pair status
    prob_cond_f <- (single_check + (1-single_check) * (1-male_previous_alive)) * phi.f[time] + 
      (1-single_check) * male_previous_alive * (sm[j,time] * joint_surv_pmf$prob.mf/phi.m[time]  + 
                                                  (1-sm[j,time]) * joint_surv_pmf$prob.f0/(1-phi.m[time]))
    
    
    sf[i,time] <- rbinom(1,1, prob = prob_cond_f * sf[i,time-1] * recruit_f[i,time-1] + (1-recruit_f[i,time-1]))
  }
  
  # Return joint survival density
  return(list(sm=sm,
              sf=sf))
}

compute_hidden_recruitment <- function(recruit_f, recruit_m, recap_f, recap_m){
  
  compute_first_cap <- function(x){
    if(sum(x) == 0) return(0)
    return(min(which(x == 1)))
  }
  
  # First captures
  first_female <- sapply(1:(nrow(recap_f)-1), function(y) compute_first_cap(recap_f[y,]))
  first_male   <- sapply(1:(nrow(recap_m)-1), function(y) compute_first_cap(recap_m[y,]))
  
  for(i in 1:length(first_female)){
    
    if(first_female[i] == 0){
      recruit_f[i,1:(ncol(recruit_f)-1)] <- NA
      next
    }
    
    if(first_female[i] == 1) next
    
    recruit_f[i, 1:(first_female[i]-1)] <- NA
    
  }
  
  
  for(i in 1:length(first_male)){
    if(first_male[i] == 0){
      recruit_m[i,1:(ncol(recruit_m)-1)] <- NA
      next
    }
    if(first_male[i] == 1) next
    
    recruit_m[i, 1:(first_male[i]-1)] <- NA
    
  }
  
  # Must have been recruited by the end of the study otherwise code will break
  # Some animals are not observed so this coding is necessary 
  recruit_f[,ncol(recruit_f)] <- 1
  recruit_m[,ncol(recruit_m)] <- 1
  
  return(list(recruit_f = recruit_f,
              recruit_m = recruit_m))
  
  
}


compute_hidden_survival <- function(recap_f,
                                    recap_m, 
                                    first_capture_f,
                                    first_capture_m,
                                    sf,
                                    sm, 
                                    k, 
                                    sex){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"])  
  n  <- length(sex)
  
  # Produce inferred survival states
  af <- matrix(NA, nrow = nrow(sf), ncol = ncol(sf))
  am <- matrix(NA, nrow = nrow(sm), ncol = ncol(sm))
  
  # If seen set survival to actual status
  # Otherwise unknown 
  
  # Females
  for(i in 1:nf){
    # Add observed status iteratively--------------------------------------
    for(t in first_capture_f[i]:k){
      if(recap_f[i,t]==1){
        af[i,t] <- sf[i,t]
      } 
    }
    
    # Go back and add inferred states based on time-------------------------
    # Last time seen alive
    last_alive <- max(which(af[i,]==1))
    
    # If the last time they were seen alive was time 1 then next 
    if(last_alive == 1) next
    
    # Change all other states prior to last alive to alive
    if(first_capture_f[i] > 1)  af[i,  1:(first_capture_f[i]-1)] <- 0   
    af[i,  first_capture_f[i]:last_alive] <- 1    
    
    rm(last_alive)
    
  }
  
  rm(i)
  
  # Males
  for(j in 1:nm){
    
    # Add observed status iteratively--------------------------------------
    for(t in first_capture_m[j]:k){
      if(recap_m[j,t]==1){
        am[j,t] <- sm[j,t]
      } 
    }
    
    # Last time seen alive
    last_alive <- max(which(am[j,]==1))
    
    # If the last time they were seen alive was time 1 then skip 
    if(last_alive == 1) next 
    
    # Change all other states prior to last alive to alive
    if(first_capture_m[j] > 1)  am[j,  1:(first_capture_m[j]-1)] <- 0   
    am[j, first_capture_m[j]:last_alive] <- 1    
    
    rm(last_alive)
  }

  # Set dummy states to known ones (used for indexing in nimble model)
  af[(nf+1),] <- 0
  am[(nm+1),] <- 0
  
  # Store results in a list
  # add dummy row 1 for the repartner mechanism (in model index t is t+1)
  state_list <- list(af = af,  
                     am = am) 
  
  # Return List
  return(state_list)
}


filter_impossible_matchs <- function(psi,
                                     apairs_f,
                                     pairs_f,
                                     k,
                                     recap_f,
                                     recap_m,
                                     nf,
                                     nm){
  for(t in 1:k){
    for(i in 1:nf){
      for(j in 1:nm){
        if(recap_f[i,t]==1 & recap_m[j,t]==1 & is.na(apairs_f[i,t])){
          psi[i,j,t] <- 1 * (pairs_f[i,t] == j)
        }  
      }
    }
  }
  
  
  return(psi)
  
}

compute_hidden_pairs <- function(pairs_f,
                                 pairs_m, 
                                 recap_f,
                                 recap_m, 
                                 recruit_f,
                                 recruit_m,
                                 first_capture_f,
                                 first_capture_m,
                                 k,
                                 sex, 
                                 repartner, 
                                 mating_f,
                                 mating_m, 
                                 show_unmated = F){
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"])  
  n <- length(sex)
  
  # Produce inferred survival states
  apairs_f   <- matrix(NA, nrow = nrow(pairs_f), ncol = ncol(pairs_f))
  apairs_m   <- matrix(NA, nrow = nrow(pairs_m), ncol = ncol(pairs_m))
  apairs     <- array(NA, dim =c(nf+1,nm+1,k+1))
  amating_f  <- matrix(NA, nrow = nf+1, ncol = k)
  amating_m  <- matrix(NA, nrow = nm+1, ncol = k) 
  arepartner <- matrix(NA, nrow = nf, ncol = k)
  
  # Add observations directly from joint recapture matrix
  
  # Females -----------------------------------------------------------------------------------
  
  for(i in 1:nf){
    
    for(t in 1:k){
      
      if(t < first_capture_f[i]){
        amating_f[i,t] <- 0
        apairs_f[i,t] <- nm+1
        arepartner[i,t] <- 0
      }
      
      # If female in unobserved, move to next time step
      if(recap_f[i,t] == 0) next 
      
      # Partner Index
      j <- pairs_f[i, t]
      single_check <- (j == (nm+1))
      
      # Scenario 1, we know if you're unmated when caught
      if(show_unmated){
        
        amating_f[i,t] <- mating_f[i,t]
        
        # Female Single and Observed
        if(single_check|amating_f[i,t] == 0){ 
          if(j != (nm+1)) browser()
          arepartner[i,t] <- 0
          apairs_f[i,t] <- j
          # Female with partner
        } else {
          # Check if male was observed
          if(recap_m[j,t]==1){
            if(j == (nm+1)) browser()
            apairs_f[i,t] <- j
            
            if(t == 1) next # no repartnership at time 1
            
            # If both seen last time, 
            if(recap_m[j,t-1] == 1 & recap_f[i,t-1] == 1){
              pairs_f_na <- is.na(apairs_f[i,t-1])
              pairs_m_na <- (recap_f[pairs_m[j,t-1],t-1] == 0)  
              if(pairs_f_na & pairs_m_na) arepartner[i,t] <- NA
              if(pairs_f_na & !pairs_m_na) arepartner[i,t] <- 0
              if(!pairs_f_na){ # needs to be nested o/w NAs will cause issues
                if(apairs_f[i,t-1] == j) arepartner[i,t] <- 1
              }
            } 
          } 
        }
        
        # Scenario 2, we dont know if you're unmated when caught  
      } else {
        
        # If the partner is caught then 
        if(recap_m[j,t] == 1){
          
          apairs_f[i,t] <- j
          amating_f[i,t] <- mating_f[i,t]
          
          if(t == 1) next # no repartnership at time 1
          
          # If both seen last time, 
          if(recap_m[j,t-1] == 1 & recap_f[i,t-1] == 1){
            pairs_f_na <- is.na(apairs_f[i,t-1])
            pairs_m_na <- (recap_f[pairs_m[j,t-1],t-1] == 0)  
            if(pairs_f_na & pairs_m_na) arepartner[i,t] <- NA
            if(pairs_f_na & !pairs_m_na) arepartner[i,t] <- 0
            if(!pairs_f_na){
              if(apairs_f[i,t-1] == j) arepartner[i,t] <- 1
            }
          } 
        } 
      }
    }
  }
  
  # Males -----------------------------------------------------------------------------------------------
  
  for(j in 1:nm){
    for(t in 1:k){
      
      if(t < first_capture_m[j]){
        amating_m[j,t] <- 0
        apairs_m[j,t] <- nf+1
      }
      
      # If male in unobserved, move to next time step
      if(recap_m[j,t] == 0) next 
      
      # Partner Index
      i <- pairs_m[j, t]
      single_check <- (i == (nf+1))
      
      # Scenario 1, we know if you're unmated when caught
      if(show_unmated){
        
        amating_m[j,t] <- mating_m[j,t]
        
        # Male Single and Observed
        if(single_check|amating_m[j,t] == 0){ 
          if(i != (nf+1)) browser()
          apairs_m[j,t] <- i
          # Male with partner
        } else {
          # Check if female was observed
          if(recap_f[i,t]==1){
            apairs_m[j,t] <- i
          } 
        }
        # Scenario 2, we dont know if you're unmated when caught  
      } else {
        # If the partner is caught then 
        if(recap_f[i,t] == 1){
          apairs_m[j,t] <- i
          amating_m[j,t] <- mating_m[j,t]
        } 
      }
    }
  }
  
  # #
    # apairs_f <- pairs_f
    # apairs_m <- pairs_m
    # arepartner <- repartner
    # amating_f <- mating_f
    # amating_m <- mating_m
  
  
  # Add details to 2d apairs matrix ------------------------------------------------------------------------------------------
  
  # Set up female exclusions and male partnerships
  for(i in 1:nf){
    
    # If all unknown skip
    if(all(is.na(apairs_f[i,]))) next
    
    # We can only update known information
    time_index <- which(!is.na(apairs_f[i,]))
    
    # Add known pairings and designate known non-pairs
    for(time in time_index){
      # Assign states 
      apairs[i+1,,time+1] <- 0 # make all pairs unmated
      if(apairs_f[i,time] == (nm + 1)){next} # if single next 
      apairs[,apairs_f[i, time]+1, time+1] <- 0 # make all male partner pairs unavailable
      apairs[i+1,apairs_f[i, time]+1, time+1] <- 1 # assign pairing to the correct male/female combo at time t
      
    }
  }
  
  # Set up male exclusions and female partnerships
  for(j in 1:nm){
    
    # If all unknown skip
    if(all(is.na(apairs_m[j,]))) next
    
    # We can only update known information
    time_index <- which(!is.na(apairs_m[j,]))
    
    # Add known pairings and designate known non-pairs
    for(time in time_index){
      # Assign states 
      apairs[,j+1,time+1] <- 0 # make all pairs unmated
      if(apairs_m[j,time] == (nf + 1)){next} # if single next 
      apairs[apairs_m[j, time]+1, , time+1] <- 0 # make all female partner pairs unavailable
      apairs[apairs_m[j, time]+1,j+1, time+1] <- 1 # assign pairing to the correct male/female combo at time t
      
    }
  }
  
  
  # Build arepartner matrix -----------------------------------------------------------------------------------------
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
        previous_partners <- apairs_f[(1:nf)[-i],time-1]  # Grab data 
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
        current_partners <- apairs_f[(1:nf)[-i],time]  # Grab data 
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
  
  # Build partnership grid
  psi <- apairs
  psi[is.na(psi)] <- 1
  psi <- psi[1:nf+1,1:nm+1,]
  psi <- psi[,,1:k+1] 
  psi_array <- array(NA,dim = c(nf+1,nm+1,k))
  
  for(t in 1:k){
    psi_array[,,t] <- rbind(cbind(psi[,,t],rep(1,(nf))),rep(1, nm+1))
  }
  
  psi <- psi_array
  
  # Special Arepartner Logic ------------------------------------------------------------------------------------------
  # Dummy states 
  amating_f[(nf+1),1:k] <- 0
  amating_m[(nm+1),1:k] <- 0
  apairs[1:(nf+1),1,] <- 0
  apairs[1,1:(nm+1),] <- 0
  apairs[,,1] <- 0 
  arepartner[,1] <- 0 
  
  if(show_unmated){
    psi <- filter_impossible_matchs(psi      = psi,
                                    apairs_f = apairs_f,
                                    pairs_f  = pairs_f,
                                    k        = k,
                                    recap_f  = recap_f,
                                    recap_m  = recap_m,
                                    nf       = nf,
                                    nm       = nm)
  }
  # browser()
  
  # Filter Down to Pairs that were observed 
  possible_partners <- lapply(1:nf, function(i) unique(c(nm+1, apairs_f[i, 1:k][!is.na(apairs_f[i,1:k])])))

  for(i in 1:nf){
    if(length(possible_partners[[i]]) == 0){
      psi[i,1:nm, 1:k] <- 0
    } else {
      psi[i, -possible_partners[[i]], 1:k] <- 0
    }
  }
  
  # Clear out single prob 
  for(i in 1:nf){
    for(t in 1:k){
      if(!is.na(apairs_f[i,t])){
        j <- apairs_f[i,t]
        psi[i,1:(nm+1),t] <- 0
        if(j != (nm+1)) psi[1:(nf+1),j,t] <- 0
        psi[i,j,t] <- 1
      }
    }
  }
  
  last_partner_matrix <- matrix(NA, nrow = nf, ncol = k)
  next_partner_matrix <- matrix(NA, nrow = nf, ncol = k)
  # 
  # # Apply Interval Filter (can only be with next/previous/single)-----------------------------------------------------------------------------------------------
  for(i in 1:nf){

    partners_i <- unique(apairs_f[i,][!is.na(apairs_f[i,])])
    first_partner <- partners_i[1]

    # If no pairs observed make no possible combos
    if(length(partners_i) == 0) next

    for(t in 1:k){
      # Who was individual i's last partner?
      last_partner <- unique(apairs_f[i,1:t][!is.na(apairs_f[i,1:t])])
      if(length(last_partner)==0){
        last_partner <- nm+1
      }
      last_partner <- last_partner[length(last_partner)]

      # Who was individual i's next partner?
      # If t == k and unknown then just grab last partner
      if(t == k){
        next_partner <-  unique(apairs_f[i,k][!is.na(apairs_f[i,k])])
        if(length(next_partner) == 0){
          next_partner <- nm+1
        }
      } else {
        next_partner <- unique(apairs_f[i,(t+1):k][!is.na(apairs_f[i,(t+1):k])])
        if(length(next_partner)==0){
          next_partner <- nm+1
        }
        next_partner <- next_partner[1]
      }

      # If observed unmated set to unmated
      if(!is.na(amating_f[i,t])){
        if(amating_f[i,t] == 0){
          apairs_f[i,t] <- (nm+1)
          psi[i,-(nm+1),t] <- 0
          next
        }
      }

      last_partner_matrix[i,t] <- last_partner
      next_partner_matrix[i,t] <- next_partner

      if(is.na(apairs_f[i,t])){

        if(next_partner == last_partner & next_partner != nm+1){
          potential_partners_i <- c(next_partner)
          psi[i,-potential_partners_i,t] <- 0
          apairs_f[i,t] <- next_partner
        } else{
          potential_partners_i <- c(nm+1, last_partner, next_partner)
          psi[i,-potential_partners_i,t] <- 0
        }

      }
    }
  }

  # -----------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # Store results in a list
  pairs_list <- list(apairs_m            = apairs_m,
                     apairs_f            = apairs_f,
                     apairs              = apairs,
                     amating_f           = amating_f,
                     amating_m           = amating_m,
                     arepartner          = arepartner,
                     psi                 = psi,
                     last_partner_matrix = last_partner_matrix,
                     next_partner_matrix = next_partner_matrix)
  
  # Return List
  return(pairs_list)
}

simulate_recapture <- function(recap_f, 
                               recap_m, 
                               first_capture_f, 
                               first_capture_m, 
                               sf, 
                               sm, 
                               pairs_f, 
                               pairs_m,
                               k,
                               sex,
                               p.f,
                               p.m, 
                               rho){
  
  
  # Number of females and males 
  nf <- length(sex[sex == "F"]) 
  nm <- length(sex[sex == "M"]) 
  
  # Generate Male History 
  for(j in 1:nm){
    if(first_capture_m[j] > 1) recap_m[j, 1:(first_capture_m[j]-1)] <- 0
    recap_m[j, first_capture_m[j]] <- 1
    for(t in (first_capture_m[j]+1):k){
      recap_m[j,t] <- rbinom(1,1,prob = p.m[t] * sm[j,t])
    }
  }
  
  # Generate female History
  for(i in 1:nf){
    if(first_capture_f[i] > 1) recap_f[i, 1:(first_capture_f[i]-1)] <- 0
    recap_f[i, first_capture_f[i]] <- 1
    for(t in (first_capture_f[i]+1):k){
      # Joint density 
      joint_recap_pmf <- compute_jbin_cjs(p.f[t],p.m[t],rho[t])
      
      # Find partner info
      j <- pairs_f[i, t]
      single_check <- 1*(j == (nm + 1))
      recap_prob_cond_f <- single_check * p.f[t] + 
                           (1-single_check) * (recap_m[j,t] * joint_recap_pmf$prob.mf/p.m[t] + 
                                              (1-recap_m[j,t]) * joint_recap_pmf$prob.f0/(1-p.m[t]))
      
      
      recap_f[i,t] <- rbinom(1,1, prob = recap_prob_cond_f * sf[i,t])
    }
  }
  
  
  return(list(recap_f = recap_f,
              recap_m = recap_m))
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
                             show_unmated,
                             data_aug){
  
  # Make sure the number of individuals simulated is even
  if(!n %% 2 == 0){
    n <- n + 1
    if(!is.null(init)){
      init[n] <- sample(k-1,1)
    }
  }
  
  # Generate SKeleton Data Structures
  sex             <- construct_sexes(n = n, prop.female = prop.female)
  initial_entry   <- construct_init_entry(n = n, k = k, random = rand_init, init = init) 
  first_capture_f <- initial_entry[sex == "F"]
  first_capture_m <- initial_entry[sex == "M"]
  recruit_list    <- construct_recruit(n = n, k = k, sex = sex)
  recruit_f       <- recruit_list[["recruit_f"]]
  recruit_m       <- recruit_list[["recruit_m"]]
  mating_list     <- construct_mated(n=n, k=k, sex = sex)
  mating_f        <- mating_list[["mating_f"]] 
  mating_m        <- mating_list[["mating_m"]]  
  repartner       <- construct_repartner(sex = sex, k = k)
  pairs_list      <- construct_pairs(sex = sex, k = k)
  pairs_f         <- pairs_list[["pairs_f"]]
  pairs_m         <- pairs_list[["pairs_m"]]
  pairs           <- pairs_list[["pairs"]]
  histories       <- construct_histories(sex = sex, k = k)
  survival_list   <- construct_survival(sex = sex, k = k)
  sf              <- survival_list[["sf"]]
  sm              <- survival_list[["sm"]]
  recap_list      <- construct_recapture(sex = sex, k = k)
  recap_f         <- recap_list[["recap_f"]]
  recap_m         <- recap_list[["recap_m"]]
  
  # Initialize Data 
  initial_time     <- min(initial_entry)
  recruit_list     <- initialize_recruit(sex = sex, k = k, initial_entry = initial_entry, recruit_f = recruit_f, recruit_m = recruit_m)
  recruit_f        <- recruit_list[["recruit_f"]]
  recruit_m        <- recruit_list[["recruit_m"]]
  mating_list_init <- initialize_mating_choice(n = n, k = k, initial_entry = initial_entry,delta = delta, mating_f = mating_f, mating_m = mating_m, sex = sex)
  mating_f         <- mating_list_init[["mating_f"]] 
  mating_m         <- mating_list_init[["mating_m"]]  
  pairs            <- initialize_partner_status(n = n, pairs = pairs, mating_f = mating_f, mating_m = mating_m, initial_entry = initial_entry, sex = sex)
  histories        <- update_history(histories = histories, 
                                     pairs = pairs, time = initial_time + 1, sex = sex)
  pairs_ind_list <- propogate_partner_state(pairs = pairs, sex =  sex, pairs_f = pairs_f, pairs_m = pairs_m, time = initial_time)
  pairs_f        <- pairs_ind_list[["pairs_f"]]
  pairs_m        <- pairs_ind_list[["pairs_m"]]
  init_surv_list <-  initialize_survival_status(sex = sex, k =  k, initial_entry = initial_entry, sf = sf, sm = sm)
  sf             <- init_surv_list[["sf"]]
  sm             <- init_surv_list[["sm"]]
  
  # Survival must happen at least on first entry, debug if not
  if(any(rowSums(sf, na.rm = T) < 1)) browser()
  if(any(rowSums(sm, na.rm = T) < 1)) browser()
  
  # Simulate Fates
  for(time in (initial_time+1):k){
    
    # Compute survival probability at t based on partners at t 
    surv_list <- compute_survival(sf        = sf, 
                                  sm        = sm, 
                                  pairs_f   = pairs_f,
                                  pairs_m   = pairs_m, 
                                  recruit_f = recruit_f, 
                                  recruit_m = recruit_m,
                                  time      = time,
                                  phi.m     = phi.m, 
                                  phi.f     = phi.f,
                                  gam       = gam)
    
    sf <- surv_list[["sf"]]
    sm <- surv_list[["sm"]]
    
    # Survival must happen at least on first capture, debug if not
    if(any(rowSums(sf, na.rm = T) < 1)) browser()
    if(any(rowSums(sm, na.rm = T) < 1)) browser()
    
    # Compute mating status at t 
    mating_list_t <- compute_mating(sex       = sex, 
                                    time      = time, 
                                    delta     = delta, 
                                    recruit_f = recruit_f,
                                    recruit_m = recruit_m, 
                                    mating_f  = mating_f, 
                                    mating_m  = mating_m, 
                                    sf        = sf, 
                                    sm        = sm)
    
    mating_f <- mating_list_t[["mating_f"]] 
    mating_m <- mating_list_t[["mating_m"]]  
    
    # Will previous pairs from time t-1 reform? 
    repartner <- compute_re_partner(repartner = repartner, 
                                    sex       =  sex, 
                                    histories = histories,
                                    betas     = betas, 
                                    pairs     = pairs, 
                                    mating_f  = mating_f, 
                                    mating_m  = mating_m,
                                    time      = time)
    
    # Compute partnership probability at t based on surv at t-1
    pairs <- compute_partnerships(sex       = sex,
                                  pairs     = pairs,
                                  mating_f  = mating_f,
                                  mating_m  = mating_m,
                                  time      = time,
                                  repartner = repartner)
    # Update partner histories going into next survival check (at time k we don't need this)
    if(time < k){
      histories <- update_history(histories = histories,
                                  pairs =  pairs,
                                  time =  time + 1,
                                  sex = sex)
    }
    
    pairs_ind_list <- propogate_partner_state(pairs   = pairs, 
                                              sex     =  sex, 
                                              pairs_f = pairs_f, 
                                              pairs_m = pairs_m, 
                                              time    = time)
    pairs_f <- pairs_ind_list[["pairs_f"]]
    pairs_m <- pairs_ind_list[["pairs_m"]]
    
    
  }
  
  # Update animals without observed histories 
  recap_ind_list <- simulate_recapture(recap_f         = recap_f, 
                                       recap_m         = recap_m, 
                                       first_capture_f = first_capture_f, 
                                       first_capture_m = first_capture_m, 
                                       sf              = sf, 
                                       sm              = sm, 
                                       pairs_f         = pairs_f, 
                                       pairs_m         = pairs_m,
                                       k               = k,
                                       sex             = sex,
                                       p.f             = p.f,
                                       p.m             = p.m, 
                                       rho             = rho)
  
  recap_f <- recap_ind_list[["recap_f"]] 
  recap_m <- recap_ind_list[["recap_m"]] 
  # browser()
  # Hidden Survival
  asurv_list <- compute_hidden_survival(recap_f         = recap_f,
                                        recap_m         = recap_m,
                                        first_capture_f = first_capture_f,
                                        first_capture_m = first_capture_m,
                                        sf              = sf, 
                                        sm              = sm, 
                                        k               = k, 
                                        sex             = sex)
  af <- asurv_list[["af"]]
  am <- asurv_list[["am"]]
  
  # Hidden Partnerships and mate choice
  apairs_list <- compute_hidden_pairs(pairs_f            = pairs_f,
                                      pairs_m            = pairs_m, 
                                      recap_f            = recap_f,
                                      recap_m            = recap_m,
                                      recruit_f          = recruit_f,
                                      recruit_m          = recruit_m,
                                      first_capture_f    = first_capture_f,
                                      first_capture_m    = first_capture_m,
                                      k                  = k,
                                      sex                = sex,
                                      repartner          = repartner,
                                      mating_m           = mating_m,
                                      mating_f           = mating_f,
                                      show_unmated       = show_unmated)
  
  apairs              <- apairs_list[["apairs"]]
  amating_f           <- apairs_list[["amating_f"]]
  amating_m           <- apairs_list[["amating_m"]]
  arepartner          <- apairs_list[["arepartner"]]
  apairs_f            <- apairs_list[["apairs_f"]]
  apairs_m            <- apairs_list[["apairs_m"]]
  psi                 <- apairs_list[["psi"]]
  next_partner_matrix <- apairs_list[["next_partner_matrix"]]
  last_partner_matrix <- apairs_list[["last_partner_matrix"]]
  # nf and nm
  sex_counts <- table(sex)
  nf         <- sex_counts[1]
  nm         <- sex_counts[2]
  
  # FOR TESTING
  psi_check  <- psi 
  # # browser()
  # apairs_f <- pairs_f
  # am <- sm
  
  
  az <- matrix(NA, nrow = nf, ncol = k)
  ar <- matrix(NA, nrow = nf, ncol = k)
  
  # Add AR/AZ known
  for(i in 1:nf){
    for(t in first_capture_f[i]:(k-1)){
      if(!is.na(apairs_f[i,t])){
        if(apairs_f[i,t] == (nm+1)){
          az[i,t] <- 1
        } else  {
          az[i,t] <- 2 + am[apairs_f[i,t],t+1]
        }
      }
    }
  }
  
  for(i in 1:nf){
    for(t in (first_capture_f[i]):k){
      if(!is.na(apairs_f[i,t])){
        if(apairs_f[i,t] == (nm+1)){
          ar[i,t] <- 1
        } else  {
          j <- apairs_f[i,t]
          ar[i,t] <- 2 + recap_m[j,t]
        }
      }
    }
  }
  
  ar[is.na(ar)] <- 0
  az[is.na(az)] <- 0
  
  # Return JAGS/NIBMLE (and true) Data
  model_data <- list(
    
    # Known data 
    n               = n, # Number of animals sampled
    k               = k,  # Number of occasions
    nf              = nf, # Number of females
    nm              = nm, # Number of males
    sex             = sex, # Sex of sampled individuals
    initial_entry   = initial_entry, # When did they enter the population
    first_capture_m = first_capture_m, # When did they enter the population
    first_capture_f = first_capture_f, # When did they enter the population
    # Latent States (true values - hidden in real data)
    mating_f        = mating_f, # Mate status of females at t (+ dummy)
    mating_m        = mating_m, # Mate status of males at t (+ dummy)
    pairs_f         = pairs_f, # partners of females 
    pairs_m         = pairs_m, # partners of males
    pairs           = pairs, # pair histories
    histories       = histories, # number of occasions a pair occurred
    sf              = sf, # true survival of females
    sm              = sm, # true survival of males
    repartner       = repartner, # whether partner from time t-1 was repicked (female x time)
    
    # Observed /Inferred states (Missing Values are possible)
    psi                 = psi, # Pairs that may exist (not excluded due to already formed pairs)
    next_partner_matrix = next_partner_matrix,
    last_partner_matrix = last_partner_matrix,
    ar_known            = ar,
    az_known            = az,
    af                  = rbind(af[1:nf,1:k],rep(0,k)),  # Female Survival with missing values
    am                  = rbind(am[1:nm,1:k],rep(0,k)),  # Male Survival with missing values
    apairs              = apairs, # Joint Pairs Matrices (array across time)
    apairs_f            = apairs_f,
    apairs_m            = apairs_m[1:nm, 1:k],
    arepartner          = arepartner, # repartner with inferred states 
    na_repartner        = 1*is.na(arepartner), # na spots for special cases repartner with inferred states  
    amating_f           = rbind(amating_f[1:nf,1:k],rep(0,k)), # Mating Status Females at T
    amating_m           = rbind(amating_m[1:nm,1:k],rep(0,k)),  # Mating Status Males at T
    recap_f             = rbind(recap_f[1:nf,1:k],rep(1,k)), #recap_f[1:nf,1:k], # Observed Recapture of Females
    recap_m             = rbind(recap_m[1:nm,1:k],rep(1,k))  # Observed Recapture of Males
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
  initial_entry <- c(model_data$first_capture_f,model_data$first_capture_m)
  
  
  # Number of females and males
  nf <- model_data$nf
  nm <- model_data$nm
  
  # Drop animals observed at final occasion (condition on first capture)
  if(any(initial_entry == model_data$k)){
    observed_at_k <- which(initial_entry == model_data$k)
    x <- x[-observed_at_k,]
    a <- a[-observed_at_k,]
    initial_entry <- initial_entry[-observed_at_k]
    
    nf2 <- nf
    nm2 <- nm
    
    # Lower amount of individuals for each one dropped
    for(i in observed_at_k){
      if(i > nf){
        nm2 <- nm2- 1
      } else {
        nf2 <- nf2 - 1
      }
    }
    
  } else {
    nf2 <- nf
    nm2 <- nm
  }
  
  # Store results in list
  results <- list(n = nf2 + nm2,
                  k = model_data$k,
                  female = c(rep(1, nf2),rep(0, nm2)), 
                  initial_entry = initial_entry,
                  x = x,
                  a = a)
  
  # Return Standard CJS Data
  return(results)
}

#----------------------------------------------------------------------------------------------------------------------