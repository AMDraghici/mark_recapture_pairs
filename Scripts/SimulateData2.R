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

mate <- function(previous_state, sex, j, delta){
  
  # Chose to mate at time t? Previous state is 0/1 here
  obs <- rbinom(1, 1, delta[j]*previous_state)

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

construct_init_entry <- function(n, k, random = F, inits = NULL){
  
  # Select initial entry randomly?
  if(random == T){
    initial_entry = sort(as.integer(sample(c(1:round(k/2)),size=n,replace=T)))
  } else {
    # If not do we have pre-specified values?
    if(!is.null(inits)){
      initial_entry <- inits 
    # If not just do fixed intervals across n and k 
    } else {
      initial_entry <- sort(rep(1:(k-1),ceiling(n/(k-1))))
      initial_entry <- initial_entry[1:n]
    }
  }
  # Return initial entries
  return(initial_entry)
}

#Skeleton Matrix for mate choice at k 
construct_mated <- function(n,k){
  # Mate status at time t for each animal 
  mating <- matrix(NA,nrow=length(sex),ncol = k)
  # Return Matrix
  return(mating)
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
  coef_list <- list(intercept =intercept, histories =histories)
  return(coef_list)
}

# Survival matrices for males, females, and 
construct_survival <- function(sex, k){
  #Females
  sf <- matrix(NA, nrow = length(sex[sex == "F"]), ncol = k)
  #Males
  sm <- matrix(NA, nrow = length(sex[sex == "M"]), ncol = k)
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

# 6. Update Data structures from [t,t+1)------------------------------------------------------------------------

# 7. Simulate Data ---------------------------------------------------------------------------------------------


simulate_cr_data <- function(n,
                             k, 
                             prop.female,
                             phi.f, 
                             phi.m, 
                             gam, 
                             p.f, 
                             p.m, 
                             rho, 
                             beta0,
                             beta1,
                             rand_sex = F,
                             rand_init = F){
  
  # Generate SKeleton Data Structures
  sex <- construct_sexes(n, prop.female, rand_sex)
  initial_entry <- construct_init_entry(n, k, rand_init) 
  mating <- construct_mated(n, k)
  pairs_list <- construct_pairs(n, k)
  pf <- pairs_list[["pf"]]
  pm <- pairs_list[["pm"]]
  pairs <- pairs_list[["pairs"]]
  coef_list <- construct_coef(sex, k)
  intercept <- coef_list[["intercept "]]
  histories <- coef_list[["histories"]]
  survival_list <- construct_survival(sex, k)
  sf <- survival_list[["sf"]]
  sm <- survival_list[["sm"]]
  spair <- survival_list[["spair"]]
  rpair <- construct_recapture(n,k)
  
  # Initialize Data 
  
  # Simulate Fates
  
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














initiate_population <- function(n,k, prop.female = 0.5){

  ###Sort Population over entire study and compute initial entry time
  population <- data_frame(animal_id = 1:n,
                           sex = construct_gender(n, prop.female = prop.female)
                           ) %>%
    mutate(initial_entry = as.integer(sample(c(1:round(k/2)),size=n,replace=T)))
  
  #Assign Initial Pairs based on entry time (partner must enter together)
  init_pop_list <- list()
  
  # Index of initial entries
  initial_entries <- sort(unique(population$initial_entry))
  
  for(j in initial_entries){
    
    #Split by gender and initial entry
    temp.m <- population %>% filter(sex=="M") %>% filter(initial_entry == j)
    temp.f <- population %>% filter(sex=="F") %>% filter(initial_entry == j)
    
    #Assign Pairs to males (categorical distribution)
    temp.m <- temp.m %>% 
      mutate(partner_id = as.integer(
        sample(
          c(rep(0,max(nrow(temp.m)-nrow(temp.f),0)),
            unique(temp.f$animal_id)
            ),
          size=nrow(temp.m),
          replace = F)
        ))
    
    #Update Female data to reflect assignment
    temp.f <- temp.f %>% 
      left_join(temp.m,by=c("animal_id"="partner_id","initial_entry")) %>% 
      select(-sex.y) %>% 
      rename(partner_id = animal_id.y,
             sex = sex.x) %>% 
      replace(., is.na(.), as.integer(0))
    
    #Join tables and store in list indexed by initial entry
    init_pop_list[[j]] <- rbind(temp.f,temp.m)
    
    rm(temp.f,temp.m)
    
  }  
  
  #Compute population matrix
  population <- bind_rows(init_pop_list)
  
  #Generate Empty Results Data Frame
  results <- merge(population$animal_id,1:k) %>% 
    rename("animal_id"="x","time"="y") %>%  
    inner_join(population,by=c("animal_id")) %>% 
    arrange(animal_id) %>% 
    mutate(partner_id = ifelse(time==initial_entry,partner_id,
                               ifelse(time<initial_entry,0,NA)),
           mating = ifelse(initial_entry==time,ifelse(partner_id==0,0,1),
                           ifelse(initial_entry>time,0,NA)),
           surv = ifelse(initial_entry==time,ifelse(partner_id==0,
                                                    ifelse(sex=="M",3,2),4),
                         ifelse(initial_entry>time,1,NA)),
           recapture = ifelse(initial_entry==time,ifelse(partner_id==0,
                                                         ifelse(sex=="M",3,2),4),
                              ifelse(initial_entry>time,1,NA)),
           surv_confounded = ifelse(initial_entry==time,surv,
                                    ifelse(initial_entry>time,1,NA)),
           partner_id_confounded = ifelse(initial_entry==time,partner_id,
                                          ifelse(initial_entry>time,0,NA))) 
  
  # Return Results Template and Population Index
  return(list(population = population, 
              results = results))
}

###########################
#Partnership Probabilities#
###########################

#########################################
#Functions for Time-dependent parameters#
#########################################

#Set up Intercept
add_intercept <- function(reg.array){
  time <- dim(reg.array)[3]
  chooser <- dim(reg.array)[4]
  
  for(i in 1:time){
    for(j in 1:chooser){
      reg.array[,1,i,j] <- 1 
    }
  }
  
  return(reg.array)
  
} 

#Update History
#Make sure current iteration of results is up to date
update_hist <- function(reg.array,gender,t,results){
  
  #Pull Current Partnership History
  temp.hist <- results %>% 
    filter(sex==gender,time==t,initial_entry <= t) %>% 
    select(animal_id,partner_id)
  
  #If there are individuals then update
  if(nrow(temp.hist)!=0){
    #Sort through each entity
    for(i in 1:nrow(temp.hist)){
      #Pick the chooser and choosen 
      chooser <- as.character(temp.hist[i,1])
      chosen <- as.character(temp.hist[i,2])
      
      #If there is no partner skip to the next iteration
      if(chosen == "0"){
        next
        #Otherwise from this time point into the future increase the partnership hist by 1
      } else {
        for(j in t:k){
          reg.array[chosen,2,j,chooser] <- reg.array[chosen,2,j,chooser] + 1
        }
      }
    }
  }
  return(reg.array)
}

###############################################
#Function to Generate dyad choice Probabilities
###############################################

#Function to Generate Pair Probabilities
gen_psi <- function(reg.array,individual,time,covariates){
  individual <- as.character(individual)
  unweighted <- covariates %*% t(reg.array[,,time,individual])
  psi <- softmax(unweighted)
  return(as.vector(psi))
} 

# Update probability Array
update_prob_array <- function(prob.array,reg.array,t,covariates){
  choosers <- rownames(prob.array[,,t])
  for(i in choosers){
    prob.array[i,,t] <- gen_psi(reg.array,i,t,covariates)
  }
  return(prob.array)
}

construct_data_template <- function(n, 
                                    k, 
                                    prop.female = 0.5, 
                                    n.covariates = 2, 
                                    beta.m = c(0.5,3),
                                    beta.f = c(0.5,3),
                                    effect.name.m = c("Intercept","History"),
                                    effect.name.f = c("Intercept","History")
                                    ){
  
  #Pull population template
  init_pop <- initiate_population(n,k, prop.female)
  population <- init_pop[[1]]
  results <- init_pop[[2]]
  
  #Number of males, females
  n.males <- population %>% filter(sex == "M") %>% nrow() 
  n.females <- population %>% filter(sex == "F") %>% nrow() 
  
  #Generate Zero Vectors 
  z.males <- rep(0,n.males)
  z.females <- rep(0,n.females)
  z.covariates <- rep(0,n.covariates-1)
  z.time <- rep(0,k)
  
  #Organize vectors
  #Male Array
  z.vec.m <- c(z.females,z.covariates,z.time,z.males)
  dim.vec.m <- c(n.females,n.covariates,k,n.males)
  
  #Female Array
  z.vec.f <- c(z.males,z.covariates,z.time,z.females)
  dim.vec.f <- c(n.males,n.covariates,k,n.females)
  
  #Generate Empty Arrays
  reg.array.m <- array(z.vec.m,dim = dim.vec.m) %>% provideDimnames()
  reg.array.f <- array(z.vec.f,dim = dim.vec.f) %>% provideDimnames()
  
  #Rename Male and Female Index on Arrays
  male.id <- population %>% 
    filter(sex == "M") %>% 
    select(animal_id) %>% 
    t() %>% 
    as.vector() %>% 
    sort()
  
  female.id <- population %>%
    filter(sex == "F") %>% 
    select(animal_id) %>% 
    t() %>% 
    as.vector() %>% 
    sort()
  
  #Apply the names
  
  #Male Choice array
  dimnames(reg.array.m)[[4]] <- male.id
  dimnames(reg.array.m)[[3]] <- 1:k
  dimnames(reg.array.m)[[2]] <- effect.name.f
  dimnames(reg.array.m)[[1]] <- female.id
  
  #Female Choice array
  dimnames(reg.array.f)[[4]] <- female.id
  dimnames(reg.array.f)[[3]] <- 1:k
  dimnames(reg.array.f)[[2]] <- effect.name.m
  dimnames(reg.array.f)[[1]] <- male.id
  
  # Add intercepts
  reg.array.f <- add_intercept(reg.array.f)
  reg.array.m <- add_intercept(reg.array.m)
  
  #Initialize Regression array with first time period 
  reg.array.f <- update_hist(reg.array.f,"F",1,results)
  reg.array.m <- update_hist(reg.array.m,"M",1,results)
  
  #Choice Probability arrays
  
  #Organize vectors
  #Male array
  z.vec.m <- c(z.males,z.females,z.time)
  dim.vec.m <- c(n.males,n.females,k)
  
  #Female array
  z.vec.f <- c(z.females,z.males,z.time)
  dim.vec.f <- c(n.females,n.males,k)
  
  #Generate Empty arrays
  prob.array.m <- array(z.vec.m,dim = dim.vec.m) %>% provideDimnames()
  prob.array.f <- array(z.vec.f,dim = dim.vec.f) %>% provideDimnames()
  
  #Add Names
  #Male Choice array
  dimnames(prob.array.m)[[1]] <- male.id
  dimnames(prob.array.m)[[2]] <- female.id
  dimnames(prob.array.m)[[3]] <- 1:k
  
  #Female Choice array
  dimnames(prob.array.f)[[1]] <- female.id
  dimnames(prob.array.f)[[2]] <- male.id
  dimnames(prob.array.f)[[3]] <- 1:k

  prob.array.f <- update_prob_array(prob.array.f,reg.array.f,1,beta.m)
  prob.array.m <- update_prob_array(prob.array.m,reg.array.m,1,beta.f)
  
  # Gather important data
  out.list <- list(population = population,
       results = results, 
       reg.array.f = reg.array.f,
       reg.array.m = reg.array.m,
       prob.array.f = prob.array.f,
       prob.array.m = prob.array.m)
  
  # Return Results
  return(out.list)
}

simulate_cr_data <- function(n, 
                             k, 
                             prop.female = 0.5, 
                             n.covariates = 2, 
                             beta.m = c(0.5,3),
                             beta.f = c(0.5,3),
                             effect.name.m = c("Intercept","History"),
                             effect.name.f = c("Intercept","History"),
                             phi.m, 
                             phi.f,
                             gam,
                             p.m,
                             p.f,
                             rho,
                             delta
                             ){
 
   # Generate initial data 
   initial_data_list <- 
     construct_data_template(n,k,
                             prop.female,n.covariates,
                             beta.m,beta.f,
                             effect.name.m,effect.name.f)

   
   # Extract Results
   population <- initial_data_list[[1]]
   results <- initial_data_list[[2]]
   reg.array.f <- initial_data_list[[3]]
   reg.array.m <- initial_data_list[[4]]
   prob.array.f <- initial_data_list[[5]]
   prob.array.m <-  initial_data_list[[6]]
   
   ##################
   #Simulate History#
   ##################
   
   for(i in 1:(k-1)){
     
     #Split Individuals by Gender or Partner Status
     it <- results %>% filter(time == i,sex=="M",partner_id !=0,initial_entry <= i)
     jt <- results %>% filter(time == i,sex=="F",partner_id !=0,initial_entry <= i)
     kt <- results %>% filter(time == i,partner_id==0,initial_entry <= i)
     mt <- rbind(it,kt)
     nt <- rbind(it,jt,kt)  
     
     ntt <- results %>%
       filter(time == i+1,initial_entry <= i)
     
     ntx <- ntt %>% inner_join(nt,by=c("animal_id","sex","initial_entry"))
     
     ntt$mating <- mapply(breed,ntx$surv.y,ntx$sex,MoreArgs = list(j=i+1,delta=delta))
     
     ############################
     #Seperate Variables present#
     ############################
     
     itt <- ntt %>% filter(mating == 1,sex == "M") 
     jtt <- ntt %>% filter(mating == 1,sex == "F")  
     ktt <- ntt %>% filter(mating == 0)
     
     ################
     #Mate Selection#
     ################
     
     #Order the choosing gender 
     chosers <- if(length(jtt$animal_id) > 1) base::sample(jtt$animal_id,length(jtt$animal_id),replace=FALSE) else if(length(jtt$animal_id)==1) jtt$animal_id else 0
     chosen <- itt$animal_id
     
     #Truncate Choosers
     if(length(chosen) < length(chosers)){
       chosers <- chosers[1:length(itt$animal_id)]
     }

     #Check if there are any breeders
     if(length(chosen) != 0 && length(chosers)!=0 && chosers != 0){
       #Simulate Selection
       for(j in chosers){
         #Extract Probability and reweight
         psi <- prob.array.f[as.character(j),,i]
         psi <- psi[as.character(chosen)]/sum(psi[as.character(chosen)])
         #Make Partner Choice
         choice_index <- which(rmultinom(1,1,psi)==1)
         choice <- names(psi)[choice_index]
         #Assign Partners
         ntt[ntt$animal_id == as.integer(choice),]$partner_id <- j #Add to choice
         ntt[ntt$animal_id == j,]$partner_id <- as.integer(choice) #Add to chosen
         #Remove chosen from options
         chosen <- chosen[!chosen %in% choice]
       }
     }
     

     #Zero out Non-maters
     if(length(ntt[ntt$mating==0,]$partner_id) != 0){
       ntt[ntt$mating==0,]$partner_id <- 0
     }
     
     if(length(ntt[ntt$mating==1 & is.na(ntt$partner_id),]$partner_id) != 0){
       ntt[ntt$mating==1 & is.na(ntt$partner_id),]$partner_id <- 0
     }
     
     ########################
     #Survival and Recapture#
     ########################
     
     for(j in sort(unique(ntt$animal_id))){
       
       #Make sure not to re-model partners who already have their fates
       if(!is.na(ntt[ntt$animal_id == j,7])){next}
       
       #Extract Key Variables
       previous_state <- nt %>% filter(animal_id == j) %>% select(surv) %>% as.integer()
       partner_id_j <- ntt %>% filter(animal_id == j) %>% select(partner_id) %>% as.integer()
       sex <- ntt %>% filter(animal_id == j) %>% select(sex) %>% as.character()
       
       if(partner_id_j != 0){
         previous_state <- 4
       }
       
       #Simulate Survival, Recapture, and Confounding 
       ntt[ntt$animal_id == j,7] <- current_state <- survival(previous_state,i+1,phi.m, phi.f, gam)
       ntt[ntt$animal_id == j,8] <- current_observation  <- recapture(current_state,i+1, p.m, p.f, rho)
       ntt[ntt$animal_id == j,9] <- current_surv_conf <- ifelse(current_observation==1,NA,current_observation)
       ntt[ntt$animal_id == j,10] <- current_partner_conf <-  ifelse(current_observation==4,partner_id_j,NA)
       
       #Update PartnerID
       if(partner_id_j != 0){
         #Simulate Survival, Recapture, and Confounding 
         ntt[ntt$animal_id == partner_id_j,7] <- current_state 
         ntt[ntt$animal_id == partner_id_j,8] <- current_observation  
         ntt[ntt$animal_id == partner_id_j,9] <- current_surv_conf
         ntt[ntt$animal_id == partner_id_j,10]  <- ifelse(current_observation==4,j,NA)
       }
       
       rm(previous_state,partner_id_j,sex,current_state,current_observation)
       
     }
     
     #######################
     #Update Master Results#
     #######################
     
     #Filters
     animal_filter <- results$animal_id %in% ntt$animal_id
     time_filter <- results$time == i + 1
     entry_filter <- results$initial_entry < i + 1
     
     #Assign Mate Status
     results[animal_filter & time_filter & entry_filter,]$mating <- ntt$mating
     #Assign PartnerId
     results[animal_filter & time_filter & entry_filter,]$partner_id <- ntt$partner_id
     #Assign Survival
     results[animal_filter & time_filter & entry_filter,]$surv <- ntt$surv
     #Assign Recapture 
     results[animal_filter & time_filter & entry_filter,]$recapture <- ntt$recapture
     #Assign Confounded Surv
     results[animal_filter & time_filter & entry_filter,]$surv_confounded <- ntt$surv_confounded
     #Assign Confounded PartnerID
     results[animal_filter & time_filter & entry_filter,]$partner_id_confounded <- ntt$partner_id_confounded
     
     
     ###########################
     #Update Regression arrays#
     ###########################
     
     reg.array.f <- update_hist(reg.array.f,"F",i+1,results)
     reg.array.m <- update_hist(reg.array.m,"M",i+1,results)
     prob.array.f <- update_prob_array(prob.array.f,reg.array.f,i+1,beta.m)
     prob.array.m <- update_prob_array(prob.array.m,reg.array.m,i+1,beta.f)
   }
   
   out.list <- list(population,
                    results, 
                    reg.array.f)
   return(out.list)
}
