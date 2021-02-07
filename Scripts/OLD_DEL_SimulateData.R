###Setup###
rm(list=ls())
options(scipen = 999)
##Libraries
libs <- c("tidyverse","openxlsx","data.table")

##If package is not installed then do so
for(pkg in libs){
  if(!pkg %in% rownames(installed.packages()))
    install.packages(pkg)
}

##Attach our libraries
lapply(libs,require,character.only=TRUE)

##Inline object concatenation (like python)
`%+%` <- function(a,b) paste0(a,b)

setwd("/home/alexd/Documents/Research/Chapter 2 - Dyads/")

##Set Working Directory and reference paths
path2dat <-getwd() %+% "/Data/"
path2out <- getwd() %+% "/Output/"
path2scripts <- getwd() %+% "/Scripts/"

###Simulate gamma-rho-CJS Data ###
##Sample Size (individuals - will be split into groups below)
n <- 100
##Sampling Occasions
k <- 15

########################
#Breeding Probabilities#
########################
delta <- runif(1:k,0.5,1)
delta <- rep(0.8,k)

breed <- function(previous_state,sex,j){
  
  # Filters
  alive_male <- (previous_state == 3 && sex == "M")
  alive_female <- (previous_state == 2 && sex == "F")
  alive_partnered <- (previous_state == 4)
  
  mask <- (alive_male | alive_female | alive_partnered)
  
  # Check if breeding occurs
  if(mask){
    obs <- rbinom(1,1,delta[j])
  } else {
    obs <- 0
  }
  return(obs)
}

########################
#Survival Probabilities#
########################

##Marginal Probability of a female and male surviving from time j to j + 1
phi.f <- runif(1:k,0.5,1)
phi.m <- runif(1:k,0.5,1)

phi.f <- rep(0.8,k)
phi.m <- rep(0.8,k)

##Binomial Standard Deviation 
sig.phi.f <- sqrt(phi.f*(1-phi.f))
sig.phi.m <- sqrt(phi.m*(1-phi.m))

##Upper and Lower bounds for covariance parameter of survivorship based on phif and phim
gamma_upper_bound <- pmin(1,
                          (phi.f-phi.f*phi.m)/(sig.phi.f*sig.phi.m),
                          (phi.m-phi.f*phi.m)/(sig.phi.f*sig.phi.m))

gamma_lower_bound <- pmax(-1,
                          (-phi.f*phi.m)/(sig.phi.f*sig.phi.m),
                          (phi.f+phi.m - phi.f*phi.m - 1)/(sig.phi.f*sig.phi.m))

##Simulate Gamma (covariance) based on range above
gamma <- runif(1:k,gamma_lower_bound,gamma_upper_bound)
gamma <- rep(0.5,k)

###Joint Probability Distn for survivorship 
phi.mf <- gamma * sig.phi.m * sig.phi.f + (phi.f*phi.m)
phi.f0 <- phi.f - phi.mf
phi.m0 <- phi.m - phi.mf
phi.00 <- 1 - phi.f0 - phi.m0 - phi.mf 


#Survival 
survival <- function(previous_state,partner_id,sex,j){
  if(sex == "M"){
    if(previous_state == 4 || previous_state == 3){
      if(partner_id == 0){
        next_state <- 
          which(rmultinom(1,1,c((1-phi.m[j]),0,phi.m[j],0)) == 1)
      } else {
        next_state <-
          which(rmultinom(1,1,c(phi.00[j],phi.f0[j],phi.m0[j],phi.mf[j])) == 1)
      }
    } else {
      next_state = 1
    }
  } else if(sex == "F"){
    if(previous_state == 4 || previous_state == 2){
      if(partner_id == 0){
        next_state <- 
          which(rmultinom(1,1,c((1-phi.f[j]),phi.f[j],0,0)) == 1)
      } else {
        next_state <-
          which(rmultinom(1,1,c(phi.00[j],phi.f0[j],phi.m0[j],phi.mf[j])) == 1)
      }
    } else {
      next_state = 1
    }
  }
  return(next_state)
} 


########################
#Capture Probabilities#
########################

##Marginal recapture Probability of a female and male at j
p.f <- runif(1:k,0.4,1)
p.m <- runif(1:k,0.4,1)

p.f <- rep(0.8,k)
p.m <- rep(0.8,k)


##Binomial Standard Deviation 
sig.p.f <- sqrt(p.f*(1-p.f))
sig.p.m <- sqrt(p.m*(1-p.m))

##Upper and Lower bounds for covariance parameter of recapture based on phif and phim
rho_upper_bound <- pmin(1,
                        (p.f-p.f*p.m)/(sig.p.f*sig.p.m),
                        (p.m-p.f*p.m)/(sig.p.f*sig.p.m))

rho_lower_bound <- pmax(-1,
                        (-p.f*p.m)/(sig.p.f*sig.p.m),
                        (p.f+p.m - p.f*p.m - 1)/(sig.p.f*sig.p.m))

##Simulate Rho (covariance) based on range above
rho <- runif(1:k,rho_lower_bound,rho_upper_bound)
rho <- rep(0.25,k)

###Joint Probability Distn for recapture 
p.mf <- rho*sig.p.m*sig.p.f + p.f*p.m
p.f0 <- p.f - p.mf
p.m0 <- p.m - p.mf
p.00 <- 1 - p.m0 - p.f0 - p.mf


recapture <- function(current_state,j){
  if(current_state == 3){
    obs <- 
      which(rmultinom(1,1,c((1-p.m[j]),0,p.m[j],0)) == 1)  
  } else if(current_state == 2){ 
    obs <- 
      which(rmultinom(1,1,c((1-p.f[j]),p.f[j],0,0)) == 1)
  } else if(current_state == 4) { 
    obs <-
      which(rmultinom(1,1,c(p.00[j],p.f0[j],p.m0[j],p.mf[j])) == 1)
  } else { 
    obs <- 1
  } 
  return(obs)
} 


###########################
#Initialize Mating Pairs###
###########################

#Pair formation parameters
prob.monog <- rep(0.8,k)
prob.mating <- rep(0.8,k)

##Sample probability for gender information based on number of individuals n
prob.female <- 0.5

###Sort Population over entire study and compute initial entry time
population <- data_frame(animal_id = 1:n,
                         sex = ifelse(rbinom(n,1,prob.female)==1,"F","M")) %>%
  mutate(initial_entry=as.integer(sample(c(1:round(k/2)),size=n,replace=T)))

#Assign Initial Pairs based on entry time (partner must enter together)
init_pop_list <- list()

for(j in unique(population$initial_entry)){
  
  #Split by gender and initial entry
  temp.m <- population %>% filter(sex=="M") %>% filter(initial_entry == j)
  temp.f <- population %>% filter(sex=="F") %>% filter(initial_entry == j)
  
  #Assign Pairs to males (categorical distribution)
  temp.m <- temp.m %>% 
    mutate(partner_id = as.integer(
      sample(
        c(rep(0,max(nrow(temp.m)-nrow(temp.f),0)),unique(temp.f$animal_id)),
        size=nrow(temp.m),
        replace = F
      ) 
    )  
    )
  
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


###########################
#Partnership Probabilities#
###########################

#Regression Tensors

#Number of males, females, covariates
n.males <- population %>% filter(sex == "M") %>% nrow() 
n.females <- population %>% filter(sex == "F") %>% nrow() 
n.covariates <- 2

#Generate Zero Vectors 
z.males <- rep(0,n.males)
z.females <- rep(0,n.females)
z.covariates <- rep(0,n.covariates-1)
z.time <- rep(0,k)

#Organize vectors
#Male Tensor
z.vec.m <- c(z.females,z.covariates,z.time,z.males)
dim.vec.m <- c(n.females,n.covariates,k,n.males)

#Female Tensor
z.vec.f <- c(z.males,z.covariates,z.time,z.females)
dim.vec.f <- c(n.males,n.covariates,k,n.females)

#Generate Empty Tensors
reg.tensor.m <- array(z.vec.m,dim = dim.vec.m) %>% provideDimnames()
reg.tensor.f <- array(z.vec.f,dim = dim.vec.f) %>% provideDimnames()

#Regression Coefficients 
beta.m <- c(0.5,3)
beta.f <- c(0.5,3)

#Rename Male and Female Index on Tensors
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

male.cov.id <- c("Intercept","History") 
female.cov.id <- c("Intercept","History") 

#Apply the names

#Male Choice Tensor
dimnames(reg.tensor.m)[[4]] <- male.id
dimnames(reg.tensor.m)[[3]] <- 1:k
dimnames(reg.tensor.m)[[2]] <- female.cov.id
dimnames(reg.tensor.m)[[1]] <- female.id

#Female Choice Tensor
dimnames(reg.tensor.f)[[4]] <- female.id
dimnames(reg.tensor.f)[[3]] <- 1:k
dimnames(reg.tensor.f)[[2]] <- male.cov.id
dimnames(reg.tensor.f)[[1]] <- male.id

##########################
#Set up Static Parameters#
##########################

#Set up Intercept
add_intercept <- function(reg.tensor){
  time <- dim(reg.tensor)[3]
  chooser <- dim(reg.tensor)[4]
  
  for(i in 1:time){
    for(j in 1:chooser){
      reg.tensor[,1,i,j] <- 1 
    }
  }
   
  return(reg.tensor)
  
} 


reg.tensor.f <- add_intercept(reg.tensor.f)
reg.tensor.m <- add_intercept(reg.tensor.m)

#########################################
#Functions for Time-dependent parameters#
#########################################

#Update History
#Make sure current iteration of results is up to date
update_hist <- function(reg.tensor,gender,t,results){
  
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
      
      #print("row "  %+% i %+% " and " %+% chooser %+% " is the chooser" )
      #print("row "  %+% i %+% " and " %+% chosen %+% " is the chosen" )
      
      #If there is no partner skip to the next iteration
      if(chosen == "0"){
        next
      #Otherwise from this time point into the future increase the partnership hist by 1
      } else {
        for(j in t:k){
          reg.tensor[chosen,2,j,chooser] <- reg.tensor[chosen,2,j,chooser] + 1
          #print("Time is " %+% j %+% " tensor is " %+% reg.tensor[chosen,2,j,chooser])
          }
      }
    }
  }
  return(reg.tensor)
}

#Initialize Regression tensor with first time period 
reg.tensor.f <- update_hist(reg.tensor.f,"F",1,results)
reg.tensor.m <- update_hist(reg.tensor.m,"M",1,results)

###############################################
#Function to Generate dyad choice Probabilities
###############################################

#Softmax Function
softmax <- function(vector){
  p <- exp(vector)/sum(exp(vector))
  return(p)
}

#Function to Generate Pair Probabilities
gen_psi <- function(reg.tensor,individual,time,covariates){
  individual <- as.character(individual)
  unweighted <- covariates %*% t(reg.tensor[,,time,individual])
  psi <- softmax(unweighted)
  return(as.vector(psi))
} 


#Choice Probability Tensors

#Organize vectors
#Male Tensor
z.vec.m <- c(z.males,z.females,z.time)
dim.vec.m <- c(n.males,n.females,k)

#Female Tensor
z.vec.f <- c(z.females,z.males,z.time)
dim.vec.f <- c(n.females,n.males,k)

#Generate Empty Tensors
prob.tensor.m <- array(z.vec.m,dim = dim.vec.m) %>% provideDimnames()
prob.tensor.f <- array(z.vec.f,dim = dim.vec.f) %>% provideDimnames()

#Add Names
#Male Choice Tensor
dimnames(prob.tensor.m)[[1]] <- male.id
dimnames(prob.tensor.m)[[2]] <- female.id
dimnames(prob.tensor.m)[[3]] <- 1:k

#Female Choice Tensor
dimnames(prob.tensor.f)[[1]] <- female.id
dimnames(prob.tensor.f)[[2]] <- male.id
dimnames(prob.tensor.f)[[3]] <- 1:k


update_prob_tensor <- function(prob.tensor,reg.tensor,t,covariates){
  choosers <- rownames(prob.tensor[,,t])
  for(i in choosers){
    prob.tensor[i,,t] <- gen_psi(reg.tensor,i,t,covariates)
  }
  return(prob.tensor)
}

prob.tensor.f <- update_prob_tensor(prob.tensor.f,reg.tensor.f,1,beta.m)
prob.tensor.m <- update_prob_tensor(prob.tensor.m,reg.tensor.m,1,beta.f)

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
  
  ntt$mating <- mapply(breed,ntx$surv.y,ntx$sex,MoreArgs = list(j=i+1))
  
  ############################
  #Seperate Variables present#
  ############################
  
  itt <- ntt %>% filter(mating == 1,sex == "M") 
  jtt <- ntt %>% filter(mating == 1,sex == "F")  
  ktt <- ntt %>% filter(mating == 0)
  
  ################
  #Mate Selection#
  ################
  
  #Order the choosing gender (prob can be updated with eigenvector)
  chosers <- if(length(jtt$animal_id) > 1) base::sample(jtt$animal_id,length(jtt$animal_id),replace=FALSE) else if(length(jtt$animal_id)==1) jtt$animal_id else 0
  chosen <- itt$animal_id
  
  #Truncate Choosers
  if(length(chosen) < length(chosers)){
    chosers <- chosers[1:length(itt$animal_id)]
  }
  
  #Remove Dead Choosers
  
  
  #Check if there are any breeders
  if(length(chosen) != 0 && length(chosers)!=0 && chosers != 0){
    #Simulate Selection
    for(j in chosers){
      
      #Extract Probability and reweight
      psi <- prob.tensor.f[as.character(j),,i]
      psi <- psi[as.character(chosen)]/sum(psi[as.character(chosen)])
      
      #Make Partner Choice
      choice_index <- which(rmultinom(1,1,psi)==1)
      choice <- names(psi)[choice_index]
      #Assign Partners
      ntt[ntt$animal_id == as.integer(choice),]$partner_id <- j #Add to choice
      ntt[ntt$animal_id == j,]$partner_id <- as.integer(choice) #Add to chosen
      #Remove chosen from options
      chosen <- chosen[!chosen %in% choice]
      print("Female " %+% j %+% " has picked male " %+% choice %+% " at time " %+% i)
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
  
  for(j in unique(ntt$animal_id)){
    
    #Make sure not to re-model partners who already have their status
    if(!is.na(ntt[ntt$animal_id == j,7])){next}
    
    #Extract Key Variables
    previous_state <- nt %>% filter(animal_id == j) %>% select(surv) %>% as.integer()
    partner_id <- ntt %>% filter(animal_id == j) %>% select(partner_id) %>% as.integer()
    sex <- ntt %>% filter(animal_id == j) %>% select(sex) %>% as.character()
    
    #Simulate Survival, Recapture, and Confounding 
    ntt[ntt$animal_id == j,7] <- current_state <- survival(previous_state,partner_id,sex,i+1)
    ntt[ntt$animal_id == j,8] <- current_observation  <- recapture(current_state,i+1)
    ntt[ntt$animal_id == j,9] <- current_surv_conf <- ifelse(current_observation==1,NA,current_observation)
    ntt[ntt$animal_id == j,10] <- current_partner_conf <-  ifelse(current_observation==4,partner_id,NA)
    
    print("Individual " %+% j)
    print("Current state= " %+% current_state)
    print("Current Obs= " %+% current_observation)
    print("Current state conf = " %+% current_surv_conf)
    print("Current partner conf= " %+% current_partner_conf)
    
    #Update PartnerID
    if(partner_id != 0){
      #Simulate Survival, Recapture, and Confounding 
      ntt[ntt$animal_id == partner_id,7] <- current_state 
      ntt[ntt$animal_id == partner_id,8] <- current_observation  
      ntt[ntt$animal_id == partner_id,9] <- current_surv_conf
      ntt[ntt$animal_id == partner_id,10]  <- ifelse(current_observation==4,j,NA)
    }
    
    rm(previous_state,partner_id,sex,current_state,current_observation)
    
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
  #Update Regression Tensors#
  ###########################
  
  reg.tensor.f <- update_hist(reg.tensor.f,"F",i+1,results)
  reg.tensor.m <- update_hist(reg.tensor.m,"M",i+1,results)
  prob.tensor.f <- update_prob_tensor(prob.tensor.f,reg.tensor.f,i+1,beta.m)
  prob.tensor.m <- update_prob_tensor(prob.tensor.m,reg.tensor.m,i+1,beta.f)

}


########################
#Add Individual History#
########################

#Individual survival

ind_surv_fn <- function(surv,gender){
  
  if(is.na(surv)){
    state <- NA
    
  } else {
    #Survival Conditions
    case1 <- surv == 4
    case2 <- surv == 3 && gender == "M"
    case3 <- surv == 2 && gender == "F"
    
    #Death Conditions
    case4 <- surv == 1
    case5 <- surv == 2 && gender == "M"
    case6 <- surv == 3 && gender == "F"
    
    #Unknown
    case7 <- is.na(surv)
    
    #Apply Conditions
    if(case1 || case2 || case3){
      state <- 1
    } else if(case4 || case5 || case6){
      state <- 0 
    }
  }
  
  return(state)
}

ind_recap_fn <- function(recapture,gender){
  
  #Survival Conditions
  case1 <- recapture == 4
  case2 <- recapture == 3 && gender == "M"
  case3 <- recapture == 2 && gender == "F"
  
  #Death Conditions
  case4 <- recapture == 1
  case5 <- recapture == 2 && gender == "M"
  case6 <- recapture == 3 && gender == "F"
  
  
  #Apply Conditions
  if(case1 || case2 || case3){
    obs <- 1
  } else if(case4 || case5 || case6){
    obs <- 0 
  }
  
  return(obs)
}

#Zero Variables
results$ind_surv <- rep(0,nrow(results))
results$ind_recapture <- rep(0,nrow(results))
results$ind_surv_conf <- rep(0,nrow(results))

#Apply Individual Survival
for(i in 1:nrow(results)){
  surv <- results$surv[i]
  gender <- results$sex[i]
  recap <- results$recapture[i]
  ind_surv_conf <- results$surv_confounded[i]
  
  results[i,11] <- ind_surv_fn(surv,gender)
  results[i,12] <- ind_recap_fn(recap,gender)
  results[i,13] <- ifelse(ind_recap_fn(recap,gender)==1,1,NA)
}

#Fill in Inferred Survival
for(i in unique(results$animal_id)){
  index <-  which(results[results$animal_id == i,13]==1)
  start <- min(index)
  end <- max(index)
  for(j in start:end){
    results[results$animal_id == i,13][j] <- 1
  }
}

###########################
#Ready Data for simulation#
###########################

format_CR_data <- function(results,reg.tensor.f,reg.tensor.m,n,k){
  #Mate History
  V <- matrix(results$partner_id,ncol=k,nrow=n,byrow=T) 
  v0 <- matrix(results$partner_id_confounded,ncol=k,nrow=n,byrow=T) 
  #Survival History
  S <- matrix(results$surv,ncol=k,nrow=n,byrow=T) 
  S0 <- matrix(results$surv_confounded,ncol=k,nrow=n,byrow=T) 
  S_ind <- matrix(results$ind_surv,ncol=k,nrow=n,byrow=T) 
  S0_ind <- matrix(results$ind_surv_conf,ncol=k,nrow=n,byrow=T) 
  #Recapture History
  RC <- matrix(results$recapture,ncol=k,nrow=n,byrow=T) 
  RC_ind <- matrix(results$ind_recapture,ncol=k,nrow=n,byrow=T)
  #Mate Attempt Choice
  M <- matrix(results$mating,ncol=k,nrow=n,byrow=T) 
  #Genders
  g <- results %>% select(animal_id,sex) %>% distinct()
  g <- g$sex
  #Initial Entry
  E <- results %>% select(animal_id,initial_entry) %>% distinct()
  E <- E$initial_entry
  #Animal ID
  ID <- results %>% select(animal_id) %>% distinct()
  ID <- ID$animal_id
  #Return Data
  data <- list(`mate_hist` = V,
               `mate_hist_unknown`=v0,
               `surv`=S,
               `surv_unknown`=S0,
               `S_ind`=S_ind,
               `S_ind_unknown`=S0_ind,
               `recapture`=RC,
               `ind_recapture`=RC_ind,
               `mating`=M,
               `gender`=g,
               `init_entry`=E,
               `ID`=ID,
               `covariates_m` = reg.tensor.m,
               `covariates_f` = reg.tensor.f,
               `ind` = n,
               `time`=k)
  
  return(data)
}

dat <- format_CR_data(results,reg.tensor.f,reg.tensor.m,n,k)

extract_cjs <- function(data){
  #CJS Data
  a <- data$S_ind_unknown
  x <- data$ind_recapture
  first <- data$init_entry
  k <- data$time
  n <- data$ind
  #Posterior generating function
  PGF <- function(n){
    phi <- runif(n)
    p <- runif(n)
    return(list("phi"=phi,"p"=p))
  }
  #
  cjs_dat <- list("PGF"=PGF,"k"=k,"n"=n,"a"=a,"x"=x,"first"=first)
  return(cjs_dat)
}

cjs_dat <- extract_cjs(dat)

saveRDS(dat,path2out %+% "full_dat.rds")
saveRDS(cjs_dat,path2out %+% "cjs_dat.rds")
