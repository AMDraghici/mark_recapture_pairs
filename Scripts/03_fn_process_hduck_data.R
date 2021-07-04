
gather_hq_data <- function(dat_dir, 
                            book1 = "Harlequin_Capture_Records_24Aug2020.xls", 
                            book2 = "HARD_Obs_10May95_08Sep20.xls"){
  #Storage
  dat1 <- list()
  dat2 <- list()
  
  # 9 sheets in book1
  for(i in 1:9){
    dat1[[i]] <-  read_excel(dat_dir %+% book1, sheet = i)
  }
 
  # 5 sheets in book2
  for(j in 1:7){
    dat2[[j]] <- read_excel(dat_dir %+% book2, sheet = j)
  }
  
  # Return two data lists
  return(list(book1 = dat1, book2 = dat2))
}

process_bandlist <- function(hq_data){
  BandList <- hq_data[["book1"]][[2]]  %>%
    rename(year_class = `Year Class`,
           first_plastic = `1st Plastic`,
           second_plastic=`2nd Plastic`,
           third_plastic = `3rd Plastic`,
           fourth_plastic = `4th Plastic`)
  return(BandList)
}

process_capture_records <- function(hq_data){
  
  # Grab Bandlist Information
  BandList <- process_bandlist(hq_data)
  
  # Grab Capture Records/Covariates/Assign partners and build history using Animal ID only
  CaptureRecords <- hq_data[["book1"]][[3]] %>%  
    mutate(Capture.Date = as.Date(`Banding Date`,origin = "1899-12-30"),
           year = year(Capture.Date),
           month = month(Capture.Date),
           day = mday(Capture.Date)) %>% 
    rename("animal_id" = "ID",
           "sex" = "Sex",
           "capture_date" = Capture.Date,
           "date_of_birth" = `Year Class`,
           "age_estimate" = "Estimate",
           "pair_status" = "Pairstat",
           "partner_id" = "Mateleg",
           "initial_entry" = "Type",
           "first_plastic" = `1st Plastic`,
           "age" = "Age",
           "year_class" = `Year Class`,
           "mother_id"  = MotherID,
           "bird_wt" = `Bird Weight`,
           "wing_crd" = `Wing Chord`,
           "tarsus_len" = `Tarsus Length`,
           "culmen_len" = `Culmen Length`) %>%
    arrange(animal_id,capture_date,initial_entry)  %>%
    left_join(na.omit(select(BandList,ID,first_plastic)),by=c("partner_id" = "first_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,second_plastic)),by=c("partner_id" = "second_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,third_plastic)),by=c("partner_id" = "third_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,fourth_plastic)),by=c("partner_id" = "fourth_plastic")) %>% 
    mutate(partner_id = pmax(ifelse(is.na(ID.x),0,ID.x),
                             ifelse(is.na(ID.y),0,ID.y),
                             ifelse(is.na(ID.x.x),0,ID.x.x),
                             ifelse(is.na(ID.y.y),0,ID.y.y)
                             ),
           time = year - 1993 + 1
    ) %>% 
    select(-ID.x,-ID.y,-ID.x.x,-ID.y.y) %>% 
    select(capture_date,year,month,day,animal_id,sex,partner_id,
           pair_status,initial_entry, age, mother_id, time, bird_wt,             
           wing_crd,tarsus_len, culmen_len, Area) %>% 
    filter(!is.na(capture_date)) %>%
    left_join(na.omit(select(BandList,ID,first_plastic)),by=c("mother_id" = "first_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,second_plastic)),by=c("mother_id" = "second_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,third_plastic)),by=c("mother_id" = "third_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,fourth_plastic)),by=c("mother_id" = "fourth_plastic")) %>% 
    mutate(mother_id = pmax(ifelse(is.na(ID.x),0,ID.x),
                             ifelse(is.na(ID.y),0,ID.y),
                             ifelse(is.na(ID.x.x),0,ID.x.x),
                             ifelse(is.na(ID.y.y),0,ID.y.y))) %>% 
    select(-ID.x,-ID.y,-ID.x.x,-ID.y.y)
  
  # Return Processed Data
  return(CaptureRecords)
}

build_cr_df <- function(hq_data){
  CaptureRecords <- process_capture_records(hq_data)
  
  
  year.dat <- select(CaptureRecords,year, time) %>% 
    distinct() %>% 
    arrange(time)
  
  #Capture-Recapture format
  cap.data <- merge(unique(CaptureRecords$animal_id),year.dat$time) %>% 
    rename(animal_id = x,time = y) %>% 
    left_join(year.dat,by="time") %>% 
    left_join(CaptureRecords,by=c("animal_id","year", "time")) %>% 
    mutate(recapture_individual = ifelse(is.na(capture_date),0,1),
           mated = ifelse((!is.na(partner_id) & partner_id != 0) | pair_status=="P" | pair_status == "FY", 1, #partner observed
                          ifelse(pair_status=="N"|pair_status=="YOY",0, #young of the year or status as unpaired
                                 ifelse(pair_status=="U",NA,NA)))) %>% #Unknown (either alone or in a group)
    arrange(animal_id,time)
  
  
  # Add initial Capture 
  initial <- cap.data %>% 
    filter(recapture_individual == 1) %>% 
    select(animal_id,time) %>% 
    group_by(animal_id) %>% 
    dplyr::summarize(initial_entry = min(time)) %>% 
    ungroup()
  
  # Add sex across all values 
  sex <- filter(cap.data,!is.na(sex)) %>% 
    select(animal_id,sex) %>% 
    distinct()
  
  # Apply merge 
  cap.data <- cap.data %>% 
    select(-initial_entry,-sex) %>% 
    left_join(initial,by="animal_id") %>%
    left_join(sex,by="animal_id") %>% 
    select(-capture_date,-month,-day,-year,-pair_status) %>% 
    mutate(mated = ifelse(initial_entry > time,0,mated),
           partner_id = ifelse(initial_entry > time,0,ifelse(is.na(mated),NA,partner_id)),
           surv_individual_confounded = ifelse(recapture_individual == 0,ifelse(time<=initial_entry,1,NA),1),
           age = ifelse(age == "HY", 0, ifelse(age == "AHY" | age == "SY" | age == "ASY", 1, 1))) %>% 
    select(animal_id,time,initial_entry,
           sex,mated,partner_id, 
           age, mother_id, Area,
           bird_wt, wing_crd, tarsus_len, culmen_len,
           recapture_individual,surv_individual_confounded) %>% 
    distinct()
  
  # Drop duplicate captures for females
  #cap.data <- cap.data %>% group_by(animal_id, time) %>% slice(1) %>% ungroup() %>% as.data.frame()
  cap.data <- cap.data %>% group_by(animal_id, time) %>% summarize(
    animal_id = unique(animal_id),
    time = unique(time),
    initial_entry = unique(initial_entry),
    sex = unique(sex),
    mated = max(mated,na.rm=T),
    partner_id = max(partner_id, na.rm =T),
    age =  max(age,na.rm=T),
    mother_id = max(mother_id, na.rm =T),
    Area = max(Area, na.rm = T),
    bird_wt = max(bird_wt, na.rm = T),
    wing_crd = max(wing_crd, na.rm=T),
    tarsus_len = max(tarsus_len, na.rm =T),
    culmen_len = max(culmen_len, na.rm=T),
    recapture_individual = max(recapture_individual, na.rm=T),
    surv_individual_confounded = max(surv_individual_confounded, na.rm=T)
  ) %>% ungroup() %>% as.data.frame()
  
  cap.data[cap.data == -Inf] <- NA
  return(cap.data)
} 

add_implied_states <- function(cap.data){
  
  # Animal IDs and amount of time 
  animal_ids <- sort(unique(cap.data$animal_id))
  k <- max(unique(cap.data$time))
  
  # Update each animal in the dataset
  for(i in animal_ids){
    # Survival ----------------------------------------------------------------------
    
    # Survival states of animal i 
    surv_i <- cap.data[cap.data$animal_id == i,"surv_individual_confounded"]
    # When were they last seen alive
    last_alive <- max(which(surv_i == 1))
    # If they were last alive at time 1 or later then populate implied states (otherwise nothing)
    if(length(last_alive) !=0 && last_alive >= 1 & last_alive <= k){
      cap.data[cap.data$animal_id == i,"surv_individual_confounded"][1:last_alive] <- 1
    }
    
    # Age ----------------------------------------------------------------------------
    
    # Age states of animal i 
    age_i <- cap.data[cap.data$animal_id == i,"age"]
    # When were they first seen as older than 1 
    first_adult <- min(which(age_i == 1))
    if(is.infinite(first_adult)){
      first_adult <- 1
    }
    # ...or was their hatch year observed? 
    hatch_yr <- min(which(age_i == 0)) # min is for errors in the data (cant be hatch year twice but it happens)
    # If they were seen as a hatchling 
    if(length(hatch_yr) != 0 && (hatch_yr >= 1 & hatch_yr < k)){
      cap.data[cap.data$animal_id == i,"age"][(hatch_yr+1):k] <- 1
    } else if(length(first_adult) != 0 && (first_adult >= 1 & first_adult < k)){ #when were they seen as an adult? 
      cap.data[cap.data$animal_id == i,"age"][(first_adult+1):k] <- 1
    } #otherwise do nothing (only applies to those only seen once)
    
    # MotherID -----------------------------------------------------------------------
    mother_i <- cap.data[cap.data$animal_id == i,"mother_id"]
    
    # Drop NA and choose the unique ID
    mother_id_int <- max(mother_i[!is.na(mother_i)])
    
    # If there is a mother ID then assign it to all values for id i 
    # Otherwise just pass a zero rather than an NA
    if(length(mother_id_int) != 0){
      cap.data[cap.data$animal_id == i,"mother_id"] <- mother_id_int
    } else{
      cap.data[cap.data$animal_id == i,"mother_id"] <- 0
    }
  }
  
  # Update mate status to be zero if age is zero (for hatch year)
  cap.data <- cap.data %>% mutate(mated = ifelse(age == 0, 0, mated))
  
  # Return cleaned data
  return(cap.data)
}

assign_ids_bysex <- function(cap.data){
  
  # Assign new ids to females and males
  id_catalog <- cap.data %>% 
    select(animal_id, sex) %>% 
    arrange(sex, animal_id) %>% 
    distinct() %>%
    group_by(sex) %>% 
    summarize(jags_id = 1:n(), animal_id = animal_id) %>% 
    ungroup() 
  
  
  # Add gender specific ids (for jags/nimble modelling)
  cap.data <- cap.data %>% 
    left_join(id_catalog, by = c("animal_id", "sex")) %>% 
    left_join(rename(id_catalog,"jags_partner_id" = jags_id), by = c("partner_id" = "animal_id")) %>% 
    select(-sex.y) %>% rename(sex = sex.x) %>% 
    left_join(filter(rename(id_catalog,"jags_mother_id" = jags_id)), by = c("mother_id" = "animal_id")) %>% 
    select(-sex.y) %>% rename(sex = sex.x)
  
  return(cap.data)
  
}

# Build Jags Data for recruitment 
populate_recruit <- function(cap.data, nf, nm, k){
  
  # Recruitment matrices
  recruit_f <- matrix(NA, nrow = nf+1, ncol = k)
  recruit_m <- matrix(NA, nrow = nm+1, ncol = k)
  
  # Extract initial entry by id 
  female_init <- cap.data %>% select(sex, jags_id, initial_entry) %>% filter(sex == "F") %>% arrange(jags_id) %>% distinct()
  male_init <- cap.data %>% select(sex, jags_id, initial_entry) %>% filter(sex == "M") %>% arrange(jags_id) %>% distinct()
  
  # Populate Female 
  for(i in 1:nrow(female_init)){
    f_id <- female_init$jags_id[i]
    init_id <- female_init$initial_entry[i]
    recruit_f[f_id, init_id:k] <- 1
  }
  
  # Populate Male 
  for(i in 1:nrow(male_init)){
    m_id <- male_init$jags_id[i]
    init_m_id <- male_init$initial_entry[i]
    recruit_m[m_id, init_m_id:k] <- 1
  }
  
  # Return recruitment list
  return(list(recruit_f = recruit_f, recruit_m = recruit_m))
  
}


populate_mating <- function(cap.data, nf, nm, k){
  
  # Attempt to mate matrix
  amating_f <- matrix(NA, nrow = nf, ncol = k)
  amating_m <- matrix(NA, nrow = nm, ncol = k)
  
  # Mating status of females/males
  female_mating <- cap.data %>% filter(sex == "F") %>% select(time, jags_id, mated) %>% distinct()
  male_mating <- cap.data %>% filter(sex == "M") %>% select(time, jags_id, mated) %>% distinct()
  
  # Assign mate status to females
  for(i in 1:nrow(female_mating)){
    f_id <- female_mating$jags_id[i]
    t <- female_mating$time[i]
    amating_f[f_id,t] <- female_mating$mated[i]
  }
  
  # Assing mate status to males
  for(i in 1:nrow(male_mating)){
    m_id <- male_mating$jags_id[i]
    t <- male_mating$time[i]
    amating_m[m_id,t] <- male_mating$mated[i]
  }
  
  return(list(amating_f = amating_f, amating_m = amating_m))
}

populate_pairs <- function(cap.data, nf, nm, k){
  
  #Pairs matrix
  apairs <- array(NA, dim = c(nf+1, nm+1, k+1))
  
  # Dummy Index
  apairs[1:(nf+1),1,] <- 0
  apairs[1,1:(nm+1),] <- 0
  apairs[,,1] <- 0 
  
  
  # Females and partners
  female_mates <- cap.data %>% filter(sex == "F") %>% select(time, jags_id, jags_partner_id, jags_mother_id)
  male_mates <- cap.data %>% filter(sex == "M") %>% select(time, jags_id, jags_partner_id, jags_mother_id)
  
  # Assign pairs
  for(t in 1:k){
    for(i in 1:nf){
      female_id <- female_mates[female_mates$time == t,"jags_id"][i]
      female_partner_id <- female_mates[female_mates$time == t,"jags_partner_id"][i]
      if(is.na(female_partner_id)) next #skip if no partner
      for(j in 1:nm){
        # Zero out other pair formation
        apairs[female_id+1, , t+1] <- 0
        apairs[,female_partner_id+1,t+1] <- 0
        # Add formed pair 
        apairs[female_id+1, female_partner_id+1, t+1] <- 1
      }
    }
  }
  
  # Replace NA values with zero for mother-child pairs
  impossible_pairs <- male_mates %>% filter(!is.na(jags_mother_id)) %>% select(jags_id, jags_mother_id) %>% distinct()
  # Go through impossible pairings and zero them out across all time
  for(l in 1:nrow(impossible_pairs)){
    male_id <- impossible_pairs$jags_id[l]
    mother_id <- impossible_pairs$jags_mother_id[l]
    apairs[mother_id+1,male_id+1,] <- 0
  }
  
  # Return pairs list
  return(apairs)
}

# Build survival matrices
populate_surv <- function(cap.data, nf, nm, k){

  #Conditional Paired Survival matrices
  af <- matrix(NA, nrow = nf+1, ncol = k)
  am <- matrix(NA, nrow = nm+1, ncol = k)
  
  # Survival data.frames
  surv_f <- cap.data %>%
    select(sex, jags_id, time, surv_individual_confounded) %>%
    filter(sex == "F") %>% 
    arrange(jags_id) %>% 
    distinct()
  
  surv_m <- cap.data %>% 
    select(sex, jags_id, time, surv_individual_confounded) %>%
    filter(sex == "M") %>% 
    arrange(jags_id) %>%
    distinct()
  
  # Populate survival matrices for females
  for(i in 1:nrow(surv_f)){
    id <- surv_f$jags_id[i]
    t <- surv_f$time[i]
    af[id, t] <- surv_f$surv_individual_confounded[i]
  }
  
  # Populate survival matrices for males
  for(i in 1:nrow(surv_m)){
    id <- surv_m$jags_id[i]
    t <- surv_m$time[i]
    am[id, t] <- surv_m$surv_individual_confounded[i]
  }
  # Return surv data
  return(list(af = af, am = am))
}

# Build recapture matrices
populate_recap <- function(cap.data, nf, nm, k){
  
  # Assign memory 
  recap_f <- matrix(NA, nrow = nf+1, ncol = k)
  recap_m <- matrix(NA, nrow = nm+1, ncol = k)
  
  # Recapture data.frames
  recapture_f <- cap.data %>%
    select(sex, jags_id, time, recapture_individual) %>%
    filter(sex == "F") %>% 
    arrange(jags_id) %>% 
    distinct()
  
  recapture_m <- cap.data %>% 
    select(sex, jags_id, time, recapture_individual) %>%
    filter(sex == "M") %>% 
    arrange(jags_id) %>%
    distinct()
  
  # Populate recapture matrices for females
  for(i in 1:nrow(recapture_f)){
    id <- recapture_f$jags_id[i]
    t <- recapture_f$time[i]
    recap_f[id, t] <- recapture_f$recapture_individual[i]
  }
  
  # Populate recapture matrices for males
  for(i in 1:nrow(recapture_m)){
    id <- recapture_m$jags_id[i]
    t <- recapture_m$time[i]
    recap_m[id, t] <- recapture_m$recapture_individual[i]
  }
  
  # Return recap data
  return(list(recap_f = recap_f, recap_m = recap_m))
}


populate_apairs_f <- function(apairs,nf,nm,k){
  # Build data object and set dummy variable (represents single)
  apairs_f <- matrix(NA,nrow = nf,ncol=k+1)
  apairs_f[,1] <- nm + 1
  
  for(t in 1:k){
    for(i in 1:nf){k
      
      pair_ijt <- apairs[i+1,2:(nm+1),t+1]
      isPaired <- any(pair_ijt == 1,na.rm=T)
      
      if(isPaired){
        partner_id <- which(pair_ijt == 1)
        if(length(partner_id) > 1) stop("Bug found at time" %+% t %+% " and female" %+% i)
        apairs_f[i,t] <- partner_id
      } else{
        next
      }
      
    }
  }
  
  return(apairs_f)
}

populate_arepartner <- function(apairs_f, nf, nm, k){
  
  # Dummy matrix
  arepartner <- matrix(NA,nrow = nf, ncol = k)

  
  # Assign based on apairs_f going forward
  for(t in 1:k){
      arepartner[,t] <- 1*(apairs_f[,t] == apairs_f[,t+1])
  }
  
  # Set known initial case to 0
  arepartner[,1] <- 0
  
  
  # If previous partner is with another mate (that was observed) then arepartner must be zero
  for(time in 2:k){
    for(i in 1:nf){
      
      mask1 <- is.na(arepartner[i,time]) #is this case unknown
      mask2 <- !is.na(apairs_f[i,time]) # do we know their last partner
      mask3 <- ifelse(mask2,any(apairs_f[i,time] == c(nm+1,apairs_f[-i,time+1]), na.rm = T),FALSE) # if partner has a new mate or if past state was single
      
      # ...Then repartner is known zero
      if(mask1 & mask2 & mask3){
        arepartner[i,time] <- 0
      } else {
        next
      }
    }
  }
  
  return(arepartner)
}

populate_psi <- function(apairs, nf, nm, k){
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
  return(psi_array)
}


# Prepare data for jags
build_jags_data <- function(cap.data){
  
  # Index values
  nf <- cap.data %>% select(sex, jags_id) %>% filter(sex=="F") %>% distinct() %>% arrange(jags_id) %>% nrow()
  nm <- cap.data %>% select(sex, jags_id) %>% filter(sex=="M") %>% distinct() %>% arrange(jags_id) %>% nrow()
  k <- cap.data %>% select(time) %>% distinct() %>% nrow()
  
  # Recruitment matrices
  recruit_list <- populate_recruit(cap.data, nf, nm, k)
  recruit_f <- recruit_list[["recruit_f"]]
  recruit_m <- recruit_list[["recruit_m"]]
  
  # Attempt to mate matrix
  mating_list <- populate_mating(cap.data, nf,nm ,k)
  amating_f <- mating_list[["amating_f"]]
  amating_m <-  mating_list[["amating_m"]]
  
  # Joint Pairs Matrices
  apairs <- populate_pairs(cap.data, nf, nm, k)
  
  #Conditional Paired Survival matrices
  surv_list <- populate_surv(cap.data, nf, nm, k)
  af <- surv_list[["af"]]
  am <- surv_list[["am"]]
  
  # Conditional Paired Recapture Matrices
  recap_list <- populate_recap(cap.data, nf, nm, k)
  recap_f <- recap_list[["recap_f"]]
  recap_m <- recap_list[["recap_m"]]
  
  # Construct partner index by female
  apairs_f <- populate_apairs_f(apairs,nf,nm,k)
  
  # Construct partially repartnership matrix
  arepartner <- populate_arepartner(apairs_f, nf, nm, k)
  
  # Grab Possible Pairings indexed by f/m/time
  psi <- populate_psi(apairs, nf, nm, k)
  
  
  # Store results in list 
  jags_data <- list(nf = nf, 
                    nm = nm,
                    k = k,
                    recruit_f = recruit_f,
                    recruit_m = recruit_m,
                    amating_f = amating_f,
                    amating_m = amating_m,
                    apairs_f = apairs_f,
                    arepartner = arepartner,
                    psi = psi, 
                    af = af,
                    am = am, 
                    recap_f = recap_f,
                    recap_m = recap_m)
  
  # Return model data
  return(jags_data)
}


# NExt steps

# (DONE) Add individual covariates (weight/size/etc)
# (DONE) Add key identifiers for non-breeders (Age/Mother/Watershed Area)
# (DONE) Address previous states for survival
# (DONE) Map animal_ids to male_id and female_id ranging from 1:Nm and 1:Nf respectively 
# (DONE) Build out model matrices (joint states and stuff)

# Special mating cases
# Check for special cases in joint matrix (deal with known mates but unknown partners)
# ie. How to deal with those who are known to be mated (mother with YoY -- probably just data augment)
# Deal with NA values in the single slots (go back and check what the model does with these)
# Add Data augmentation (addresses one point above)

# Update model to account for new information
# Update data simulation to match what we see here
# Build vanilla M/F JS model with recruitment 
# Run both for many iterations and compare
# Build out post-processesing functions + figures + tables
# Run small simulation study to show that mod2 is identifiable 
# Figure out conditional states for pair density 
# Do we need a penalty function of some sort? 



# PAIR INCONSISTENCY TABLE
#      k    pc   pc2    del
# 1    4    9    8.5    0.5
# 2    7   15    14.0   1.0
# 3    8    6    5.5    0.5
# 4    9   12    13.0  -1.0
# 5   18    4    5.0   -1.0
# 6   20    4    4.5   -0.5
# 7   24    7    8.0   -1.0
