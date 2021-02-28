
gather_hq_data <- function(dat_dir, 
                            book1 = "Harlequin_Capture_Records_29Aug2012.xls", 
                            book2 = "HARD_Obs_10May95_31Aug2012.xls"){
  #Storage
  dat1 <- list()
  dat2 <- list()
  
  # 8 sheets in book1
  for(i in 1:8){
    dat1[[i]] <-  read_excel(dat_dir %+% book1, sheet = i)
  }
 
  # 5 sheets in book2
  for(j in 1:5){
    dat2[[j]] <- read_excel(dat_dir %+% book2, sheet = j)
  }
  
  # Return two data lists
  return(list(book1 = dat1, book2 = dat2))
}

hq_data <- gather_hq_data(dat_dir)


process_bandlist <- function(hq_data){
  BandList <- hq_data[["book1"]][[1]]  %>%
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
  CaptureRecords <- hq_data[["book1"]][[2]] %>%  
    mutate(Capture.Date = as.Date(`Capture Date`,origin = "1899-12-30"),
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
           "first_plastic" = `1stPlastic`,
           "age" = "Age",
           "year_class" = `Year Class`,
           "mother_id"  = MotherID) %>%
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
           pair_status,initial_entry, age, mother_id, time, Weight,             
           Wing,Tarsus, Culmen, Area) %>% 
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
           Weight, Wing, Tarsus, Culmen,
           recapture_individual,surv_individual_confounded) %>% 
    distinct()
  
  # Drop duplicate captures for females
  cap.data <- cap.data %>% group_by(animal_id, time) %>% slice(1) %>% ungroup() %>% as.data.frame()
  
  return(cap.data)
} 

cap.data <- build_cr_df(hq_data) %>%  add_implied_states() %>% assign_ids_bysex()

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


# NExt steps

# (DONE) Add individual covariates (weight/size/etc)
# (DONE) Add key identifiers for non-breeders (Age/Mother/Watershed Area)
# (DONE) Address previous states for survival
# Map animal_ids to male_id and female_id ranging from 1:Nm and 1:Nf respectively 
# Build out model matrices (joint states and stuff)
# Add Data augmentation 
# How to deal with those who are known to be mated (mother with YoY)
# Update model to account for new information
# Update data simulation to match what we see here
# Build vanilla M/F JS model with recruitment 
# Run both for many iterations and compare
# Build out post-processesing functions + figures + tables
# Run small simulation study to show that mod2 is identifiable 
# Figure out conditional states for pair density 
# Do we need a penalty function of some sort? 
