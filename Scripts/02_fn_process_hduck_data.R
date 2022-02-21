# Extract Data From Excel Spreadsheets
gather_hq_data <- function(dat_dir, 
                            book1 = "Harlequin_Capture_Records_24Aug2020.xls", 
                            book2 = "HARD_Obs_20Aug94_08Sep20_revised.xls"){
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

# Clean up rolling bandlists into one duck ID
process_bandlist <- function(hq_data){
  BandList <- hq_data[["book1"]][[2]]  %>%
    rename(year_class = `Year Class`,
           first_plastic = `1st Plastic`,
           second_plastic=`2nd Plastic`,
           third_plastic = `3rd Plastic`,
           fourth_plastic = `4th Plastic`)
  return(BandList)
}

# Function for observation data if desired for use in modelling
# For now this is not being utilized however the function is still available 
process_extensive_records <- function(hq_data, year.dat){
  #Grab Bandlist Information
  BandList <- process_bandlist(hq_data)
  
  # Clean Messy Data and Format to be added to capture record list
  cap2 <- hq_data$book2[[4]]%>% 
    rowid_to_column(var = "JoinID") %>% 
    rename(CodeSex = `Sex Code`,
           PrStatus4 = `PrStaus4`,
           `YBand1` = `YBand 1`,
           `YBand2` = `YBand 2`,
           `YBand3` = `YBand 3`,
           `YBand4` = `YBand 4`,
           `YBand5` = `YBand 5`,
           `YBand6` = `YBand 6`,
           `YBand7` = `YBand 7`,
           `YBand8` = `YBand 8`,
           `YBand9` = `YBand 9`,
           `YBand10` = `YBand 10`) %>% 
    mutate(PrStatus7 = as.character(PrStatus7),
           PrStatus8 = as.character(PrStatus8),
           PrStatus9 = as.character(PrStatus9),
           PrStatus10 = as.character(PrStatus10))
  
  BandDf <- cap2 %>% 
    select(
      JoinID,
      Date,
      Year,
      CodeSex,
      Area,
      `Obs Type`,
      Band1, Band2, Band3,Band4, Band5,
      Band6, Band7, Band8,Band9,Band10
    ) %>% 
    pivot_longer(
      cols = starts_with(c("Band")),
      names_to = "NumSurvey",
      names_prefix = "Band",
      values_to = "Band"
    ) 
  
  SexDf <- cap2 %>% 
    select(
      JoinID,
      Date,
      Year,
      CodeSex,
      Area,
      `Obs Type`,
      Sex1,Sex2,Sex3,Sex4,Sex5,
      Sex6,Sex7,Sex8,Sex9,Sex10
    ) %>% 
    pivot_longer(
      cols = starts_with(c("Sex")),
      names_to = "NumSurvey",
      names_prefix = "Sex",
      values_to = "Sex"
    )  
  
  
  PrDf <- cap2 %>% 
    select(
      JoinID,
      Date,
      Year,
      CodeSex,
      Area,
      `Obs Type`,
      PrStatus1,PrStatus2,PrStatus3,PrStatus4,
      PrStatus5,PrStatus6,PrStatus7,PrStatus8,
      PrStatus9,PrStatus10,
    ) %>% 
    pivot_longer(
      cols = starts_with(c("PrStatus")),
      names_to = "NumSurvey",
      names_prefix = "PrStatus",
      values_to = "PrStatus"
    )  
  
  YBandDf <- cap2 %>% 
    select(
      JoinID,
      Date,
      Year,
      CodeSex,
      Area,
      `Obs Type`,
      `YBand1`,`YBand2`,`YBand3`,`YBand4`,`YBand5`,
      `YBand6`,`YBand7`,`YBand8`,`YBand9`,`YBand10`
    ) %>% 
    pivot_longer(
      cols = starts_with(c("YBand")),
      names_to = "NumSurvey",
      names_prefix = "YBand",
      values_to = "YBand"
    ) 
  
  
  
  cap.data2 <- BandDf %>% 
    inner_join(SexDf, by = c("JoinID", "Date", "Year", "CodeSex",  "Area","Obs Type" ,"NumSurvey")) %>% 
    inner_join(PrDf, by = c("JoinID", "Date", "Year", "CodeSex",  "Area","Obs Type" ,"NumSurvey")) %>% 
    inner_join(YBandDf, by = c("JoinID", "Date", "Year", "CodeSex",  "Area","Obs Type" ,"NumSurvey")) %>% 
    mutate(AllNA = 1*(is.na(Band) & is.na(Sex) & is.na(PrStatus) & is.na(YBand))) %>% 
    filter(AllNA == 0) %>% 
    select(-AllNA) %>% 
    left_join(na.omit(select(BandList,ID,first_plastic)),by=c("Band" = "first_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,second_plastic)),by=c("Band" = "second_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,third_plastic)),by=c("Band" = "third_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,fourth_plastic)),by=c("Band" = "fourth_plastic")) %>% 
    mutate(BANDID = pmax(ifelse(is.na(ID.x),0,ID.x),
                         ifelse(is.na(ID.y),0,ID.y),
                         ifelse(is.na(ID.x.x),0,ID.x.x),
                         ifelse(is.na(ID.y.y),0,ID.y.y))) %>% 
    select(-ID.x,
           -ID.y,
           -ID.x.x,
           -ID.y.y) %>% 
    left_join(na.omit(select(BandList,ID,first_plastic)),by=c("YBand" = "first_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,second_plastic)),by=c("YBand" = "second_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,third_plastic)),by=c("YBand" = "third_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,fourth_plastic)),by=c("YBand" = "fourth_plastic")) %>% 
    mutate(YBANDID = pmax(ifelse(is.na(ID.x),0,ID.x),
                          ifelse(is.na(ID.y),0,ID.y),
                          ifelse(is.na(ID.x.x),0,ID.x.x),
                          ifelse(is.na(ID.y.y),0,ID.y.y))) %>% 
    select(-ID.x,
           -ID.y,
           -ID.x.x,
           -ID.y.y) %>% 
    filter(!(BANDID == 0 & YBANDID == 0)) %>% 
    mutate(Mated = ifelse(PrStatus %in% c("1","2","3"),1, 
                          ifelse(PrStatus == "0",0,NA)),
           age = "AHY") %>% 
    select(-NumSurvey, -CodeSex, -PrStatus) %>% 
    rename(initial_entry = "Obs Type")
  
  # Add Mother ID
  
  cap_index <- cap.data2 %>% filter(YBANDID > 0) %>% pull(JoinID) %>% unique()
  data_list <- list()
  
  for(i in 1:length(cap_index)){
    cap_id <- cap_index[i]
    temp <- cap.data2 %>% filter(JoinID == cap_id) 
    mother_id <- temp$BANDID[1]
    temp_y <- temp %>% 
      mutate(mother_id = mother_id,
             BANDID = YBANDID,
             Mated = 0,
             Band = YBand,
             Sex = NA,
             age = "HY")
    temp_mom <- temp %>% slice(1) %>% 
      mutate(YBand  = NA,
             YBANDID = NA,
             Mated = 1,
             mother_id = 0)
    
    data_list[[i]] <- rbind(temp_mom, temp_y) %>% select(-YBand,-YBANDID)
  }
  
  yband.data <- do.call(rbind,data_list)
  
  mother_index <- yband.data %>% 
    select(BANDID, mother_id) %>% 
    distinct() %>%
    filter(mother_id > 0)
  
  cap.data3 <- cap.data2 %>% select(-YBand,-YBANDID) %>% 
    filter(!(JoinID %in% cap_index)) %>% 
    left_join(mother_index, by = "BANDID") %>% 
    rbind(yband.data) %>% 
    mutate(mother_id = ifelse(is.na(mother_id),0,mother_id))
  
  #Adhoc Fixes
  cap.data3 <- cap.data3 %>% 
    mutate(BANDID = ifelse(Band == "Yg36",483,BANDID)) %>% 
    filter(BANDID > 0)
  
  
  # Add partners
  known_pairs <- hq_data$book2[[7]] %>% 
    select(1:3) %>% 
    rename(Male = `Male Band`,
           Female = `Female Band`) %>% 
    filter(tolower(trimws(Male)) != "ub",
           tolower(trimws(Male)) != "u",
           tolower(trimws(Female)) != "ub",
           tolower(trimws(Female)) != "u") %>% 
    left_join(na.omit(select(BandList,ID,first_plastic)),by=c("Male" = "first_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,second_plastic)),by=c("Male" = "second_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,third_plastic)),by=c("Male" = "third_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,fourth_plastic)),by=c("Male" = "fourth_plastic")) %>% 
    mutate(Male = pmax(ifelse(is.na(ID.x),0,ID.x),
                       ifelse(is.na(ID.y),0,ID.y),
                       ifelse(is.na(ID.x.x),0,ID.x.x),
                       ifelse(is.na(ID.y.y),0,ID.y.y))) %>% 
    select(-ID.x,
           -ID.y,
           -ID.x.x,
           -ID.y.y) %>% 
    left_join(na.omit(select(BandList,ID,first_plastic)),by=c("Female" = "first_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,second_plastic)),by=c("Female" = "second_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,third_plastic)),by=c("Female" = "third_plastic")) %>% 
    left_join(na.omit(select(BandList,ID,fourth_plastic)),by=c("Female" = "fourth_plastic")) %>% 
    mutate(Female = pmax(ifelse(is.na(ID.x),0,ID.x),
                         ifelse(is.na(ID.y),0,ID.y),
                         ifelse(is.na(ID.x.x),0,ID.x.x),
                         ifelse(is.na(ID.y.y),0,ID.y.y))) %>% 
    select(-ID.x,
           -ID.y,
           -ID.x.x,
           -ID.y.y) %>% 
    filter(Male > 0 & Female > 0)
  
  
  known_pairs_m <- known_pairs %>% 
    mutate(Sex = "M") %>% rename(BANDID = "Male", PartnerID = "Female") %>% 
    select(Year,Sex, BANDID, PartnerID)
  known_pairs_f <- known_pairs %>% 
    mutate(Sex = "F") %>% rename(BANDID = "Female", PartnerID = "Male") %>% 
    select(Year, Sex, BANDID, PartnerID)
  
  
  pair_data <- rbind(known_pairs_m,known_pairs_f)
  
  cap.data4 <- cap.data3 %>% 
    left_join(pair_data, by = c("Year","Sex","BANDID")) %>% 
    select(-JoinID) %>% 
    distinct()
  
  # Remove Duplicate information
  cap.data5 <- cap.data4 %>% 
    group_by(Year, BANDID)  %>% 
    arrange(desc(Mated), desc(PartnerID),desc(mother_id),desc(Sex),desc(Area)) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(capture_date   = as.Date(Date,origin = "1899-12-30"),
           year = year(capture_date),
           month = month(capture_date),
           day = day(capture_date ),
           pair_status = ifelse(Mated == 1, "P", ifelse(Mated == 0, "N", "U")),
           initial_entry = ifelse(initial_entry == "C","C","RC"),
           bird_wt = NA,
           wing_crd = NA,
           tarsus_len = NA,
           culmen_len = NA) %>% 
    left_join(year.dat, by = "year") %>% 
    rename(
      "animal_id" = "BANDID",
      "sex" = "Sex",
      "mated" = "Mated",
      "partner_id" = "PartnerID"
    ) %>% 
    select(capture_date,  year, month, day, animal_id, sex,  
           partner_id, pair_status, initial_entry, age, mother_id, time,
           bird_wt, wing_crd, tarsus_len, culmen_len, Area)
  
  return(cap.data5)
}

# Format capture records into long data format with all possible recapture years in data
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
  
  # Drop hatchlings that only seen once (we only want adults)
  drop_transient_hatchlings <- CaptureRecords %>%
    mutate(age2 = ifelse(age == "HY",1,0)) %>% 
    group_by(animal_id) %>% 
    summarize(newborn = sum(age2), nobs = n()) %>% 
    ungroup() %>% 
    filter(newborn > 0 & nobs <= 1) %>% 
    pull(animal_id)
  
  # Return Processed Data
  return(CaptureRecords %>% filter(!(animal_id %in% drop_transient_hatchlings)))
}

build_cr_df <- function(hq_data){
  CaptureRecords <- process_capture_records(hq_data)
  
  
  year.dat <- select(CaptureRecords, year, time) %>% 
    distinct() %>% 
    arrange(time)
  # 
  extended.cap <- process_extensive_records(hq_data, year.dat)
# 
# 
#   # Unite Capture Records
#   CaptureRecords <- rbind(CaptureRecords, extended.cap)

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

populate_missing_mate_data <- function(cap.data){
  
  animal_ids <- sort(unique(cap.data$animal_id))
  k <- max(unique(cap.data$time))
  
  for(i in 1:length(animal_ids)){
    for(t in 1:k){
      id <- animal_ids[i]
      pid <- cap.data %>% filter(time == t, animal_id == id) %>% pull(partner_id)
      
      # Skip unmated case
      if(is.na(pid)) next
      if(!is.na(pid)) if(pid == 0) next
      
      # Otherwise add partner id to other sie
      cap.data <- cap.data %>% mutate(partner_id = ifelse(animal_id == pid & time == t, id, partner_id),
                                      mated = ifelse(animal_id == pid & time == t, 1, mated),
                                      recapture_individual  = ifelse(animal_id == pid & time == t, 1, recapture_individual),
                                      surv_individual_confounded = ifelse(animal_id == pid & time == t, 1, surv_individual_confounded))
      
    }
  }
  
  # Update initial Entry
  init_df <- cap.data %>%
    group_by(animal_id) %>% 
    summarize(initial_entry = min(which(recapture_individual==1))) 
  
  cap.data <- cap.data %>% select(-initial_entry) %>% inner_join(init_df, by = "animal_id")
  
  return(cap.data)
}



# Add implied survival states (must have been alive previously if alive now etc.)
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


add_last_capture <- function(cap.data){
  
    endpoints <- cap.data %>%
    mutate(age = ifelse(is.na(age),1,age)) %>% 
    group_by(animal_id) %>% 
    summarize(initial_entry = min(which(recapture_individual==1)),
              final_entry = max(which(recapture_individual==1)),
              min_age = min(age)) %>%
    ungroup() %>% 
    mutate(known_lifespan = final_entry-initial_entry + 1,
           max_remaining = 12 - known_lifespan,
           first_possible = ifelse(min_age == 0, initial_entry + 1, ifelse(initial_entry - max_remaining <= 1, 1, initial_entry - max_remaining)),
           last_possible = ifelse(final_entry + max_remaining >= 28, 28, final_entry + max_remaining)) 
  
  cap.data <- cap.data %>% inner_join(endpoints, by = c("animal_id", "initial_entry")) 
  
  
  return(cap.data)
}

# Final Cleanup
clean_filtered <- function(cap.data){
  
  
  
  cap.data <- cap.data %>% 
    mutate(initial_entry = ifelse(first_possible > initial_entry, first_possible, initial_entry),
           mated = ifelse(time < first_possible| time > last_possible, 0, ifelse(!is.na(age), ifelse(age== 0, 0, mated), mated)),
           partner_id = ifelse(time < first_possible| time > last_possible, 0, partner_id),
           surv_individual_confounded = ifelse(time > last_possible, 0, surv_individual_confounded),
           recapture_individual = ifelse(time < first_possible,0,recapture_individual), # this one ensures hatchlings are excluded at first capture
    )
  
  # Add minimum possible age and maximum possible age
  cap.data <- cap.data %>% 
    mutate(upper_age = ifelse(time - first_possible + 1 < 12, time - first_possible + 1, 12),
           upper_age = ifelse(upper_age < 0, 0, upper_age),
           lower_age = ifelse(time - initial_entry + 1 < 12, time - initial_entry + 1, 12),
           lower_age = ifelse(lower_age < 0, 0, lower_age))
  
  return(cap.data)
}

# Assign male id and female ids
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
  recruit_f <- matrix(NA, nrow = nf, ncol = k)
  recruit_m <- matrix(NA, nrow = nm, ncol = k)
  
  # Extract initial entry by id 
  female_init <- cap.data %>% 
    select(sex, jags_id, initial_entry, first_possible) %>% 
    filter(sex == "F") %>% 
    arrange(jags_id) %>% 
    distinct()

  male_init <- cap.data %>%
    select(sex, jags_id, initial_entry, first_possible) %>% 
    filter(sex == "M") %>%
    arrange(jags_id) %>% 
    distinct()
  
  
  # Populate Female 
  for(i in 1:nrow(female_init)){
    f_id <- female_init$jags_id[i]
    init_id <- female_init$initial_entry[i]
    recruit_f[f_id, init_id:k] <- 1
    # Add limit to when they could have entered
    first_possible <- female_init$first_possible[i]
    if(first_possible <= 1){
      next
    } else {
      recruit_f[f_id, 1:(first_possible-1)] <- 0
    }

  }
  
  # Populate Male 
  for(i in 1:nrow(male_init)){
    m_id <- male_init$jags_id[i]
    init_m_id <- male_init$initial_entry[i]
    recruit_m[m_id, init_m_id:k] <- 1
    first_possible_m <- male_init$first_possible[i]
    if(first_possible_m <= 1){
      next
    } else {
      recruit_m[m_id, 1:(first_possible_m-1)] <- 0
    }
  }
  
  # Return recruitment list
  return(list(recruit_f = recruit_f, recruit_m = recruit_m))
  
}

# Build JAGS data for desire to mate at t
populate_mating <- function(cap.data, nf, nm, k){
  
  # Attempt to mate matrix
  amating_f <- matrix(NA, nrow = nf, ncol = k)
  amating_m <- matrix(NA, nrow = nm, ncol = k)
  
  # Mating status of females/males
  female_mating <- cap.data %>%
    filter(sex == "F") %>% 
    select(time, jags_id, mated, first_possible, last_possible) %>% 
    mutate(mated = ifelse(first_possible > time| time > last_possible, 0, mated)) %>% 
    distinct()
  
  male_mating <- cap.data %>% 
    filter(sex == "M") %>% 
    select(time, jags_id, mated, first_possible, last_possible) %>% 
    mutate(mated = ifelse(first_possible > time| time > last_possible, 0, mated)) %>% 
    distinct()
  
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

# Assign known pairs and possible fates for JAGS 
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
  
  # Zero-out pairs outside of the window of possibility 
  # Mating status of females/males
  female_mating <- cap.data %>%
    filter(sex == "F") %>% 
    select(time, jags_id, mated, first_possible, last_possible) %>% 
    mutate(mated = ifelse(first_possible > time| time > last_possible, 0, mated)) %>% 
    distinct()
  
  male_mating <- cap.data %>% 
    filter(sex == "M") %>% 
    select(time, jags_id, mated, first_possible, last_possible) %>% 
    mutate(mated = ifelse(first_possible > time| time > last_possible, 0, mated)) %>% 
    distinct()
  
  # Assign no mate to females who couldn't be in population based on age
  for(i in 1:nrow(female_mating)){
    f_id <- female_mating$jags_id[i]
    t <- female_mating$time[i]
    if(!is.na(female_mating$mated[i])) if(female_mating$mated[i] == 0) apairs[f_id+1,,t+1] <- 0
  }
  
  # Assign no mate to males who couldn't be in population based on age
  for(i in 1:nrow(male_mating)){
    m_id <- male_mating$jags_id[i]
    t <- male_mating$time[i]
    if(!is.na(male_mating$mated[i])) if(male_mating$mated[i] == 0) apairs[,m_id+1,t+1] <- 0
  }
  
  # Return pairs list
  return(apairs)
}

# Build survival matrices
populate_surv <- function(cap.data, nf, nm, k){

  #Conditional Paired Survival matrices
  af <- matrix(NA, nrow = nf, ncol = k)
  am <- matrix(NA, nrow = nm, ncol = k)
  
  # Survival data.frames
  surv_f <- cap.data %>%
    select(sex, jags_id, time, surv_individual_confounded, last_possible) %>%
    #mutate(surv_individual_confounded = ifelse(time > last_possible, 0, surv_individual_confounded)) %>% 
    filter(sex == "F") %>% 
    arrange(jags_id) %>% 
    distinct()
  
  last_alive_f <- surv_f %>% 
    group_by(jags_id) %>% 
    summarize(last_alive = max(which(surv_individual_confounded == 1))) %>% 
    ungroup()
  
  surv_m <- cap.data %>% 
    select(sex, jags_id, time, surv_individual_confounded, last_possible)  %>%
   # mutate(surv_individual_confounded = ifelse(time > last_possible, 0, surv_individual_confounded)) %>% 
    filter(sex == "M") %>% 
    arrange(jags_id) %>%
    distinct()
  
  last_alive_m <- surv_m %>% 
    group_by(jags_id) %>% 
    summarize(last_alive = max(which(surv_individual_confounded == 1))) %>% 
    ungroup()
  
  # Populate survival matrices for females
  for(i in 1:nrow(surv_f)){
    id <- surv_f$jags_id[i]
    t <- surv_f$time[i]
    last_seen <- last_alive_f %>% filter(jags_id == id) %>% pull(last_alive)
    last_possible <- surv_f$last_possible[i]
    
    # Must have been alive prior to being seen (or not in pop yet)
    if(t <= last_seen){
      af[id, t] <- 1
    }
    
    # Must have perished after maximal life expectancy is reached 
    if(t > last_possible){
      af[id, t] <- 0
    }
    
  }
  
  # Populate survival matrices for males
  for(i in 1:nrow(surv_m)){
    id <- surv_m$jags_id[i]
    t <- surv_m$time[i]
    last_alive <- last_alive_m %>% filter(jags_id == id) %>% pull(last_alive)
    last_possible <- surv_m$last_possible[i]
    
    # Must have been alive prior to being seen (or not in pop yet)
    if(t <= last_alive){
      am[id, t] <- 1
    } 
    
    # Must have perished after maximal life expectancy is reached 
    if(t > last_possible){
      am[id, t] <- 0
    }
  }
  # Return surv data
  return(list(af = af, am = am))
}

# Build recapture matrices
populate_recap <- function(cap.data, nf, nm, k){
  
  # Assign memory 
  recap_f <- matrix(NA, nrow = nf, ncol = k)
  recap_m <- matrix(NA, nrow = nm, ncol = k)
  
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

# Assign compact known pairs by female by year using male ids 
populate_apairs_f <- function(apairs,amating_f,nf,nm,k){
  # Build data object and set dummy variable (represents single)
  # apairs_f <- matrix(NA,nrow = nf,ncol=k+1)
  # apairs_f[,1] <- nm + 1
  
  apairs_f <- matrix(NA,nrow = nf,ncol=k)
  
  for(t in 1:k){
    for(i in 1:nf){
      
      pair_ijt <- apairs[i+1,2:(nm+1),t+1]
      isPaired <- any(pair_ijt == 1,na.rm=T)
      
      if(isPaired){
        partner_id <- which(pair_ijt == 1)
        if(length(partner_id) > 1) stop("Bug found at time" %+% t %+% " and female" %+% i)
        apairs_f[i,t] <- partner_id
      } else{
        if(!is.na(amating_f[i,t])) if(amating_f[i,t] == 0) apairs_f[i,t] <- nm + 1
      }
      
    }
  }
  
  return(apairs_f)
}

# Assign known re-partnerships by female by year
populate_arepartner <- function(apairs_f, nf, nm, k){
  
  # Dummy matrix
  arepartner <- matrix(NA,nrow = nf, ncol = k)

  
  # Assign based on apairs_f going forward
  for(t in 2:k){
      arepartner[,t] <- 1*(apairs_f[,t] == apairs_f[,t-1]) * (apairs_f[,t-1] != nm + 1) * (apairs_f[,t] != nm + 1)
  }
  
  # Set known initial case to 0
  arepartner[,1] <- 0
  
  
  # If previous partner is with another mate (that was observed) then arepartner must be zero
  for(time in 2:k){
    for(i in 1:nf){
      
      mask1 <- is.na(arepartner[i,time]) #is this case unknown
      mask2 <- !is.na(apairs_f[i,time-1]) # do we know their last partner
      mask3 <- ifelse(mask2,any(apairs_f[i,time-1] == c(nm+1,apairs_f[-i,time]), na.rm = T),FALSE) # if partner has a new mate or if past state was single
      
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

# Assign data object with possible pairings based on known fates and known survival apriori to sampling
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
  apairs_f <- populate_apairs_f(apairs,amating_f,nf,nm,k)
  
  # Construct partially repartnership matrix
  arepartner <- populate_arepartner(apairs_f, nf, nm, k)
  
  # Grab Possible Pairings indexed by f/m/time
  psi <- populate_psi(apairs, nf, nm, k)
  
  add_dummy_row <- function(mat, x = 1){
    return(rbind(mat,rep(x,ncol(mat))))
  }
  add_dummy_col <- function(mat, x = 1){
    return(cbind(rep(x,nrow(mat)),mat))
  }
  
  # Store results in list 
  jags_data <- list(nf = nf, 
                    nm = nm,
                    k = k,
                    recruit_f = add_dummy_row(recruit_f),
                    recruit_m = add_dummy_row(recruit_m),
                    amating_f = add_dummy_row(amating_f),
                    amating_m = add_dummy_row(amating_m),
                    apairs_f = apairs_f,
                    arepartner = arepartner,
                    apairs = apairs,
                    psi = psi, 
                    af = add_dummy_row(af),
                    am = add_dummy_row(am),
                    recap_f = add_dummy_row(recap_f, 0),
                    recap_m = add_dummy_row(recap_m,0))
  
  # Return model data
  return(jags_data)
}