#Capture-Recapture Records (with nets)
wb1 <- loadWorkbook(path2dat %+% "Harlequin_Capture_Records_29Aug2012.xlsx" )
sheet_names1 <- sheets(wb1)
sheet_names1

#Read in excel tables
BandList <- read.xlsx(wb1,sheet = "Band List") %>%
  rename(first = `1st.Plastic`,
         second=`2nd.Plastic`,
         third = `3rd.Plastic`,
         fourth = `4th.Plastic`)

CaptureRecords <- read.xlsx(wb1,sheet = "Capture Records McLeod Birds") %>%  
  mutate(Capture.Date = as.Date(Capture.Date,origin = "1899-12-30"),
         year = year(Capture.Date),
         month = month(Capture.Date),
         day = mday(Capture.Date)) %>% 
  rename("animal_id" = "ID",
         "sex" = "Sex",
         "capture_date" = "Capture.Date",
         "date_of_birth" = "Year.Class",
         "age_estimate" = "Estimate",
         "pair_status" = "Pairstat",
         "partner_id" = "Mateleg",
         "initial_entry" = "Type") %>%
  arrange(animal_id,capture_date,initial_entry) %>% 
  left_join(na.omit(select(BandList,ID,first)),by=c("partner_id" = "first")) %>% 
  left_join(na.omit(select(BandList,ID,second)),by=c("partner_id" = "second")) %>% 
  left_join(na.omit(select(BandList,ID,third)),by=c("partner_id" = "third")) %>% 
  left_join(na.omit(select(BandList,ID,fourth)),by=c("partner_id" = "fourth")) %>% 
  mutate(partner_id = pmax(ifelse(is.na(ID.x),0,ID.x),
                           ifelse(is.na(ID.y),0,ID.y),
                           ifelse(is.na(ID.x.x),0,ID.x.x),
                           ifelse(is.na(ID.y.y),0,ID.y.y)
                           )
  ) %>% 
  select(-ID.x,-ID.y,-ID.x.x,-ID.y.y) %>% 
  select(capture_date,year,month,day,animal_id,sex,partner_id,pair_status,initial_entry) %>% 
  filter(!is.na(capture_date))

year.dat <- select(CaptureRecords,year) %>% 
  distinct() %>% 
  t() %>% 
  as.vector() %>%
  sort() %>% 
  data.frame() %>% 
  rename("year"=".") %>% 
  mutate(time = 1:length(year))

#Capture-Recapture format
cap.data <- merge(unique(CaptureRecords$animal_id),year.dat$time) %>% 
  rename(animal_id = x,time = y) %>% 
  left_join(year.dat,by="time") %>% 
  left_join(CaptureRecords,by=c("animal_id","year")) %>% 
  mutate(recapture_individual = ifelse(is.na(capture_date),0,1),
         mated = ifelse(!is.na(partner_id)&partner_id != 0,2,
                        ifelse(pair_status=="N"|pair_status=="YOY",1,
                               ifelse(pair_status=="U",NA,
                                      ifelse(pair_status=="P"&partner_id==0,2,NA)))))  %>% 
  arrange(animal_id,time)

#Removing animals that have multiple partners (Mate for life)
x <- cap.data %>% 
  filter(!is.na(capture_date),sex=="F",partner_id>0) %>%
  group_by(animal_id) %>% 
  summarize(num = length(unique(partner_id))) %>% 
  filter(num >1)

z <- cap.data %>% 
  inner_join(x,by="animal_id") %>%
  inner_join(cap.data,by=c("partner_id"="animal_id")) %>% 
  select(animal_id,partner_id,partner_id.y) %>% 
  distinct()

exclude <- c(z[,1],z[,2],z[,3]) %>% unique() 
exclude <- exclude[!is.na(exclude)]
exclude <- exclude[exclude>0] %>% data.frame() %>% rename("animal_id"=".")

cap.data <- cap.data %>% anti_join(exclude,by="animal_id") 

initial <- cap.data %>% filter(recapture_individual == 1) %>% 
  select(animal_id,time) %>% 
  group_by(animal_id) %>% 
  summarize(initial_entry = min(time))

gender <- filter(cap.data,!is.na(sex)) %>% 
  select(animal_id,sex) %>% 
  distinct()

cap.data <- cap.data %>% 
  select(-initial_entry,-sex) %>% 
  left_join(initial,by="animal_id") %>%
  left_join(gender,by="animal_id") %>% 
  select(-capture_date,-month,-day,-year,-pair_status) %>% 
  mutate(mated = ifelse(initial_entry > time,0,mated),
         partner_id = ifelse(initial_entry > time,0,ifelse(is.na(mated),NA,partner_id)),
         surv_individual_confounded = ifelse(recapture_individual == 0,ifelse(time<initial_entry,0,NA),1),
         recapture=ifelse(recapture_individual==1,ifelse(is.na(partner_id),ifelse(sex=="M",3,2),ifelse(partner_id >0,4,ifelse(sex=="M",3,2))),1),
         surv_confounded=ifelse(recapture==4,4,ifelse(time < initial_entry,1,ifelse(recapture==1,NA,recapture)))) %>% 
  select(animal_id,time,sex,initial_entry,partner_id,mated,recapture,surv_confounded,recapture_individual,surv_individual_confounded) %>% 
  distinct()

cap.data.list <- list()

#Assign Confounded States
for(m in unique(cap.data$animal_id)){
  
  animal <- cap.data %>% filter(animal_id == m) 
  
  if(unique(animal$sex) == "M"){
    last.seen.alive <- max(which(animal$surv_confounded %in% c(3,4)))
    
    for(i in 1:nrow(animal)){
      
      initial_entry <- animal$initial_entry[i]
      temp.recap <- animal$recapture[i]
      animal$surv_confounded[i] <- ifelse(i < last.seen.alive & is.na(animal$surv_confounded[i]),ifelse(i<initial_entry,1,3),animal$surv_confounded[i])
    }
    
  } else {
    last.seen.alive <- max(which(animal$surv_confounded %in% c(2,4)))
    
    for(i in 1:nrow(animal)){
      
      initial_entry <- animal$initial_entry[i]
      temp.recap <- animal$recapture[i]
      animal$surv_confounded[i] <- ifelse(i < last.seen.alive & is.na(animal$surv_confounded[i]),ifelse(i<initial_entry,1,2),animal$surv_confounded[i])
    }
    
  }
  
  cap.data.list[[m]] <- animal %>% 
    mutate(surv_individual_confounded = 
             ifelse(is.na(surv_individual_confounded),
                    ifelse(surv_confounded %in% c(2,3),1,NA),
                    surv_individual_confounded))
  
}


block.data <- bind_rows(cap.data.list)


#Mate for life formatting 
block.male <- block.data %>% filter(sex=="M") %>% select(-recapture,-surv_confounded) %>% 
  rename("male" = animal_id,"male_entry"=initial_entry,d_male = "mated","rc_male" = recapture_individual,"sv_male" = surv_individual_confounded) 
block.female <- block.data %>% filter(sex=="F") %>% select(-recapture,-surv_confounded) %>% 
  rename("female" = animal_id,"female_entry"=initial_entry,d_female = "mated","rc_female" = recapture_individual,"sv_female" = surv_individual_confounded) 

male.partners <- block.male %>% 
  select(male,partner_id) %>% 
  mutate(partner_id =ifelse(is.na(partner_id),0,partner_id)) %>% 
  group_by(male) %>%
  summarize(partner=max(partner_id))

female.partners <- block.female %>% 
  select(female,partner_id) %>% 
  mutate(partner_id =ifelse(is.na(partner_id),0,partner_id)) %>% 
  group_by(female) %>%
  summarize(partner=max(partner_id))

groups1 <- male.partners %>% 
  inner_join(block.female,by=c("partner"="female")) %>%
  select(male,partner) %>% 
  distinct() %>% 
  rename(female = partner)
  
groups2 <- female.partners %>% 
  inner_join(block.male,by=c("partner"="male")) %>%
  select(female,partner) %>% 
  distinct() %>% 
  rename(male = partner) %>% 
  select(male,female)

#Investigate discrepancies
female.partners %>% filter(partner !=0) %>% distinct() %>% anti_join(groups1)
male.partners %>% filter(partner !=0) %>% distinct() %>% anti_join(groups2)

#Set 173,101, and 174 as singles.
#61 and 163 are a pair bond but 61 seemed to be tagged with an incorrect partner_id under 163s record

groups <- full_join(groups1,groups2) %>% filter(male != 173)
single.females <- female.partners %>% filter(partner==0|female %in% c(174,101)) %>% mutate(male=0) %>% select(male,female)
single.males <- male.partners %>% filter((partner==0|male==173)&male != 163) %>% mutate(female=0) %>% select(male,female)

hd_entities <- bind_rows(groups,single.females,single.males) %>% 
  full_join(block.male,by="male") %>% select(-sex,-partner_id) %>% 
  full_join(block.female,by=c("female","time")) %>% select(-sex,-partner_id) %>% 
  distinct() %>% filter(!is.na(time)) %>% 
  mutate(male=ifelse(is.na(male),0,male),
         female=ifelse(is.na(female),0,female))


hd_entities_m <- hd_entities %>% filter(male!=0,female==0) %>%
  mutate(initial_entry=male_entry,d=1,rc=ifelse(rc_male==1,3,1)) %>% 
  select(male,female,time,initial_entry,d,rc)

hd_entities_f <- hd_entities %>% 
  filter(male==0,female!=0) %>% 
  mutate(initial_entry=female_entry,d=1,rc=ifelse(rc_female==1,2,1))  %>% 
  select(male,female,time,initial_entry,d,rc)

hd_entities_g <- hd_entities %>%
  filter(male!=0,female!=0) %>% 
  mutate(initial_entry = pmax(male_entry,female_entry),
         rc=ifelse(rc_female==1&rc_male==1,4,ifelse(rc_female==1,2,ifelse(rc_male==1,3,1))),
         d=ifelse((time<initial_entry)|(is.na(d_male)&is.na(d_female))|rc!=4,NA,
                  ifelse(is.na(d_male),d_female,ifelse(is.na(d_female),d_male,pmax(d_male,d_female))))) %>% 
  select(male,female,time,initial_entry,rc,d) %>% distinct()

entities <- bind_rows(hd_entities_g,hd_entities_m,hd_entities_f) %>% 
  mutate(a=ifelse(time==initial_entry,4,NA))

index.set <- entities %>%
  select(male,female) %>%
  distinct() %>%
  mutate(ID = 1:n()) 

entities <- entities %>% 
  inner_join(index.set,by=c("male","female")) %>% 
  mutate(a = ifelse(time<initial_entry|rc==1,NA,rc)) %>% 
  select(ID,time,male,female,initial_entry,rc,a,d) %>% 
  distinct()

k <- max(entities$time)


#Assign values sequentially
for(i in unique(entities$ID)){
  init <- filter(entities,ID==i) %>% select(initial_entry) %>% unique() %>% t() %>% as.vector()
  #Confounded survival states at time i
  state_confounded <- entities %>% filter(ID == i) %>% select(a) %>% t() %>% as.vector()
  #Positions 
  last.female.seen <- ifelse(abs(max(which(state_confounded %in% c(2))))==Inf,0,max(which(state_confounded %in% c(2))))
  last.male.seen <- ifelse(abs(max(which(state_confounded %in% c(3))))==Inf,0,max(which(state_confounded %in% c(3))))
  last.both.seen <- ifelse(abs(max(which(state_confounded %in% c(4))))==Inf,0,max(which(state_confounded %in% c(4))))
  last.seen <- max(last.both.seen,last.male.seen,last.female.seen)
  #Update Confounded survival Information
  for(l in (init):last.seen){
    state_confounded[l] <- ifelse(last.both.seen>=l|(last.male.seen>=l&last.female.seen>=l),4,ifelse(last.male.seen>=l,3,ifelse(last.female.seen>=l,2,NA)))
    entities[(i - 1) * k + l,7] <- state_confounded[l]
  }
}

###################################
#Functions to format data for Jags#
###################################

#Format MRC data into dgr model 
cjs_dgr_data <- function(Data){
  d <- Data[,c("time","d")] %>% as.tibble()
  a <- Data[,c("time","a")] %>% as.tibble()
  X <- Data[,c("time","rc")] %>% as.tibble()
  first <- Data %>% select(ID,initial_entry) %>% distinct() %>% select(initial_entry) %>% t() %>% as.vector()
  M <- Data %>% select(ID) %>% distinct() %>% nrow()
  K <- Data %>% select(time) %>% distinct() %>% nrow()
  return(list("X" = X,"a" = a,"d"=d,"M" = M, "K" = K,"first"=first))
}

#Generate DGR Initial Conditions for n amount of parallel chains
cjs_dgr_init <- function(k,chains=1){
  initial_states <- list()
  for(i in 1:chains){
    initial_states[[i]] <- list("PhiF" = runif(k-1,0,1),
                                "PhiM" = runif(k-1,0,1),
                                "PF" = runif(k-1,0,1),
                                "PM" = runif(k-1,0,1),
                                "delta" = runif(k-1,0,1))
  }
  return(initial_states)
}


#Format MRC data into rho model 

cjs_rho_data <- function(Data){
  a <- Data[,c("time","a")] %>% as.tibble()
  X <- Data[,c("time","rc")] %>% as.tibble()
  first <- Data %>% 
    select(ID,initial_entry) %>% distinct() %>% select(initial_entry) %>% t() %>% as.vector()
  M <- Data %>% select(ID) %>% distinct() %>% nrow()
  K <- Data %>% select(time) %>% distinct() %>% nrow()
  return(list("X" = X,"a" = a,"M" = M, "K" = K,"first"=first))
}

#Generate Rho Initial Conditions for n amount of parallel chains
cjs_rho_init <- function(k,chains=1){
  initial_states <- list()
  for(i in 1:chains){
    initial_states[[i]] <- list("PhiF" = runif(k-1,0,1),
                                "PhiM" = runif(k-1,0,1),
                                "PF" = runif(k-1,0,1),
                                "PM" = runif(k-1,0,1))
  }
  return(initial_states)
}

saveRDS(entities,path2out %+% "HD-Data.rds")

