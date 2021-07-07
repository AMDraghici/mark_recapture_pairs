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
    Year,
    CodeSex,
    Area,
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
    Year,
    CodeSex,
    Area,
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
    Year,
    CodeSex,
    Area,
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
    Year,
    CodeSex,
    Area,
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
  inner_join(SexDf, by = c("JoinID", "Year", "CodeSex",  "Area", "NumSurvey")) %>% 
  inner_join(PrDf, by = c("JoinID", "Year", "CodeSex", "Area", "NumSurvey")) %>% 
  inner_join(YBandDf, by = c("JoinID", "Year", "CodeSex", "Area","NumSurvey")) %>% 
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
                        ifelse(PrStatus == "0",0,NA))) %>% 
  select(-NumSurvey, -CodeSex, -PrStatus) 

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
           Sex = NA)
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
  rename(
    "year" = "Year",
    "animal_id" = "BANDID",
    "sex" = "Sex",
    "mated" = "Mated",
    "partner_id" = "PartnerID"
  ) %>% 
  select(year,animal_id,sex,partner_id,mated,mother_id, Area)

# Map to cap.data 
# Make sure JAGS IDs get added after
# Make inclusion of this data optional
# Check special cases (like known deaths and stuff)
# WHen adding update the recapture and surival columns
