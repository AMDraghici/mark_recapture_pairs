
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

hq_data <- gather_hq_data(path2dat)


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
  
  
  BandList <- process_bandlist(hq_data)
  
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
           first_plastic = `1stPlastic`,
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
    )
    ) %>% 
    select(-ID.x,-ID.y,-ID.x.x,-ID.y.y) %>% 
    select(capture_date,year,month,day,animal_id,sex,partner_id,pair_status,initial_entry, age, mother_id) %>% 
    filter(!is.na(capture_date))
  return(CaptureRecords)
}


