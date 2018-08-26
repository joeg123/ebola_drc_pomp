############################
## Data read in file
## Reads in data and prepares it for being added into the pomp model
############################

## Not used anymore
# ebola <- read.csv("data/Ebola_outbreak_DRC_data.csv")
# ebola$Date <- chron(as.character(ebola$Date), format=c(dates = "day mon year"))
# ebola$Cases[1] <- 0
# data <- na.omit(ebola[c("Date","Cases")])
# names(data) <- c("times","cases")
# begin <- chron("26 Jul 2014", format=c(dates = "day mon year"))
# data$times <- data$times - data$times[1]
# data$times <- as.numeric(data$times)
require(tidyverse)
require(lubridate)


drc <- read_csv("data/elife-09015-supp1-v1.csv")

drc <- drc %>% 
  mutate_at(vars(Date_of_onset_symp, Date_of_Hospitalisation, Date_of_notification, Date_hospital_discharge, Date_disease_ended, Date_of_Death), .funs = mdy) %>% 
  mutate(date_infection = pmin(Date_of_onset_symp, Date_of_Hospitalisation, Date_of_notification, Date_hospital_discharge, Date_disease_ended, Date_of_Death, na.rm=T),
         outbreak=Outbreak) %>% 
  # mutate(# date_used = if_else(!is.na(Date_of_onset_symp), "Symptom Onset", if_else(!is.na(Date_of_Hospitalisation), "Hospitalization", "Notification")),
         # date_infection = if_else(!is.na(Date_of_onset_symp), Date_of_onset_symp, if_else(!is.na(Date_of_Hospitalisation), Date_of_Hospitalisation, Date_of_notification)),
         # outbreak = Outbreak) %>% 
  filter(!is.na(date_infection)) %>% 
  dplyr::select(outbreak, date_infection) %>%
  group_by(outbreak, date_infection) %>%
  summarize(cases = n()) %>%
  #pad() %>%
  #replace_na(replace = list(cases=0)) %>%
  mutate(times = as.numeric((date_infection - (min(date_infection)))+1)) %>% 
  filter(times < 500) ## Removes the single case that is in kikwit from 1996 notification date


drc <- bind_rows(drc, read_csv("data/drc_2018_data.csv") %>% 
  mutate(event_date = ymd(event_date),
         report_date = ymd(report_date)) %>% 
  group_by(event_date, health_zone) %>% 
  filter(report_date == max(report_date, na.rm = T)) %>% 
  group_by(event_date) %>% 
  summarize(conf = sum(confirmed_cases, na.rm = T)) %>% 
  mutate(cases = c(0, diff(conf)),
         outbreak = "Equator") %>% 
  filter(cases !=0) %>% 
  mutate(times = as.numeric(event_date - min(event_date) + 1) ) %>% 
  dplyr::select(-conf) %>% 
  rename(date_infection = event_date)
) 


