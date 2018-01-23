############################
## Data read in file
## Reads in data and prepares it for being added into the pomp model
############################

ebola <- read.csv("data/Ebola_outbreak_DRC_data.csv")
ebola$Date <- chron(as.character(ebola$Date), format=c(dates = "day mon year"))
ebola$Cases[1] <- 0
data <- na.omit(ebola[c("Date","Cases")])
names(data) <- c("times","cases")
begin <- chron("26 Jul 2014", format=c(dates = "day mon year"))
data$times <- data$times - data$times[1]
data$times <- as.numeric(data$times)



drc <- read_csv("data/elife-09015-supp1-v1.csv")

drc <- drc %>% 
  mutate(date_infection = if_else(!is.na(Date_of_onset_symp), Date_of_onset_symp, if_else(!is.na(Date_of_Hospitalisation), Date_of_Hospitalisation, Date_of_notification)),
         outbreak = Outbreak) %>% 
  filter(!is.na(date_infection)) %>% 
  select(outbreak, date_infection) %>%
  mutate(date_infection = mdy(date_infection)) %>% 
  group_by(outbreak, date_infection) %>%
  summarize(cases = n()) %>%
  #pad() %>%
  #replace_na(replace = list(cases=0)) %>%
  mutate(times = as.numeric((date_infection - (min(date_infection)))+1))
  
  



# drc %>%
#   ggplot(aes(times, cases)) + facet_wrap(~outbreak, scales="free_y", nrow=6) + geom_bar(stat="identity")
