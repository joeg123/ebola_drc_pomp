## Investigating 2018 ebola drc outbreak data
library(tidyverse)
library(lubridate)

drc <- read_csv("data/drc_2018_data.csv")


drc %>% 
  
  # group_by(date) %>% 
  # summarize(conf = sum(confirmed_cases), prob = sum(probabl_cases))
  group_by(health_zone) %>% 
  mutate(new_col = new_confirmed + new_probable) %>% 
  filter(!is.na(new_col)) %>% 
  mutate(cum_new = cumsum(new_col)) %>% 
  ggplot(aes(date, cum_new, color = health_zone)) + geom_line()


sum_drc_df <- drc %>% 
  mutate(event_date = ymd(event_date),
         report_date = ymd(report_date)) %>% 
  group_by(event_date, health_zone) %>% 
  filter(report_date == max(report_date, na.rm = T)) %>% 
  group_by(event_date) %>% 
  summarize(conf = sum(confirmed_cases, na.rm = T),
            prob = sum(probable_cases, na.rm = T))

## Focus on the confirmed cases
sum_drc_df %>%   
  ggplot(aes(event_date, conf)) + geom_line()

sum_drc_df %>%   
  mutate(inc = c(0, diff(sum_drc_df$conf))) %>% 
  ggplot(aes(event_date, inc)) + geom_bar(stat="identity")

  
