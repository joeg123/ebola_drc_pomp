##
# Script to calculate the time dependent R0s from each outbreak for plotting
##

library(R0)
library(tidyverse)
library(lubridate)
library(cowplot)

source("R/read_in_drc_data.R")



# Calculating the Observed R0 ---------------------------------------------

calc_td_rnot <- function(outbreak_name, df){
  ## Function to calculate the time-dependent R0 for the Ebola outbreaks
  ## df is assumed to be the drc data_frame
  ## outbreak should correspond to a single outbreak
  df <- df %>% filter(outbreak == outbreak_name)
  
  # serial interval of 15.3 days with a standard deviation of 9.3 days
  ebola_gt <- generation.time(type = "gamma", val = c(15.3, 9.3))  
  td_rnot <- est.R0.TD(epid =df$cases, t = df$date_infection, GT = ebola_gt, begin = 1, end = max(df$times), nsim = 10000)
  
  data_frame(outbreak = outbreak_name,
             date= td_rnot$epid$t, 
             day = seq_along(td_rnot$epid$t),
             cases = td_rnot$epid$incid, 
             rnot = td_rnot$R, 
             ci_lower = td_rnot$conf.int$lower, 
             ci_upper = td_rnot$conf.int$upper)
}

outbreaks <- unique(drc$outbreak)


test <- outbreaks %>% 
          map(calc_td_rnot, df = drc) %>% 
          bind_rows()



test %>% 
  filter(cases != 0, ci_upper<50) %>% 
  group_by(outbreak, week = week(date)) %>% 
  summarize(date = min(date),day =min(day), rnot = mean(rnot), ci_lower = mean(ci_lower), ci_upper = mean(ci_upper)) %>% 
  ggplot(aes(day, rnot)) + 
    facet_wrap(~outbreak, ncol = 1) +
    geom_point() +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper))



# Calculating our estimated R0 --------------------------------------------


