##
# Script to calculate the time dependent R0s from each outbreak for plotting
##

library(R0)
library(tidyverse)
library(lubridate)
library(cowplot)
library(magrittr)
source("R/read_in_drc_data.R")

if(!dir.exists("ms_figs")){
  dir.create("ms_figs")
}

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


obs_rnots <- outbreaks %>% 
          map(calc_td_rnot, df = drc) %>% 
          bind_rows()

save(obs_rnots, file="data_produced/fig_results/obs_rnots.rda")

# Calculating our estimated R0 - superspreading --------------------------------------------

calc_ss_rnot <- function(beta, p, times, gamma = 1/7.411374){
  rep(beta*p/gamma, length(times))
}

calc_outbreak_ss_rnot <- function(outbrk, parms_df, outbreak_df){
  o_df <- outbreak_df %>% filter(outbreak == outbrk)
  p_df <- parms_df %>% filter(outbreak == outbrk)
  data_frame(outbreak = outbrk,
             date_infection = o_df$date_infection,
             times = o_df$times,
             rnot = calc_ss_rnot(p_df$beta_estimate, p_df$p_estimate, o_df$times))
}

load("data_produced/outbreak_rda/parm_est_ss.rda")
outbreaks <- unique(drc$outbreak)

## First setup output for use in functions
ss_parms <- ss_output %>% 
  gather(key, value, estimate:upper) %>% 
  unite(key, parameter, key) %>% 
  spread(key, value)

## calculate the superspreading R0 estimate for each time point
ss_rnots <- outbreaks %>% map(calc_outbreak_ss_rnot, parms_df = ss_parms, outbreak_df=drc) %>% 
  bind_rows()  %>% 
  mutate(week_times = times/7) 

# save(ss_rnots, file = "data_produced/fig_results/ss_rnots.rda")


# Calculate R0s for intervention model ------------------------------------

calc_int_rnot <- function(beta, k, tau, times, gamma = 1/7.411374){
  betas <- if_else(times < tau, beta, (beta)*(exp(-k*(times-tau))))
  betas/gamma
}

calc_outbreak_int_rnot <- function(outbrk, parms_df, outbreak_df){
  o_df <- outbreak_df %>% filter(outbreak == outbrk)
  p_df <- parms_df %>% filter(outbreak == outbrk)
  data_frame(outbreak = outbrk,
             date_infection = o_df$date_infection,
             times = o_df$times,
             rnot = calc_int_rnot(p_df$beta_estimate, p_df$k_estimate, p_df$tau_estimate, o_df$times))
}

load("data_produced/outbreak_rda/parm_est_int.rda")
int_parms <- int_output %>% gather(key, value, estimate:upper) %>% 
  unite(key, parameter, key) %>% 
  spread(key, value)

int_rnots <- outbreaks %>% map(calc_outbreak_int_rnot, parms_df = int_parms, outbreak_df=drc) %>% 
  bind_rows() %>% 
  mutate(week_times = times/7) 

# save(int_rnots, file = "data_produced/fig_results/int_rnots.rda")


ss_rnots <- ss_rnots %>% mutate(model = "ss")
int_rnots <- int_rnots %>% mutate(model = "int")
mod_rnots <- bind_rows(ss_rnots, int_rnots)
save(mod_rnots, file = "data_produced/fig_results/mod_rnots.rda")
# Old code for bounds, tabled for now -------------------------------------
# calc_ss_rnot_bounds <- function(beta, p, gamma = 1/7.411374){
#   df = data.frame(matrix(quantile(rgeom(1000000, prob = gamma / (beta+gamma)) * rbernoulli(1000000, p = p), 
#                                   probs = c(0.025, 0.975)), nrow = 1))
#   names(df) <- c("rnot_lwr", "rnot_upr")
#   df
# }
# 
# calc_ss_rnot_bounds_stat <- function(interval, beta_estimate, beta_lower, beta_upper, p_estimate, p_lower, p_upper, gamma = 1/7.411374){
#   var_p <- ((p_upper - p_lower) / 3.92)^2
#   var_beta <- ((beta_upper - beta_lower) / 3.92)^2
#   
#   ## Equation assumes independence between variables...
#   var_rnot <- beta_estimate^2 * var_p + p_estimate^2 * var_beta 
#   
#   med_rnot <- beta_estimate * p_estimate / gamma
#   
#   if(interval == "lower"){
#     med_rnot - 1.96 * sqrt(var_rnot)
#   } else{
#     med_rnot + 1.96 * sqrt(var_rnot)
#   }
# }
# ss_parms <- ss_parms %>% mutate(med_rnot = calc_ss_rnot(beta_estimate, p_estimate)) %>% 
#   mutate(result = map2(beta_estimate, p_estimate, calc_ss_rnot_bounds)) %>% 
#   unnest() %>% 
#   mutate(rnot_lwr_stat = calc_ss_rnot_bounds_stat("lower", beta_estimate, beta_lower, beta_upper, p_estimate, p_lower, p_upper),
#          rnot_upr_stat = calc_ss_rnot_bounds_stat("upper", beta_estimate, beta_lower, beta_upper, p_estimate, p_lower, p_upper),
#          outbreak = factor(outbreak, levels = outbreaks_by_type))
