#####
## Manuscript Figure Creator
#####

library(tidyverse)
library(cowplot)
library(lubridate)
library(magrittr)
sapply(c("R/read_in_drc_data.R"), source)



# Setup objects used in many figures --------------------------------------
outbreaks_by_date <- c("Yambuku", "Kikwit", "Mweka2007","Mweka2008", "Isiro", "Boende")
outbreaks_by_size <- c("Mweka2008", "Isiro", "Boende", "Mweka2007", "Kikwit", "Yambuku" )

outbreaks_by_type <- c("Yambuku", "Mweka2008", "Boende", "Kikwit", "Mweka2007", "Isiro")




# Figure 1 ----------------------------------------------------------------
# Description of problem and how we did solved it. 


# Figure 2 ----------------------------------------------------------------
# Parameter estimate figure

## Read in the parameter estimates
load("data_produced/outbreak_rda/parm_est_int.rda")
load("data_produced/outbreak_rda/parm_est_ss.rda")

ss_output %>% 
  group_by(parameter) %>% 
  mutate(outbreak_means = mean(estimate)) %>% 
  mutate(outbreak = factor(outbreak, levels = outbreaks_by_type)) %>% 
  ggplot(aes(estimate, outbreak)) + 
  geom_point() +
  facet_wrap(~parameter, scales = "free_x", ncol=1) +
  geom_errorbarh(aes(xmin = lower, xmax=upper)) +
  background_grid(major = "x") +
  geom_vline(aes(xintercept = outbreak_means), lty=2, col = "red") +
  labs(y = "", x = "Parameter Estimate (95% CI)")

  

int_output %>% 
  group_by(parameter) %>% 
  mutate(outbreak_means = mean(estimate)) %>% 
  mutate(outbreak = factor(outbreak, levels = outbreaks_by_type)) %>% 
  ggplot(aes(estimate, outbreak)) + 
  geom_point() +
  facet_wrap(~parameter, scales = "free_x", ncol=1) +
  geom_errorbarh(aes(xmin = lower, xmax=upper)) +
  background_grid(major = "x") +
  geom_vline(aes(xintercept = outbreak_means), lty=2, col = "red") +
  labs(y = "", x = "Parameter Estimate (95% CI)")
# theme(axis.text.x = element_text(angle = 45))




# Rnot comparison figure --------------------------------------------------
load("data_produced/fig_results/obs_rnots.rda")
load("data_produced/fig_results/mod_rnots.rda")

obs_rnots <- obs_rnots %>% 
  filter(cases != 0, ci_upper<50) %>%
  # group_by(outbreak, week = week(date)) %>%
  group_by(outbreak, week = floor_date(date, unit = "week")) %>%
  summarize(date = min(week), rnot = mean(rnot), ci_lower = mean(ci_lower), ci_upper = mean(ci_upper)) %>%
  mutate(week = (week - min(week))/7)
  

mod_rnots %>% 
  ggplot(aes(week_times, rnot)) +
  geom_line() +
  facet_grid(outbreak ~ model) +
  geom_point(data=obs_rnots, aes(week, rnot), inherit.aes=FALSE) +
  geom_errorbar(data = obs_rnots, aes(x = week, ymin = ci_lower, ymax = ci_upper), inherit.aes = FALSE)

save_plot("ms_figs/ss_outbrk_rnot_plot.pdf", ss_outbrk_rnot_plot, base_height = 10, base_aspect_ratio = .5)


## Figure X: Error Bars x = parms, color = outbreak will be two separate figures for the two different models
# facet by parameter, x = outbreak, may not need to color it, but may need to for clearity

## Figure X: Calculate Ro over time. ss is const., int changes with function e^kxt 
# compare distribution of Ro over time, raw data, int, and ss 
# plot digitizer






