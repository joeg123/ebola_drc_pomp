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
short_outbreaks <- c("Yambuku", "Mweka2008", "Boende")



# Figure 1 ----------------------------------------------------------------
# Description of problem and how we did solved it. 
drc %>% 
  ungroup() %>% 
  mutate(outbreak = factor(outbreak, levels = outbreaks_by_type),
         outbreak_type = if_else(outbreak %in% short_outbreaks, "short", "long")) %>% 
  ggplot(aes(times, cases)) + 
    geom_bar(stat="identity",fill="black", color="black") + 
    facet_wrap(~outbreak, ncol = 2, dir = "v") +
    labs(y = "Cases", x = "Day of Outbreak") +
    panel_border(colour = "black") +
    theme(strip.background = element_rect(fill = NA) ) -> cases_init_plot

save_plot("ms_figs/cases_init_plot.pdf", cases_init_plot, base_height = 5, base_aspect_ratio = 1.8)




# Figure 2 ----------------------------------------------------------------
# Parameter estimate figure

## Read in the parameter estimates
load("data_produced/outbreak_rda/parm_est_int.rda")
load("data_produced/outbreak_rda/parm_est_ss.rda")

ss_output %>% 
  group_by(parameter) %>% 
  mutate(outbreak_means = mean(estimate)) %>% 
  mutate(outbreak = factor(outbreak, levels = rev(outbreaks_by_type))) %>% 
  ungroup() %>% 
  mutate(parameter = if_else(parameter == "beta", "beta", "rho")) %>% 
  ggplot(aes(estimate, outbreak)) + 
    geom_point() +
    facet_wrap(~parameter, scales = "free_x", ncol=2, labeller = label_parsed) +
    geom_errorbarh(aes(xmin = lower, xmax=upper)) +
    background_grid(major = "x") +
    geom_vline(aes(xintercept = outbreak_means), lty=2, col = "red") +
    panel_border(colour = "black", size=0.7) + 
    theme(strip.text = element_text(size=20), strip.background = element_rect(fill=NA)) +
    labs(y = "", x = "Parameter Estimate (95% CI)") -> ss_parm_plot


int_output %>% 
  group_by(parameter) %>% 
  mutate(outbreak_means = mean(estimate)) %>% 
  mutate(outbreak = factor(outbreak, levels = rev(outbreaks_by_type))) %>% 
  ungroup() %>% 
  mutate(parameter = if_else(parameter == "beta", "beta[0]", parameter)) %>% 
  ggplot(aes(estimate, outbreak)) + 
    geom_point() +
    facet_wrap(~parameter, scales = "free_x", ncol=3, labeller = label_parsed) +
    geom_errorbarh(aes(xmin = lower, xmax=upper)) +
    background_grid(major = "x") +
    geom_vline(aes(xintercept = outbreak_means), lty=2, col = "red") +
    panel_border(colour = "black", size=0.7) +
    theme(strip.text = element_text(size=20), strip.background = element_rect(fill=NA)) +
    labs(y = "", x = "Parameter Estimate (95% CI)") -> int_parm_plot


parm_estimates_plot <- plot_grid(ss_parm_plot, int_parm_plot, nrow = 2, labels = "AUTO")

save_plot("ms_figs/parm_estimates_plot.pdf", parm_estimates_plot, base_height = 8, base_aspect_ratio = 1)




# Figure 3 ----------------------------------------------------------------
# Rnot comparison figure 
load("data_produced/fig_results/obs_rnots.rda")
load("data_produced/fig_results/mod_rnots.rda")

obs_rnots <- obs_rnots %>% 
  filter(cases != 0, ci_upper<50) %>%
  # group_by(outbreak, week = week(date)) %>%
  group_by(outbreak, week = floor_date(date, unit = "week")) %>%
  summarize(date = min(week), rnot = mean(rnot), ci_lower = mean(ci_lower), ci_upper = mean(ci_upper)) %>%
  mutate(week = (week - min(week))/7) %>% 
  ungroup() %>% 
  mutate(outbreak = factor(outbreak, levels = outbreaks_by_type))
  
mod_labeller <- function(variable,value){
  mod_names <- list(
    'int'="Behavior Change",
    'ss'="Superspreading"
  )
  return(mod_names[value])
}

mod_names <- c(
  'int'="Behavior Change",
  'ss'="Superspreading"
)

mod_rnots %>% 
  group_by(outbreak, model) %>% 
  summarize(lwr = min(rnot_lwr), upr =max(rnot_upr)) %>% 
  ggplot(aes(model, upr)) + geom_point()
  


mod_rnots %>% 
  mutate(outbreak = factor(outbreak, levels = outbreaks_by_type)) %>% 
  ggplot(aes(week_times, rnot)) +
    facet_grid(outbreak ~ model, labeller = labeller(model = mod_names), scales = "free_y") +
    geom_point(data=obs_rnots, aes(week, rnot), color = "grey", inherit.aes = FALSE) +
    geom_errorbar(data = obs_rnots, aes(x = week, ymin = ci_lower, ymax = ci_upper), color = "grey", inherit.aes = FALSE) +
    geom_line(color = "black", lty=2) +
    geom_ribbon(aes(ymin = rnot_lwr, ymax = rnot_upr), alpha=.2)+
    panel_border(colour = "black") +
    theme(strip.background = element_rect(fill=NA))

save_plot("ms_figs/ss_outbrk_rnot_plot.pdf", ss_outbrk_rnot_plot, base_height = 10, base_aspect_ratio = .5)


## Figure X: Error Bars x = parms, color = outbreak will be two separate figures for the two different models
# facet by parameter, x = outbreak, may not need to color it, but may need to for clearity

## Figure X: Calculate Ro over time. ss is const., int changes with function e^kxt 
# compare distribution of Ro over time, raw data, int, and ss 
# plot digitizer


# load("data_produced/fig_results/obs_rnots.rda")
# obs_rnots %>% 
#   filter(cases != 0, ci_upper<50) %>%
#   # group_by(outbreak, week = week(date)) %>%
#   group_by(outbreak, week = floor_date(date, unit = "week")) %>%
#   summarize(date = min(week), rnot = mean(rnot), ci_lower = mean(ci_lower), ci_upper = mean(ci_upper)) %>%
#   mutate(week = (week - min(week))/7) %>% 
#   ungroup() %>% 
#   mutate(outbreak = factor(outbreak, levels = outbreaks_by_type)) %>% 
#   ggplot(aes(week*7, rnot)) + 
#     facet_grid(outbreak~.) +
#     geom_point() +
#     geom_errorbar(aes(x = week*7, ymin = ci_lower, ymax = ci_upper)) +
#     labs(y = expression("R"[0]), x = "Day of Outbreak") +
#     theme(strip.background = element_rect(fill = NA)) -> rnot_init_plot
# 
# 
# fig1_cases_rnot <- plot_grid(cases_init_plot, rnot_init_plot, align="v", labels = "AUTO")
# 
# save_plot("ms_figs/fig1_cases_rnot.pdf", fig1_cases_rnot, base_height = 8, base_aspect_ratio = 1)




