#####
## Manuscript Figure Creator
#####

library(tidyverse)
library(padr)
library(dplyr)
library(knitr)
library(magrittr)
library(ggplot2)
library(xtable)
library(grid)
library(gridExtra)

results_ss <- read_csv("SS_results.csv")

results_int <- read_csv("int_results.csv")

## Figure 1:
# Description of problem and how we did solved it. 

## Figure X: Parameter Estimates
# SS vs Intervention (panel for ss and panel for int)
# 6 outbreaks
# 2 or 3 parameters (we dont want to directly compare the different parameters)

pd <- position_dodge(0.1) # move them .05 to the left and right

results %>% ggplot(aes(x=X1, y=p)) + 
  geom_errorbar(aes(ymin=p_lower, ymax=p_upper), width=.1, position=pd) +
  geom_point(position=pd) + labs(xlab("Outbreak"), ylab("Proportion of Spreaders")) -> p_plot

results %>% ggplot(aes(x=X1, y=beta)) +
  labs(xlab("Outbreak"), ylab("Rate of Infection")) +
  geom_errorbar(aes(ymin=beta_lower, ymax=beta_upper), width=.1, position=pd) +
  geom_point(position=pd) -> beta_plot

grid.arrange(p_plot, beta_plot, top="Super-Spreading Model Estimated Parameters", ncol=2)


## Intervention Model

results_int %>% ggplot(aes(x=X1, y=k)) + 
  geom_errorbar(aes(ymin=k_lower, ymax=k_upper), width=.1, position=pd) +
  geom_point(position=pd) + labs(xlab("Outbreak"), ylab("Rate of Decay in Beta")) -> k_plot

results_int %>% ggplot(aes(x=X1, y=beta)) +
  labs(xlab("Outbreak"), ylab("Rate of Infection")) +
  geom_errorbar(aes(ymin=beta_lower, ymax=beta_upper), width=.1, position=pd) +
  geom_point(position=pd) -> beta_plot

results_int %>% ggplot(aes(x=X1, y=tau)) +
  labs(xlab("Outbreak"), ylab("Time Until Intervention")) +
  geom_errorbar(aes(ymin=tau_lower, ymax=tau_upper), width=.1, position=pd) +
  geom_point(position=pd) -> tau_plot

grid.arrange(k_plot, beta_plot, tau_plot, top="Intervention Model Estimated Parameters", ncol=2)


## Figure X: Error Bars x = parms, color = outbreak will be two separate figures for the two different models
# facet by parameter, x = outbreak, may not need to color it, but may need to for clearity

## Figure X: Calculate Ro over time. ss is const., int changes with function e^kxt 
# compare distribution of Ro over time, raw data, int, and ss 
# plot digitizer
  

