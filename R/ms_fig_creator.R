#####
## Manuscript Figure Creator
#####
rm(list = ls())


library(tidyverse)
library(padr)
library(chron)
library(dplyr)
library(knitr)
library(magrittr)
library(ggplot2)
library(xtable)
library(grid)
library(gridExtra)

sapply(c("R/helper_functions.R"), source)

outbrk_list <- c("Yambuku", "Kikwit", "Mweka2007","Mweka2008", "Isiro", "Boende")

results_ss <- read_csv("SS_results.csv")

results_int <- read_csv("int_results.csv")

## Figure 1:
# Description of problem and how we did solved it. 

## Figure X: Parameter Estimates
# SS vs Intervention (panel for ss and panel for int)
# 6 outbreaks
# 2 or 3 parameters (we dont want to directly compare the different parameters)

pd <- position_dodge(0.1) # move them .05 to the left and right

results_ss %>% ggplot(aes(x=X1, y=p)) + 
  geom_errorbar(aes(ymin=p_lower, ymax=p_upper), width=.1, position=pd) +
  geom_point(position=pd) + labs(xlab("Outbreak"), ylab("Proportion of Spreaders")) -> p_plot

results_ss %>% ggplot(aes(x=X1, y=beta)) +
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

int_beta_calc <- function(outbreak, tau1, k, beta0) {
  df <- data.frame(matrix(ncol = 3))
  names(df) <- c("outbreak", "time", "r_naut")
  total_time = 400
  for (t in 1:total_time) {
    if (t < tau1) {
      beta = beta0
    } else{
      x = -k*(t-tau1)
      beta = (beta0)*(exp(x))
    }
    r_naut = beta/1/7.411374
    df_new <- data.frame(outbreak, t, r_naut)
    names(df_new) <- c("outbreak", "time", "r_naut")
    df <- rbind(df, df_new)
  }
  df <- na.omit(df)
  return(df)
}





