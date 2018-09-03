#######
## Calls all outbreaks and applies pomp, tmatch, mif2, and log profile functions
######
rm(list = ls())

library(deSolve)
library(mvtnorm)
library(pomp)
library(chron)
library(lubridate)
library(tidyverse)
library(cowplot)
library(foreach)
library(doMC)
library(doParallel)
library(padr)
library(reshape2)
library(dplyr)
library(knitr)
library(DT)
library(xtable)

sapply(c("R/read_in_drc_data.R","R/ss_combo_pomp_mod.R", "R/helper_functions.R"), source)


# # outbrk_list <- c("Yambuku", "Kikwit", "Mweka2007","Mweka2008", "Isiro", "Boende")
# 
# # Profile Likelihood Bounds
# # Parameters are listed in the following order
# # beta_lower, beta_upper,p0_lower, p0_upper
# bounds <- list(Yambuku=c(1,15,.001,.1),
#                Kikwit=c(1,4,.01,.15),
#                Mweka2007=c(.5,3.5,.001,.13),
#                Mweka2008=c(.1,4,.001,.8),
#                Boende=c(.2,4,.001,.6),
#                Isiro=c(.2,3,.001,.6))

# mod_runner <- function(outbrk_list,dat) {
#   
#   results_df <- data.frame()
#   
#   for (outbreak in outbrk_list) {
    # Settings
    # 1: Number of cores
    # 2: MIF2 Np
    # 3: MIF2 Nmif
    # 4: Profile Likelihood Np
    # 5: Slice Length
    # 6: Slice Each
    # 7: Outbreak
    # 8: Model Used
    # 9: Model Parameters
    # 10: Parameter Standard Deviation
    # 11: Intensive Profile Likelihood?
    # 12: Bounds for the the Profile Likelihood when intensive
settings <- list(num_cores = 1,
                 mif_nparticles = 2000, 
                 mif_niter = 2000,
                 prof_lik_nparticles = 1000,
                 slice_length = 100, 
                 slice_reps = 50,
                 outbreak = "combo",
                 model_used = "ss",
                 est_parms = c("beta0", "p0"),
                 parms_sd = rw.sd(beta0=.1, p0=0.02),
                 intensive = TRUE,
                 bounds = c(0.5,5, .01,.2))
  
print("Generating the pomp model...")
## First generate the pomp model for the outbreak
pomp_mod <- generate_pomp_model(drc)
  
print("Fitting the model to the outbreak data...")
## Now iteratively filter to find MLE
mif_runs <- mif2_multirun(pomp_obj = pomp_mod, 
                          settings = settings, 
                          refresh = T)

print("Starting the profile likelihood...")
## Extract best fit model
max_mif <- find_max_ll_mif(mif_runs)
# print(outbreak)
## For best fit parameter estimates, calculate the likelihood profile
prof_lik <- prof_lik_run(mif2_obj = max_mif,
                         settings=settings,
                         refresh = T)

print("Calculating the parameter confidence intervals...")
plot_prof_lik(prof_lik, max_mif, settings)
conf_int <- conf_interval(prof_lik, max_mif, settings)

#Store results in data frame
results <- ss_results(settings$outbreak, max_mif, conf_int)


# Write results out to rda
save(results, file = "data_produced/outbreak_rda/parm_est_ss_combo.rda")



# Investigate fit with the data -------------------------------------------

sim_epidemics <- simulate(max_mif, nsim = 1000, as.data.frame=T) 

outbrk_data <- drc %>% padr::pad() %>% 
  ungroup() %>% 
  mutate(day = seq_along(times)) %>% 
  filter(!is.na(times)) %>% 
  #group_by(outbreak) %>% 
  mutate(cum_cases = cumsum(cases)) %>% 
  ungroup() %>% 
  select(outbreak, day, cum_cases, times)

## Looks reasonable
sim_epidemics %>% 
  as_tibble() %>% 
  left_join(outbrk_data, by = c("time" = "day")) %>% 
  group_by(sim) %>% 
  mutate(cum_sim_cases = cumsum(cases)) %>% 
  group_by(time) %>% 
  summarize(lwr = quantile(cum_sim_cases, 0.025),
            upr = quantile(cum_sim_cases, 0.975),
            times = first(times),
            cum_cases = first(cum_cases)) %>% 
  ggplot(aes(x = time, cum_cases)) + 
  geom_ribbon(aes(ymin= lwr, ymax = upr), alpha = .2) +  
  geom_line() #+
  #facet_wrap(~outbreak, ncol = 1, scales = "free_y")



outbrk_data <- drc %>% padr::pad() %>% 
  ungroup() %>% 
  mutate(day = seq_along(times)) %>% 
  filter(!is.na(times)) %>% 
  group_by(outbreak) %>% 
  mutate(cum_cases = cumsum(cases)) %>% 
  ungroup() %>% 
  select(outbreak, day, cum_cases, times)
  
sim_epidemics %>% 
  as_tibble() %>% 
  left_join(outbrk_data, by = c("time" = "day")) %>% 
  group_by(sim,outbreak) %>% 
  mutate(cum_sim_cases = cumsum(cases)) %>% 
  group_by(time, outbreak) %>% 
  summarize(lwr = quantile(cum_sim_cases, 0.025),
            upr = quantile(cum_sim_cases, 0.975),
            times = first(times),
            cum_cases = first(cum_cases)) %>% 
  ggplot(aes(x = times, cum_cases)) + 
    geom_ribbon(aes(ymin= lwr, ymax = upr), alpha = .2) +  
    geom_line() +
    facet_wrap(~outbreak, ncol = 1, scales = "free_y")


