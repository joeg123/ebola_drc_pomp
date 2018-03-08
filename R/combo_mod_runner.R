#######
## Calls all outbreaks and applies pomp, tmatch, mif2, and log profile functions
######


# First setup the environment ---------------------------------------------
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

sapply(c("R/read_in_drc_data.R", "R/combo_pomp_mod.R", "R/helper_functions.R"), source)

# Setup MIF settings ------------------------------------------------------

# Settings
# 1: Number of cores
# 2: MIF2 Np
# 3: MIF2 Nmif
# 4: Log Lik Np
# 5: Profile Likelihood Np
# 6: Slice Length
# 7: Slice Each

settings <- list(num_cores = 1,
              mif_nparticles = 10, 
              mif_niter = 10,
              prof_lik_nparticles = 20,
              slice_length = 50, 
              slice_reps = 10,
              outbreak = "Boende",
              model_used = "combo",
              est_parms = c("beta0", "p0", "k", "tau1"),
              parms_sd = rw.sd(beta0=.1, p0=0.02, k=0.02, tau1=.1))


# Running mif on all analysis (not yet) ---------------------------------------------

# outbrk_list <- c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende")

## First generate the pomp model for the outbreak
pomp_mod <- generate_pomp_model("Boende", drc)

## Now iteratively filter to find MLE
mif_runs <- mif2_multirun(pomp_obj = pomp_mod, 
                          settings = settings, 
                          refresh = T)

## Extract best fit model
max_mif <- find_max_ll_mif(mif_runs)

## For best fit parameter estimates, calculate the likelihood profile
prof_lik <- prof_lik_run(mif2_obj = max_mif,
                         settings=settings,
                         refresh = T)

prof_lik %>% gather(key,val, est_parms) %>% 
  filter(key == slice) %>% 
  group_by(key, val) %>% 
  summarize(avg_ll = mean(ll)) %>% 
  mutate(avg_ll = avg_ll - max(avg_ll)) %>% 
  ggplot(aes(val, avg_ll)) + facet_wrap(~key, scales="free_x") + 
  geom_line() +
  coord_cartesian(ylim = c(-100,0)) +
  stat_smooth()


## This doesn't work and shouldn't be run yet.
calc_all_parm_ci(prof_lik, est_parms, outbreak, model_used)



# Testing trajectory match ------------------------------------------------

combined <- generate_pomp_model("Mweka2007", drc)
traj.match(combined, 
           method = c("Nelder-Mead"),
           est=c("tau1", "beta0","k", "p0"),
           transform=TRUE) -> t_match

logLik(t_match)

t_match@params
# # Not looked at yet -------------------------------------------------------
# run_fitting <- function(outbreak = c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende"), 
#                         dat) {
#     pomp_mod <- generate_pomp_model("Mweka2007", drc)
#     par <- mif2_run(pomp_mod, outbreak, settings)
#     # bounds <- prof_lik_run(pomp_mod, outbreak, par, settings)
#     browser()
#     
#   return(list(par, bounds))
# }
# 
# 
# # run mif on single outbreak ----------------------------------------------
# outbreak <- "Boende"
# 
# test <- run_fitting(outbreak, drc)
# 
# max_mif2 <- multi_drc(outbrk_list, drc)
# 
# dimnames(max_mif2) <- list(outbrk_list,c('p','B_0','R_0', 'CV', 'k', 'LL'))
# 
# as.data.frame(max_mif2) %>% xtable(display = c("s","s","s", "fg", "fg", "fg", "fg"))
# 
# 
