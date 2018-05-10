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

sapply(c("R/read_in_drc_data.R","R/ss_pomp_mod.R", "R/helper_functions.R"), source)


outbrk_list <- c("Yambuku", "Kikwit", "Mweka2007","Mweka2008", "Isiro", "Boende")

# Profile Likelihood Bounds
# Parameters are listed in the following order
# beta_lower, beta_upper,p0_lower, p0_upper
bounds <- list(Yambuku=c(1,15,.001,.1),
               Kikwit=c(1,4,.01,.15),
               Mweka2007=c(.5,3.5,.001,.13),
               Mweka2008=c(.1,4,.001,.8),
               Boende=c(.2,4,.001,.6),
               Isiro=c(.2,3,.001,.6))

mod_runner <- function(outbrk_list,dat) {
  
  results_df <- data.frame()
  
# Setup MIF settings ------------------------------------------------------

  for (outbreak in outbrk_list) {
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
                   outbreak = outbreak,
                   model_used = "ss",
                   est_parms = c("beta0", "p0"),
                   parms_sd = rw.sd(beta0=.1, p0=0.02),
                   intensive = TRUE,
                   bounds = bounds[outbreak])
  
  ## First generate the pomp model for the outbreak
    pomp_mod <- generate_pomp_model(outbreak, drc)
  
  ## Now iteratively filter to find MLE
    mif_runs <- mif2_multirun(pomp_obj = pomp_mod, 
                            settings = settings, 
                            refresh = F)

  ## Extract best fit model
    max_mif <- find_max_ll_mif(mif_runs)
    # print(outbreak)
  ## For best fit parameter estimates, calculate the likelihood profile
    prof_lik <- prof_lik_run(mif2_obj = max_mif,
                           settings=settings,
                           refresh = F)

    #plot_prof_lik(prof_lik, settings)
    conf_int <- conf_interval(prof_lik, settings)
    
    #Store results in data frame
    results <- ss_results(outbreak, max_mif, conf_int)
    results_df <- rbind(results_df, results)
  }
  return(results_df)
  }

ss_results <- mod_runner(outbrk_list,drc)

# Write results out to rda
dest <- "data_produced/outbreak_rda/ss_results.rda"
save(ss_results, file = dest)

