#######
## Intervention Model
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




sapply(c("R/read_in_drc_data.R","R/int_pomp_mod.R", "R/helper_functions.R"), source)

outbrk_list <- c("Yambuku", "Kikwit","Mweka2007","Mweka2008","Isiro","Boende")


mod_runner <- function(outbrk_list,dat) {
  
  i <- 1
  
  for (outbreak in outbrk_list) {
    # Settings
    # 1: Number of cores
    # 2: MIF2 Np
    # 3: MIF2 Nmif
    # 4: Log Lik Np
    # 5: Profile Likelihood Np
    # 6: Slice Length
    # 7: Slice Each
    
    settings <- list(num_cores = 1,
                     mif_nparticles = 2000, 
                     mif_niter = 2000,
                     prof_lik_nparticles = 1000,
                     slice_length = 100, 
                     slice_reps = 100,
                     outbreak = outbreak,
                     model_used = "int",
                     est_parms = c("beta0", "k", "tau1"),
                     parms_sd = rw.sd(beta0=.1, k=0.02, tau1=.1))
    
    
    ## First generate the pomp model for the outbreak
    pomp_mod <- generate_pomp_model(outbreak, drc)
    
    ## Now iteratively filter to find MLE
    mif_runs <- mif2_multirun(pomp_obj = pomp_mod, 
                              settings = settings, 
                              refresh = F)

    ## Extract best fit model
    max_mif <- find_max_ll_mif(mif_runs)
    
    print(outbreak)
    print(max_mif@params)
    print(max_mif@loglik)
    ## For best fit parameter estimates, calculate the likelihood profile
    prof_lik <- prof_lik_run(mif2_obj = max_mif,
                             settings=settings,
                             refresh = F)
    
    prof_lik %>% gather(key,val, settings$est_parms) %>% 
      filter(key == slice) %>% 
      group_by(key, val) %>% 
      summarize(avg_ll = mean(ll)) %>% 
      mutate(avg_ll = avg_ll - max(avg_ll)) %>% 
      ggplot(aes(val, avg_ll)) + facet_wrap(~key, scales="free_x") + 
      geom_line() +
      coord_cartesian(ylim = c(-100,0)) +
      stat_smooth()

  }
  return(mif_runs)
  }

max_mif2_alt <- mod_runner(outbrk_list,drc)



#dimnames(max_mif2_alt) <- list(outbrk_list,c('B_0','R_0', 'CV', 'log'))

#as.data.frame(max_mif2_alt) %>% xtable(display = c("fg","fg", "fg", "fg", "fg"))
