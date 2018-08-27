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

outbrk_list <- c("Yambuku", "Kikwit", "Mweka2007", "Isiro", "Boende", "Mweka2008", "Equator")

est_parms <- c("beta0", "k", "tau1")

# Profile Likelihood Bounds
# Parameters are listed in the following order
# beta_lower, beta_upper,k_lower, k_upper, tau_lower, tau_upper

bounds <- list(Yambuku=c(.2,1.8,.08,.2,8,20),
               Kikwit=c(.25,.75,.06,.09,100,130),
               Mweka2007=c(.1,.3, 0, 1,125,175),
               Mweka2008=c(.1,3,.01,.15,0,20),
               Boende=c(.1,1.8,.03,.15,.01,18),
               Isiro=c(.05,.45,.001,.02,.001,50),
               Equator=c(.5,5,.001, 2.5, .001, 10))


mod_runner <- function(outbrk_list,dat) {
  
  results_df <- data.frame()
  
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
    print(outbreak)
    
    settings <- list(num_cores = 1,
                     mif_nparticles = 2000, 
                     mif_niter = 2000,
                     prof_lik_nparticles = 1000,
                     slice_length = 100, 
                     slice_reps = 50,
                     outbreak = outbreak,
                     model_used = "int",
                     est_parms = est_parms,
                     parms_sd = rw.sd(beta0=.1, k=0.02, tau1=.1),
                     intensive=TRUE,
                     bounds=bounds[outbreak])
    
    print("Generating the pomp model...")
    ## First generate the pomp model for the outbreak
    pomp_mod <- generate_pomp_model(outbreak, drc)
    
    print("Fitting the model to the outbreak data...")
    ## Now iteratively filter to find MLE
    mif_runs <- mif2_multirun(pomp_obj = pomp_mod, 
                              settings = settings, 
                              refresh = F)
    
    print("Starting the profile likelihood...")
    ## Extract best fit model
    max_mif <- find_max_ll_mif(mif_runs)
    # print(outbreak)
    # browser()
    ## For best fit parameter estimates, calculate the likelihood profile
    prof_lik <- prof_lik_run(mif2_obj = max_mif,
                             settings=settings,
                             refresh = F)
    
    print("Calculating the parameter confidence intervals...")
    plot_prof_lik(prof_lik, max_mif, settings)
    conf_int <- conf_interval(prof_lik, max_mif, settings)

    #Store results in data frame
    results <- int_results(outbreak, max_mif, conf_int)
    results_df <- rbind(results_df, results)
  }
  return(results_df)
}

int_output <- mod_runner(outbrk_list,drc)

# Write results out to rda
save(int_output, file = "data_produced/outbreak_rda/parm_est_int.rda")


##extra code
# testeroni <- mod_runner(c("Mweka2007", "Kikwit", "Mweka2008"), drc)
# 
# testeroni <- mod_runner(c("Mweka2007", "Isiro"), drc)
# 
# 
# load("data_produced/outbreak_rda/int_Mweka2007_prof.rda")
# prof_lik %>% gather(key,val, beta0, tau1, k) %>% 
#   filter(key == slice) %>% 
#   group_by(key, val) %>% 
#   summarize(avg_ll = mean(ll)) %>%
#   mutate(avg_ll = avg_ll - max(avg_ll)) %>%
#   # mutate(avg_ll = ll-max(ll)) %>%
#   ggplot(aes(val, avg_ll)) + 
#   facet_wrap(~key, scales="free_x") + 
#   geom_point(alpha=1,size=.5) +
#   coord_cartesian(ylim = c(-30,0)) +
#   geom_hline(yintercept=-1.96)
