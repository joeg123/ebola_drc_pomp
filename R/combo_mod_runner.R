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



# Testing trajectory match ------------------------------------------------

combined <- generate_pomp_model("Mweka2007", drc)
traj.match(combined, 
           method = c("Nelder-Mead"),
           est=c("tau1", "beta0","k", "p0"),
           transform=TRUE) -> t_match

logLik(t_match)

t_match@params

# Setup MIF settings ------------------------------------------------------

# Settings
# 1: Number of cores
# 2: MIF2 Np
# 3: MIF2 Nmif
# 4: Log Lik Np
# 5: Profile Likelihood Np
# 6: Slice Length
# 7: Slice Each

settings <- c(2,
              10, 10,
              20, 20,
              50, 50)


# Running mif on all analysis ---------------------------------------------

outbrk_list <- c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende")



outbreak <- "Boende"
model_used <- "combo"
est_parms <- c("beta0", "p0", "k", "tau1")
parms_sd <- rw.sd(beta0=.1, p0=0.02, k=0.02, tau1=.1)

## First generate the pomp model for the outbreak
pomp_mod <- generate_pomp_model("Boende", drc)

## Now iteratively filter to find MLE
mif_runs <- mif2_multirun(pomp_obj = pomp_mod, 
                outbreak = outbreak, 
                settings = settings, 
                est_parms = est_parms, 
                parms_sd = parms_sd, 
                model_used = model_used,
                refresh = T)

## Extract best fit model
max_mif <- find_max_ll_mif(mif_runs)



# None of this works well... I think maybe we should do it manually -------


## For best fit parameter estimates, calculate the likelihood profile
prof_lik <- prof_lik_run(mif2_obj = max_mif, 
                         est_parms = est_parms, 
                         outbreak = outbreak, 
                         settings=settings, 
                         model_used = model_used,
                         refresh = T)

## Plot the likelihood profile
plot_prof_lik(prof_lik, est_parms)


find_scam_root <- function(df, lower = T){
  solver_fxn <- function(x, mod, mod_val_at_max){
    ## Add in the 1.96, because likelihoods are negative, and we need to find -1.96
    unname(predict(mod, data_frame(var = x))) - mod_val_at_max + 1.96
  }
  print(lower)
  if(lower){
    df <- df %>% filter(var <= df$var[which.max(df$rel_ll)] & rel_ll >= -100)  
    browser()
    mod <- scam(rel_ll ~ s(var, k = nrow(df), bs="mpd"), data = df)
    root <- try(uniroot(f = solver_fxn, 
                        interval = c(min(df$var), max(df$var)), 
                        mod = mod, 
                        mod_val_at_max = unname(predict(mod, data_frame(var = min(df$var))))))
    
  } else{
    df <- df %>% filter(var >= df$var[which.max(df$rel_ll)] & rel_ll >= -100)  
    mod <- scam(rel_ll ~ s(var, k = nrow(df), bs="mpd"), data = df)
    root <- try(uniroot(f = solver_fxn, 
                        interval = c(min(df$var), max(df$var)), 
                        mod = mod, 
                        mod_val_at_max = unname(predict(mod, data_frame(var = min(df$var))))))
  }
  
  if(class(root) == "try-error"){
    warning("May need to increase the range of parameters investigated, check plots")
    if_else(lower, 0, Inf)
  } else{
    root$root  
  }
}

calc_single_parm_ci <- function(parm, prof_lik){
  library(scam)
  ## First average across pfilter likelihoods, and find the relative ll
  df <- prof_lik %>% filter(slice==parm) %>% 
    group_by_at(vars(-slice, -ll)) %>% 
    summarize(avg_ll = mean(ll)) %>% 
    ungroup() %>% 
    mutate(rel_ll = avg_ll - max(avg_ll, na.rm=T)) %>% 
    select(var = parm, rel_ll)
  
  ## Now find the lowerbound
  lower <- find_scam_root(df, lower=TRUE)
  
  ## Now find upperbound
  upper <- find_scam_root(df, lower=FALSE)
  
  return(data_frame(parm=parm, lower = lower, upper = upper))
}

calc_all_parm_ci <- function(prof_lik, est_parms, outbreak, model_used){
  ## Calculates each parameters confidence interval from the prof_lik
  est_parms %>% map(calc_single_parm_ci, prof_lik) %>% bind_rows()
}

calc_all_parm_ci(prof_lik, est_parms, outbreak, model_used)


prof_lik %>% gather(key,val, est_parms) %>% 
  filter(key == slice) %>% 
  group_by(key, val) %>% 
  summarize(avg_ll = mean(ll)) %>% 
  mutate(avg_ll = avg_ll - max(avg_ll)) %>% 
  ggplot(aes(val, avg_ll)) + facet_wrap(~key, scales="free_x") + 
  geom_line() +
  #  coord_cartesian(ylim = c(-100,0)) +
  stat_smooth()



# Not looked at yet -------------------------------------------------------
run_fitting <- function(outbreak = c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende"), 
                        dat) {
    pomp_mod <- generate_pomp_model("Mweka2007", drc)
    par <- mif2_run(pomp_mod, outbreak, settings)
    # bounds <- prof_lik_run(pomp_mod, outbreak, par, settings)
    browser()
    
  return(list(par, bounds))
}


# run mif on single outbreak ----------------------------------------------
outbreak <- "Boende"

test <- run_fitting(outbreak, drc)

max_mif2 <- multi_drc(outbrk_list, drc)

dimnames(max_mif2) <- list(outbrk_list,c('p','B_0','R_0', 'CV', 'k', 'LL'))

as.data.frame(max_mif2) %>% xtable(display = c("s","s","s", "fg", "fg", "fg", "fg"))


