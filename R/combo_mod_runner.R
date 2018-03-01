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
              2000, 
              1000,
              1000, 
              1000,
              250, 
              150)


# Running mif on all analysis ---------------------------------------------

outbrk_list <- c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende")


mif2_run <- function(pomp_obj, 
                     outbreak, 
                     settings,
                     est_parms, 
                     parms_sd,
                     model_used,
                     refresh=FALSE) {
  
  
  
  parallel_vars <- new.env()
  assign("pomp_obj", pomp_obj, envir= parallel_vars)
  assign("num_particle1", settings[2], envir= parallel_vars)
  assign("num_mif", settings[3], envir= parallel_vars)
  
  traj.match(pomp_obj, 
             method = c("Nelder-Mead"),
             est= est_parms,
             transform=TRUE) -> t_match
  assign("t_match", t_match, envir= parallel_vars)
  
  seir_parm <- sub_parms(ref_parms = t_match@params)
  seir_parm <- get_parms(seir_parm)
  assign("seir_parm", seir_parm, envir= parallel_vars)
  
  dest <- paste0("/", model_used,"_", outbreak, "_mif.rda")
  
  cl <- makeCluster(settings[1])
  clusterExport(cl,c("pomp_obj","t_match","seir_parm","num_particle1","num_mif", "parms_sd"), envir = parallel_vars)
  registerDoParallel(cl)
  

  # Run the global mif step to find the best fit parameters -----------------
  if(refresh){
    file.remove(paste0("data_produced/outbreak_rda", dest))
  }
  stew(file=paste0("data_produced/outbreak_rda", dest),{
    t_local <- system.time({
      mifs_global <- foreach(i=1:10,.packages='pomp', .combine=c, .options.multicore=list(set.seed=TRUE)) %dopar% {
        mif2(
          pomp_obj,
          start=seir_parm,
          Np=num_particle1,
          Nmif=num_mif,
          cooling.type="geometric",
          cooling.fraction.50=0.6,
          transform=TRUE,
          rw.sd= parms_sd
        )
      }
    })
  },seed=80101,kind="L'Ecuyer")
  
  stopCluster(cl)
  stopImplicitCluster()
  
  dest <- paste0("/", model_used,"_", outbreak, "_pf.rda")
  
  parallel_vars <- new.env()
  assign("mifs_global", mifs_global, envir = parallel_vars)
  assign("pomp_obj", pomp_obj, envir= parallel_vars)
  assign("num_particle2", settings[4], envir= parallel_vars)
  
  cl <- makeCluster(settings[1])
  clusterExport(cl,c("pomp_obj","mifs_global", "num_particle2"), envir = parallel_vars)
  registerDoParallel(cl)
  

  # What is point of this section? ------------------------------------------

  if(refresh){
    file.remove(paste0("data_produced/outbreak_rda", dest))
  }
  stew(file=paste0("data_produced/outbreak_rda", dest),{
    t_local_eval <- system.time({
      liks_global <- foreach(i=1:10,.packages='pomp',.combine=rbind) %dopar% {
        evals <- replicate(10, logLik(pfilter(pomp_obj, params=coef(mifs_global[[i]]),Np=num_particle2)))
        logmeanexp(evals, se=TRUE)
      }
    })
  },seed=900242057,kind="L'Ecuyer")
  
  stopCluster(cl)
  stopImplicitCluster()
  closeAllConnections()
  
  return(mifs_global)
  
  
  # mif2_best_match <- mifs_global[[which.max(map(mifs_global, logLik) %>% flatten_dbl())]]
  # LL <- mif2_best_match$loglik
  # cv <- calc_cv(mif2_best_match$params)
  # r_0 <- calc_rnot(mif2_best_match$params)
  # k <- calc_k(mif2_best_match$params)
  # out_par <- ((mif2_best_match$params))
  # out_par <- c(unname(out_par['p0']),unname(out_par['beta0']), r_0, cv, k, LL)
  # return(out_par)
}


pomp_mod <- generate_pomp_model("Boende", drc)

est_parms <- c("beta0", "p0", "k", "tau1")
parms_sd <- rw.sd(beta0=.1, p0=0.02, k=0.02, tau1=.1)
model_used <- "combo"
par <- mif2_run(pomp_obj = pomp_mod, 
                outbreak = "Mweka2007", 
                settings = settings, 
                est_parms = est_parms, 
                parms_sd = parms_sd, 
                model_used = model_used,
                refresh = T)

  

## Not looked at yet


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


