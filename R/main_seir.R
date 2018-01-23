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




sapply(c("R/read_in_drc_data.R","R/seir_pomp_mod.R", "R/main_func_seir.R"), source)

#

outbrk_list_full <- c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende")

outbrk_list <- c("Yambuku","Mweka2007","Mweka2008","Isiro","Boende")
multi_drc_seir <- function(outbrk_list,dat) {
  
  mif2_results <- matrix(0,5,4)
  i <- 1
  
  for (outbreak in outbrk_list) {
    i
    mif2_results
    print(outbreak)
    seir_pomp <- ebola_seir_model(outbreak,dat)
    seir_pomp
    par <- mif2_run_seir(seir_pomp,outbreak)
    bounds <- prof_lik_run(seir_pomp, outbreak, par[1])
    print(bounds[1])
    beta_out <- paste0(round(par[1], digits = 2)," (",bounds[1],"-",bounds[2], ")")
    mif2_results[i,1] <- beta_out
    mif2_results[i,2] <- round(par[2], digits = 2)
    mif2_results[i,3] <- round(par[3], digits = 2)
    mif2_results[i,4] <- round(par[4], digits = 2)
    i <- i + 1
    
    
    
  }
  return(mif2_results)
  
}

max_mif2_seir <- multi_drc_seir(outbrk_list,drc)


dimnames(max_mif2_seir) <- list(outbrk_list,c('B_0','R_0', 'CV', 'LL'))

as.data.frame(max_mif2_seir) %>% xtable(display = c("fg","fg", "fg", "fg", "fg"))



