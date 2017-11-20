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

sapply(c("R/read_in_drc_data.R","R/ss_pomp_mod.R", "R/ss_test_script.R"), source)

outbrk_list <- c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende")

multi_drc <- function(outbrk_list,dat,mif2=FALSE,log_prof=FALSE) {

for (outbreak in outbrk_list) {
  if (mif2 == TRUE) {
    ss_seir_pomp <- ebola_ss_model(outbreak,dat)
    t_match_func(ss_seir_pomp)
    #mif2_run(ss_seir_pomp)
  }
  if (log_prof == TRUE) {
    prof_lik_run(ss_seir_pomp)
  }
  return(t_match)
}
}



multi_drc(outbrk_list,drc,mif2 = TRUE)

