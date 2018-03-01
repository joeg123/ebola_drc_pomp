############################
## Trajectory Matching Althaus Model
## Uses trajectory matching to fit Althaus pomp to data
############################
rm(list=ls())

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

sapply(c("R/read_in_drc_data.R","R/althaus_pomp_mod_code.R"), source)


traj.match(althaus_seir_pomp, 
           start=seir_parm,
           method = c("Nelder-Mead"),
           est=c("tau1", "beta0","k"),
           transform=TRUE) -> t_match

t_match$params

