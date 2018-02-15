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
library(dplyr)
library(magrittr)
library(pryr)


sapply(c("R/read_in_drc_data.R", "R/main_functions.R"), source)

outbrk_list <- c("Yambuku","Mweka2007","Mweka2008","Isiro","Boende")

########
# Settings
# 1: Number of cores
# 2: MIF2 Np
# 3: MIF2 Nmif
# 4: Log Lik Np
# 5: Profile Likelihood Np
# 6: Slice Length
# 7: Slice Each
########

settings <- c(2,
              100, 50,
              100, 100,
              150, 100)


simulation_generator <- function(model,num_sims) {
  model <- model
  sims <- simulate(model,
                   nsim=num_sims,as.data.frame=TRUE,include.data=FALSE)
  return(sims)
}

create_data <- function(size) {
  times <- seq(1,size,1)
  cases <- integer(size)
  data <- data.frame(times=times,cases=cases)
}

sim_seir_model <- function(times) {
  dat <- create_data(times)
  source("R/seir_pomp_mod.R")
  seir_model <- ebola_seir_model("Yambuku",dat,TRUE)
}

sim_ss_model <- function(dat) {
  source("R/ss_pomp_mod.R")
  ss_model <- ebola_ss_model("Yambuku",dat)
}

threshold_sim_gen <- function(pomp_mod, case_threshold){
  while (TRUE) {
    sims <- simulation_generator(pomp_mod,1)
    if (sum(sims$cases) >= case_threshold) {
      sims <- sims %>% select(times=time, cases)
      return(sims)
    }
  }
}

sim_data_gen <- function(pomp_mod,num_sims, case_threshold) {
  sim_keep <- vector("list", num_sims)
  for(i in 1:num_sims){
    sim_keep[[i]] <- threshold_sim_gen(pomp_mod, case_threshold)
  }
  return(sim_keep)
}

sim_study <- function(num_sims) {
  sim_data <- sim_data_gen(sim_seir_model(200),num_sims,30)
  for (j in 1:num_sims) {
    out_file <- paste0("sim_", j)
    pomp_mod <- sim_ss_model(sim_data[[j]])
    par <- mif2_run(pomp_mod,out_file,settings)
    bounds <- prof_lik_run(pomp_mod, out_file,par,settings)
}
}

sim_study(20)







