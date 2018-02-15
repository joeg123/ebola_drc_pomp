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
              10000, 500,
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

sim_seir_model <- function(times,outbreak="Yambuku") {
  dat <- create_data(times)
  source("R/seir_pomp_mod.R")
  seir_model <- ebola_seir_model(outbreak,dat,TRUE)
}

sim_ss_model <- function(times,outbreak="Yambuku") {
  dat <- create_data(times)
  source("R/ss_pomp_mod.R")
  seir_model <- ebola_seir_model(outbreak,dat,TRUE)
}


sim_data_gen <- function(seir_model,num_sims) {
  sim_keep <- data.frame(times=numeric(0),cases=numeric(0))
  while (TRUE) {
  sims <- simulation_generator(seir_model,10)
  sims %>% mutate(times=time) %>% group_by(sim) %>%
    filter(sum(cases) >= 30) -> sim_temp
    if (length(sim_keep) < 20) {
      sim_keep <- cbind(sim_keep,sim_temp)
    } else {break}
  }
  sim_keep <- as.data.frame(sim_keep)
  View(sim_keep)
}

sim_study <- function(sims, ss_model) {
  for (j in 1:num_sims) {
    out_file <- paste0("sim_", j)
    sims %>% filter(sim==j) -> data
    #par <- mif2_run(ss_mod,out_file,settings)
    #bounds <- prof_lik_run(ss_mod, out_file,par,settings)
  break
}
}


#seir_mod <- sim_seir_model(200)
sim_data <- sim_data_gen(sim_seir_model(200),20)
ss_model <- sim_ss_model(200)
sim_study(sim_data,ss_model)





sim_study_wrong <- function(outbrk_list,dat) {
  num_sims <- 1
  for (outbreak in outbrk_list) {
    source("R/seir_pomp_mod.R")
    seir_model <- ebola_seir_model(outbreak,drc)
    seir_model@params["beta0"] <- unname(1.5 * seir_model@params["gamma"])
    sims <- simulation_generator(seir_model,num_sims)
    sims %>% mutate(times=time) %>% select(times,cases,sim) %>% 
      group_by(cases, times) -> sims
    sims %>% group_by(sim) %>% 
      summarise(outbreak_size = sum(cases)) %>% 
      ggplot(aes(outbreak_size)) + geom_histogram()
    
    for (j in 1:num_sims) {
      source("R/ss_pomp_mod.R")
      out_file <- paste0(outbreak, "_sim_", j)
      sims %>% filter(sim==j) -> data
      data <- as.data.frame(data)
      ss_seir_mod <- ebola_ss_model(outbreak,data)
      #par <- mif2_run(ss_seir_mod,out_file,settings)
      #bounds <- prof_lik_run(ss_seir_mod, out_file,par,settings)
      break
    }
    break
  }
}







