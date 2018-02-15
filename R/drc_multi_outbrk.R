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




sapply(c("R/read_in_drc_data.R","R/ss_pomp_mod.R", "R/main_functions.R"), source)

outbrk_list <- c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende")
multi_drc <- function(outbrk_list,dat) {

mif2_results <- matrix(0,6,6)
i <- 1

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
              2000, 200,
              1000, 1000,
              250, 150)

for (outbreak in outbrk_list) {
    i
    mif2_results
    print(outbreak)
    ss_seir_pomp <- ebola_ss_model(outbreak,dat)
    par <- mif2_run(ss_seir_pomp,outbreak,settings)
    bounds <- prof_lik_run(ss_seir_pomp, outbreak,par,settings)
    p_out <- paste0(round(par[1], digits = 2)," (",bounds[3],"-",bounds[4], ")")
    beta_out <- paste0(round(par[2], digits = 2)," (",bounds[1],"-",bounds[2], ")")
    mif2_results[i,1] <- p_out
    mif2_results[i,2] <- beta_out
    mif2_results[i,3] <- round(par[3], digits = 2)
    mif2_results[i,4] <- round(par[4], digits = 2)
    mif2_results[i,5] <- round(par[5], digits = 2)
    mif2_results[i,6] <- round(par[6], digits = 2)
    i <- i + 1


}
return(mif2_results)

}

max_mif2 <- multi_drc(outbrk_list,drc)

dimnames(max_mif2) <- list(outbrk_list,c('p','B_0','R_0', 'CV', 'k', 'LL'))

as.data.frame(max_mif2) %>% xtable(display = c("s","s","s", "fg", "fg", "fg", "fg"))



