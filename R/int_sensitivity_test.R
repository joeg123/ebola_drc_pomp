#####
## Sensitivity Testing 
#####

rm(list = ls())

library(pomp)
library(chron)
library(lubridate)
library(tidyverse)
library(magrittr)
library(cowplot)
library(padr)
library(reshape2)
library(dplyr)
library(knitr)
library(DT)
library(xtable)
library(scam)

sapply(c("R/read_in_drc_data.R","R/ss_pomp_mod.R", "R/helper_functions.R"), source)

outbrk_list <- c("Yambuku", "Kikwit", "Mweka2007", "Isiro", "Boende", "Mweka2008")
est_parms <- c("beta0", "k", "tau1")

mod_runner <- function(outbreak,dat) {
  bounds <- list(Yambuku=c(1,15,.001,.1),
                 Kikwit=c(1,4,.01,.15),
                 Mweka2007=c(.5,3.5,.001,.13),
                 Mweka2008=c(.1,4,.001,.8),
                 Boende=c(.2,4,.001,.6),
                 Isiro=c(.2,3,.001,.6))
  results_df <- data.frame()
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

  return(max_mif)
}

#### Boende Sensitivity Test

outbrk <- "Boende"
outbreak <- "Boende"

#10 Additional cases + 25 days
shift <- 25

drc %>% filter(outbreak==outbrk) %>%
  mutate(times=times+shift) %>%
  ungroup() %>% 
  add_row(outbreak=outbrk,
          date_infection = min(drc$date_infection[drc$outbreak==outbrk])-shift,
          cases = 5, times = 1, .before=1)-> df

#### Mweka2008 Sensitivity Test

outbrk2 <- "Mweka2008"
outbreak2 <- "Mweka2008"

#10 Additional cases + 25 days
shift <- 25

drc %>% filter(outbreak==outbrk2) %>%
  mutate(times=times+shift) %>%
  ungroup() %>% 
  add_row(outbreak=outbrk2,
          date_infection = min(drc$date_infection[drc$outbreak==outbrk2])-shift,
          cases = 5, times = 1, .before=1)-> df2

#### Yambuku Sensitivity Test

outbrk3 <- "Yambuku"
outbreak3 <- "Yambuku"

#10 Additional cases + 25 days
shift <- 25

drc %>% filter(outbreak==outbrk3) %>%
  mutate(times=times+shift) %>%
  ungroup() %>% 
  add_row(outbreak=outbrk3,
          date_infection = min(drc$date_infection[drc$outbreak==outbrk3])-shift,
          cases = 5, times = 1, .before=1)-> df3

boende_result <- mod_runner(outbrk, df)
Mweka_result <- mod_runner(outbrk2, df2)
Yambuku_result <- mod_runner(outbrk3, df3)

b_ll <- boende_result@loglik
b_aic <- 6 - 2*b_ll
m_ll <- Mweka_result@loglik
m_aic <- 6 - 2*m_ll
y_ll <- Yambuku_result@loglik
y_aic <- 6 - 2*y_ll

b_ll
b_aic
m_ll
m_aic 
y_ll
y_aic

