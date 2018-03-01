##############
# Purpose: This script reads in mif2 and profile likelihood rda files and plots them
# Author: Joseph Garcia
# Date: 03/1/18
# Notes:
##############
rm(list = ls())

library(pomp)
library(lubridate)
library(tidyverse)
library(cowplot)
library(padr)
library(reshape2)
library(dplyr)
library(knitr)
library(xtable)
library(ggplot2)

setwd("/Users/Joseph/Desktop/Epi Lab/ebola_drc_pomp/")

load(file = "mif2_results.rda")

load_obj <- function(f)
{
  setwd("/Users/Joseph/Desktop/Epi Lab/ebola_drc_pomp/")
  env <- new.env()
  rda_name <- load(f, env)
  browser()
  params=coef(rda_name)
  rm(params)
}

outbrk_list <- c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende")

mif_rda_data <- list()
mif_pf_rda_data <- list()
prof_rda_data <- list()

i <- 1

for (outbreak in outbrk_list) {

mif_data_name <- paste0(outbreak, "_mif.rda")
mif_rda_data[[i]] <- load_obj(mif_data_name)

mif_pf_data <- paste0(outbreak, "_mif_pf.rda")
mif_pf_rda_data[[i]] <- load_obj(mif_pf_data)

prof_data <- paste0(outbreak, "_prof_lik.rda")
prof_rda_data[[i]] <- load_obj(prof_data)

break
}


# Table: Maximum MIF2 Results

dimnames(max_mif2) <- list(outbrk_list,c('p','B_0','R_0', 'CV', 'k', 'LL'))

as.data.frame(max_mif2) %>% xtable(display = c("s","s","s", "fg", "fg", "fg", "fg"))

# Plot: Plots MIF2 Global
plot(mifs_global)

# Plot: Profile likelihood of Beta and p

prof_lik %>% group_by(beta0, p0, slice)  %>%
  summarize(loglik = mean(loglik)) %>%
  ungroup() %>%
  mutate(loglik = loglik - max(loglik)) %>%
  gather(key, value, beta0:p0) %>%
  filter(slice == key, loglik > -10) %>%    
  ggplot(aes(x=value,y=loglik))+
  geom_point()+ geom_line() +
  facet_wrap(~key, scales = "free_x") +
  geom_hline(yintercept = -1.92, lty=2)

file_name <- paste0(outbreak,".pdf")

#ggsave(file_name)

