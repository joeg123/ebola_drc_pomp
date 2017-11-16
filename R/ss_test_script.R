###################################
## Testing superspreading model
###################################


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

sapply(c("R/read_in_drc_data.R","R/ss_pomp_mod.R"), source)


traj.match(ss_seir_pomp, 
           start=seir_parm,
           method = c("Nelder-Mead"),
           est=c("beta0","p0"),
           transform=TRUE) -> t_match

t_match$params

trajectory(pomp(ss_seir_pomp, params = t_match@params), as.data.frame=T) %>%
  ggplot(aes(time, C)) + geom_line()

logLik(pfilter(pomp(ss_seir_pomp, params = t_match@params), Np=1000))
simulate(pomp(ss_seir_pomp, params = t_match@params), nsim = 100, as.data.frame=T) %>%
  ggplot(aes(time, C, group = sim)) + geom_line()

calc_rnot <- function(fit_parms){
  unname(fit_parms["beta0"] * fit_parms["p0"] / fit_parms["gamma"])
}

calc_cv <- function(fit_parms){
  p <- unname(fit_parms["p0"])
  rnot <- calc_rnot(fit_parms)
  sqrt( (rnot * (2 / p - 1) + 1) / rnot)
}

# MIF2 Analysis -----------------------------------------------------------

sub_parms <- function(sub_parms=NULL,
                      ref_parms) {
  for(nm in names(sub_parms)) {
    ref_parms[nm] <- sub_parms[nm]
  }
  ref_parms
}
get_parms <- function(ref_parms){
  parm_box <- rbind(
    beta0=c(1e-4,10),
    p0=c(1e-8, 1)
  )
  rand_parms <- apply(parm_box, 1, function(x) exp(runif(1, log(x[1]), log(x[2]))))
  sub_parms(rand_parms, ref_parms)
}



get_most_recent_date <- function(dir){
  dirs <- list.dirs(dir, full.names = F)
  dirs[length(dirs)]
}
date_time <- get_most_recent_date("data_produced")


date_time <- Sys.time()
date_time <- paste0("ss-", format(date_time, "%Y-%m-%d"), "-t", format(date_time, "%H%M"))
dir.create(paste0("data_produced/", date_time))




registerDoMC(cores=2) 

mcopts <- list(set.seed=TRUE)
seir_parm <- sub_parms(ref_parms = t_match@params)

stew(file=paste0("data_produced/", date_time, "/global_mif_seir.rda"),{
  t_local <- system.time({
    mifs_global <- foreach(i=1:10,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  {
      mif2(
        ss_seir_pomp,
        start=get_parms(seir_parm),
        Np=10000,
        Nmif=500,
        cooling.type="geometric",
        cooling.fraction.50=0.6,
        transform=TRUE,
        rw.sd=rw.sd(
          beta0=.1,
          p0=0.02)
      )
    }
  })
},seed=80101,kind="L'Ecuyer")


plot(mifs_global)

stew(file=paste0("data_produced/", date_time, "/global_mif_pf.rda"),{
  t_local_eval <- system.time({
    liks_global <- foreach(i=1:10,.packages='pomp',.combine=rbind) %dopar% {
      evals <- replicate(10, logLik(pfilter(ss_seir_pomp,params=coef(mifs_global[[i]]),Np=2000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=900242057,kind="L'Ecuyer")



data.frame(logLik=liks_global[,1], logLik_se=liks_global[,2], t(sapply(mifs_global,coef))) %>%
  filter(logLik > -100) %>%
  gather(key,value, beta0, p0) %>%
  ggplot(aes(value, logLik)) + geom_point() + facet_wrap(~key, scales="free_x")



best_match <- mifs_global[[which.max(map(mifs_global, logLik) %>% flatten_dbl())]]

calc_rnot(best_match$params)
calc_cv(best_match$params)



## Slice design

sliceDesign(
  center=c(beta0 = 1.464, p0=.0642080, sigma = 1/9.312799, 
           gamma = 1/7.411374, ff = 49/69),
  beta0=rep(seq(from=0.5,to=5,length=20),each=10),
  p0=rep(seq(from=.0001,to=.15,length=20),each=10)
) -> p

registerDoParallel()
set.seed(108028909,kind="L'Ecuyer")

foreach (theta=iter(p,"row"),.combine = rbind,
         .inorder = FALSE,
         .options.multicore=list(set.seed=TRUE)
) %dopar% {
  pfilter(ss_seir_pomp,params=unlist(theta),Np=1000) -> pf
  theta$loglik <- logLik(pf)
  theta
} -> p


p %>% group_by(beta0, p0, slice) %>%
  summarize(mll = mean(loglik)) %>%
  gather(parm, value, beta0, p0) %>%
  filter(slice == parm) %>%
  ggplot(aes(x=value,y=mll, color=parm))+
  geom_point()+
  facet_wrap(~parm, scales = "free_x") +
  coord_cartesian(ylim = c(-100, -92))





