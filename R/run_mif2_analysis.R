############################
## Run MIF2 Althaus model
## Runs MIF2 analysis
############################
source("R/althaus_traj_match.R")



times <- 1:100
rnot <- if_else(times<t_match$params["tau1"], 
                t_match$params["beta0"]/t_match$params["gamma"], 
                exp(-t_match$params["k"]*(times - t_match$params["tau1"]))*t_match$params["beta0"]/t_match$params["gamma"])
plot(times, rnot, type = "l") # Plot matches althaus fig 2

## Do particle filtering of parameters
seir_parm["k"] = t_match$params["k"]
seir_parm["tau1"] = t_match$params["tau1"]
seir_parm["beta0"] = t_match$params["beta0"]

trajectory(pomp(althaus_seir_pomp, params = seir_parm), as.data.frame=TRUE)


pfilter(althaus_seir_pomp, params=seir_parm, Np=1000, filter.mean=TRUE) -> test
logLik(test)


#################################################
## mif2 global search
#################################################
sub_parms <- function(sub_parms=NULL,
                      ref_parms) {
  for(nm in names(sub_parms)) {
    ref_parms[nm] <- sub_parms[nm]
  }
  ref_parms
}
get_parms <- function(ref_parms){
  parm_box <- rbind(
    beta0=c(1e-4,1),
    tau1=c(1e-4,50),
    k=c(1e-4, 1)
  )
  rand_parms <- apply(parm_box, 1, function(x) exp(runif(1, log(x[1]), log(x[2]))))
  sub_parms(rand_parms, ref_parms)
}

date_time <- Sys.time()
dir.create(paste0("data_produced/", date_time))

registerDoMC(cores=2) 

mcopts <- list(set.seed=TRUE)
rw_sd <- 0.002

stew(file=paste0("data_produced/", date_time, "/global_mif_seir.rda"),{
  t_local <- system.time({
    mifs_global <- foreach(i=1:50,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  {
      mif2(
        althaus_seir_pomp,
        start=get_parms(seir_parm),
        Np=100,
        Nmif=10,
        cooling.type="geometric",
        cooling.fraction.50=0.6,
        transform=TRUE,
        rw.sd=rw.sd(
          beta0=rw_sd,
          tau1=rw_sd,
          k=rw_sd)
      )
    }
  })
},seed=808,kind="L'Ecuyer")


stew(file=paste0("data_produced/", date_time, "/global_mif_pf.rda"),{
  t_local_eval <- system.time({
    liks_global <- foreach(i=1:50,.packages='pomp',.combine=rbind) %dopar% {
      evals <- replicate(50, logLik(pfilter(althaus_seir_pomp,params=coef(mifs_global[[i]]),Np=500)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=900242057,kind="L'Ecuyer")

results_global <- data.frame(logLik=liks_global[,1],logLik_se=liks_global[,2],t(sapply(mifs_global,coef)))
summary(results_global$logLik,digits=7)
#pairs(~logLik+beta_I+beta_W, data=results_local)
pairs(~logLik+tau1+k+beta0, data=results_global)




# #################################################
# ## Pfilter tester
# ## Can delete later
# #################################################
# stew(file="data_produced/like_seir.rda",{
#   t_test_eval <- system.time({
#     liks_test <- foreach(i=1:10,.packages='pomp',.combine=rbind) %dopar% {
#       evals <- replicate(10, logLik(pfilter(althaus_seir_pomp,params=seir_parm,Np=1000)))
#       test_mean<-logmeanexp(evals, se=TRUE)
#     }
#   })
# },seed=39048,kind="L'Ecuyer")
# 
# results_test <- data.frame(logLik=liks_test[,1],logLik_se=liks_test[,2])
# summary(results_test$logLik,digits=8)
# 
# #################################################
# ## mif2 Local search
# #################################################
# registerDoMC(cores=2) 
# 
# mcopts <- list(set.seed=TRUE)
# rw_sd <- 0.002
# 
# stew(file="data_produced/local_mif_seir2.rda",{
#   t_local <- system.time({
#     mifs_local <- foreach(i=1:100,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  {
#       mif2(
#         althaus_seir_pomp,
#         start=seir_parm,
#         Np=1000,
#         Nmif=10,
#         cooling.type="geometric",
#         cooling.fraction.50=0.6,
#         transform=TRUE,
#         rw.sd=rw.sd(
#           beta0=rw_sd,
#           tau1=rw_sd,
#           k=rw_sd)
#       )
#     }
#   })
# },seed=808,kind="L'Ecuyer")
# 
# 
# stew(file="likes_mif_loca2l.rda",{
#   t_local_eval <- system.time({
#     liks_local <- foreach(i=1:100,.packages='pomp',.combine=rbind) %dopar% {
#       evals <- replicate(10, logLik(pfilter(althaus_seir_pomp,params=coef(mifs_local[[i]]),Np=1000)))
#       logmeanexp(evals, se=TRUE)
#     }
#   })
# },seed=900242057,kind="L'Ecuyer")
# results_local <- data.frame(logLik=liks_local[,1],logLik_se=liks_local[,2],t(sapply(mifs_local,coef)))
# summary(results_local$logLik,digits=7)
# #pairs(~logLik+beta_I+beta_W, data=results_local)
# pairs(~logLik+tau1+k+beta0, data=results_local)
# 
