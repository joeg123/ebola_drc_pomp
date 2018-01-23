###################################
## Testing superspreading model
###################################

#rm(list=ls())
#sapply(c("R/read_in_drc_data.R","R/ss_pomp_mod.R"), source)

t_match_func <<- function(ss_seir_pomp, graph_traj_sim = FALSE) {
  # browser()
traj.match(ss_seir_pomp, 
           method = c("Nelder-Mead"),
           est=c("beta0","p0"),
           transform=TRUE) ->> t_match



if (graph_traj_sim == TRUE) {
  t_match$params
  logLik(pfilter(pomp(ss_seir_pomp, params = t_match@params), Np=1000))
  trajectory(pomp(ss_seir_pomp, params = t_match@params), as.data.frame=T) %>%
  ggplot(aes(time, C)) + geom_line()
}

if (graph_traj_sim == TRUE) {
simulate(pomp(ss_seir_pomp, params = t_match@params), nsim = 100, as.data.frame=T) %>%
  ggplot(aes(time, C, group = sim)) + geom_line()
}

}

calc_rnot <- function(fit_parms){
  unname(fit_parms["beta0"] * fit_parms["p0"] / fit_parms["gamma"])
}

calc_cv <- function(fit_parms){
  p <- unname(fit_parms["p0"])
  rnot <- calc_rnot(fit_parms)
  sqrt( (rnot * (2 / p - 1) + 1) / rnot)
}

# MIF2 Analysis -----------------------------------------------------------

mif2_run <<- function(ss_seir_pomp, outbreak, graph_mif2 = TRUE) {


ss_seir_pomp <<- ss_seir_pomp
t_match_func(ss_seir_pomp)

sub_parms <<- function(sub_parms=NULL,
                      ref_parms) {
  for(nm in names(sub_parms)) {
    ref_parms[nm] <<- sub_parms[nm]
  }
  ref_parms <<- ref_parms
}

get_parms <<- function(ref_parms){
  parm_box <<- rbind(
    beta0=c(1e-4,10),
    p0=c(1e-8, 1)
  )
  rand_parms <<- apply(parm_box, 1, function(x) exp(runif(1, log(x[1]), log(x[2]))))
  sub_parms(rand_parms, ref_parms)
}



get_most_recent_date <<- function(dir){
  dirs <<- list.dirs(dir, full.names = F)
  dirs[length(dirs)]
}
date_time <- get_most_recent_date("data_produced")


date_time <- Sys.time()
date_time <- paste0("ss-", format(date_time, "%Y-%m-%d"), "-t", format(date_time, "%H%M"))
#dir.create(paste0("data_produced/", date_time))

registerDoMC(cores=2) 

mcopts <<- list(set.seed=TRUE)

seir_parm <<- sub_parms(ref_parms = t_match@params)

outbrk <<- outbreak

dest <- paste0("/", outbrk, "_mif_seir.rda")

stew(file=paste0("data_produced/outbreak_rda", dest),{
  t_local <<- system.time({
    mifs_global <<- foreach(i=1:10,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  {
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

if (graph_mif2 == TRUE) {
  plot(mifs_global)
}

dest <- paste0("/", outbrk, "_mif_pf.rda")

stew(file=paste0("data_produced/outbreak_rda", dest),{
  t_local_eval <- system.time({
    liks_global <- foreach(i=1:10,.packages='pomp',.combine=rbind) %dopar% {
      evals <- replicate(10, logLik(pfilter(ss_seir_pomp,params=coef(mifs_global[[i]]),Np=2000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=900242057,kind="L'Ecuyer")

if (graph_mif2 == TRUE) {
data.frame(logLik=liks_global[,1], logLik_se=liks_global[,2], t(sapply(mifs_global,coef))) %>%
  filter(logLik > -100) %>%
  gather(key,value, beta0, p0) %>%
  ggplot(aes(value, logLik)) + geom_point() + facet_wrap(~key, scales="free_x")
  }

mif2_best_match <- mifs_global[[which.max(map(mifs_global, logLik) %>% flatten_dbl())]]

calc_rnot(mif2_best_match$params)
calc_cv(mif2_best_match$params)

return(mif2_best_match)
}

## Slice design

prof_lik_run <<- function(ss_seir_pomp, beta_naut, p_naut,outbreak) {
outbrk <<- outbreak
beta_naut <<- beta0
p_naut <<- p0

#b0 = 1.464, p0=.0642080 

sliceDesign(
  center=c(beta0 = beta_naut, p0=p_naut, sigma = 1/9.312799, 
           gamma = 1/7.411374, ff = 49/69),
  beta0=rep(seq(from=0.3,to=5,length=20),each=10),
  p0=rep(seq(from=.005,to=.35,length=20),each=10)
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
} -> prof_lik


  prof_lik %>% group_by(beta0, p0, slice) %>%
    summarize(mll = mean(loglik)) %>%
    gather(parm, value, beta0, p0) %>%
    filter(slice == parm) %>%
    ggplot(aes(x=value,y=mll, color=parm))+
    geom_point()+
      facet_wrap(~parm, scales = "free_x") +
      coord_cartesian(ylim = c(-100, -92))
  
  file_name <<- paste0(outbrk,".pdf")
  ggsave(file_name)

}


