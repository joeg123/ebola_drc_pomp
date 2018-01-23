####
### Creates mif2 and pf rda files
####

calc_interval <- function(val) {
  lower_lim <- 0
  upper_lim <- val*2.5
  lims <- c(lower_lim, upper_lim)
  return(lims)
}


calc_k <- function(parms) {
  p <- unname(parms["p0"])
  R <- calc_rnot(parms)
  p <- 1-p
  f <- function(x)  ((1+(R/x))^(-x)-p)
  uniroot(f, lower=0.00000001, upper= 2)$root
}


conf_int_find <- function(x, y, log_val) {
  # beta0/p0 = x
  # loglik = y
  lower_bound <- 0
  upper_bound <- 0
  index <- length(y)
  max_log <- which.max(y)
  for (i in 2:max_log) {
    y1 <- y[i-1]
    y2 <- y[i]
    if ((y2 > log_val) && (y1 < log_val)){
      x1 <- x[i-1]
      x2 <- x[i]
      lower_bound <- x1 + abs((log_val - y1)*((x2-x1)/(y2-y1)))
    }
  }

  for (i in max_log+1:index) {
    y1 <- y[i-1]
    y2 <- y[i]
    if (is.na(y1) | is.na(y2)) {break}
    if ((y1 > log_val) && (y2 < log_val)){
      x1 <- x[i-1]
      x2 <- x[i]
      upper_bound <- x1 + abs((log_val - y1)*((x2-x1)/(y2-y1)))
    }
  }
  bounds <- c(lower_bound,upper_bound)
  return(bounds)
}

sub_parms <- function(sub_parms=NULL,
                      ref_parms) {
  for(nm in names(sub_parms)) {
    ref_parms[nm] <- sub_parms[nm]
  }
  ref_parms <- ref_parms
}

get_parms <- function(ref_parms){
  parm_box <- rbind(
    beta0=c(1e-4,10),
    p0=c(1e-8, 1)
  )
  rand_parms <- apply(parm_box, 1, function(x) exp(runif(1, log(x[1]), log(x[2]))))
  sub_parms(rand_parms, ref_parms)
}

calc_rnot <- function(fit_parms){
  unname(fit_parms["beta0"] * fit_parms["p0"] / fit_parms["gamma"])
}

calc_cv <- function(fit_parms){
  p <- unname(fit_parms["p0"])
  rnot <- calc_rnot(fit_parms)
  sqrt( (rnot * (2 / p - 1) + 1) / rnot)
}


mif2_run <- function(ss_seir_pomp, outbreak, graph_mif2 = TRUE) {
  ss_seir_pomp <<- ss_seir_pomp
  
  calc_cv
  calc_rnot
  traj.match(ss_seir_pomp, 
             method = c("Nelder-Mead"),
             est=c("beta0","p0"),
             transform=TRUE) -> t_match

  
  mcopts <- list(set.seed=TRUE)
  
  outbrk <- outbreak
  
  dest <- paste0("/", outbrk, "_mif_seir.rda")
  
  registerDoMC(cores=2)

  stew(file=paste0("data_produced/outbreak_rda", dest),{
    t_local <- system.time({
      mifs_global <- foreach(i=1:10,.packages='pomp', .combine=c, .options.multicore='mcopts') %dopar% {
        seir_parm <- sub_parms(ref_parms = t_match@params)
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
  
  mifs_global <<- mifs_global
  
  #plot(mifs_global)

  dest <- paste0("/", outbrk, "_mif_pf.rda")
  
  stew(file=paste0("data_produced/outbreak_rda", dest),{
    t_local_eval <- system.time({
      liks_global <- foreach(i=1:10,.packages='pomp',.combine=rbind) %dopar% {
        evals <- replicate(10, logLik(pfilter(ss_seir_pomp,params=coef(mifs_global[[i]]),Np=2000)))
        logmeanexp(evals, se=TRUE)
      }
    })
  },seed=900242057,kind="L'Ecuyer")
  
  liks_global
  
  mif2_best_match <- mifs_global[[which.max(map(mifs_global, logLik) %>% flatten_dbl())]]
  
  LL <- mif2_best_match$loglik
  
  cv <- calc_cv(mif2_best_match$params)
  r_0 <- calc_rnot(mif2_best_match$params)
  k <- calc_k(mif2_best_match$params)
  #print(k)
  out_par <- ((mif2_best_match$params))
  out_par <- c(unname(out_par['p0']),unname(out_par['beta0']), r_0, cv, k, LL)
  
  return(out_par)

}

sim_study <- function(ss_seir_pomp, outbrk,beta0,p0) {
  ss_seir_pomp <- ss_seir_pomp
  b <- beta0
  p <- p0
  sims <- simulate(ss_seir_pomp,params=c(beta0 = b , p0=p, sigma = 1/9.312799, 
                   gamma = 1/7.411374, ff = 49/69),
                   nsim=40,as.data.frame=TRUE,include.data=TRUE)
  print(sims)
  
}


prof_lik_run <- function(ss_seir_pomp, outbrk, beta0, p0) {
  outbrk
  ss_seir_pomp
  conf_int_find
  b <- beta0
  p <- p0
  prof_file <- paste0(outbrk, "_prof_lik.rda")
  #b0 = 1.464, p0=.0642080 
  b_lims <- calc_interval(b)
  p_lims <- calc_interval(p)
  
  if (file.exists(prof_file)) {
    load(file = prof_file)
  } else {
  sliceDesign(
    center=c(beta0 = b , p0=p, sigma = 1/9.312799, 
             gamma = 1/7.411374, ff = 49/69),
    beta0=rep(seq(from=b_lims[1],to=b_lims[2],length=250),each=100),
    p0=rep(seq(from=p_lims[1],to=p_lims[2],length=250),each=100)
  ) -> p
  
  # length = 80
  # each = 150
  
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
  prof_lik
  save(prof_lik, file = prof_file)
  }
  
  prof_lik %>% group_by(beta0, p0, slice)  %>%
    summarize(loglik = mean(loglik)) %>%
    ungroup() %>%
    mutate(loglik = loglik - max(loglik)) %>%
    gather(key, value, beta0:p0) %>%
    filter(slice == key, loglik > -10) %>%
    filter(slice == "beta0") -> spline_dat_beta0
  
  beta_results <- conf_int_find(spline_dat_beta0$value,spline_dat_beta0$loglik,-1.92)
  
  prof_lik %>% group_by(beta0, p0, slice)  %>%
    summarize(loglik = mean(loglik)) %>%
    ungroup() %>%
    mutate(loglik = loglik - max(loglik)) %>%
    gather(key, value, beta0:p0) %>%
    filter(slice == key, loglik > -10) %>%
    filter(slice == "p0") -> spline_dat_p0
  
  p_results <- conf_int_find(spline_dat_p0$value,spline_dat_p0$loglik,-1.92)
  
  prof_out <- c(round(beta_results[1],digits = 2),round(beta_results[2], digits = 2),
                round(p_results[1], digits = 2),round(p_results[2], digits = 2))
  
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
  
  file_name <- paste0(outbrk,".pdf")

  ggsave(file_name)
  rm(prof_file)
  return(prof_out)
  
}









