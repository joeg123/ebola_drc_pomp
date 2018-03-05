####
### Functions that can be used by any pomp models
####

pull_outbreak_data <- function(outbreak=c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende"),
                               data = NULL){
  ## Sets up the data based on the given outbreak name
  outbrk <- match.arg(outbreak)
  data %>%
    ungroup() %>%
    filter(outbreak==outbrk) %>%
    select(times,cases) -> data
  
  data <- as.data.frame(data)  
}

randomize_parms <- function(pomp_obj, model_used){
  ## Randomizes parameters for mif2 starting iteration
  if(model_used == "combo"){
    parm_box <- rbind(
      beta0 = c(1e-4, 10),
      p0 = c(1e-8, 1),
      k = c(1e-5, 1),
      tau1 = c(1, max(pomp_mod@times))
    )
  } else if (model_used == "ss"){
    parm_box <- rbind(
      beta0 = c(1e-4, 10),
      p0 = c(1e-8, 1)
    )
  } else if (model_used == "int"){
    parm_box <- rbind(
      beta0 = c(1e-4, 10),
      k = c(1e-5, 1),
      tau1 = c(1, max(pomp_mod@times))
    )
  } else{
    stop("Do not recognize the model_used parameter. Either a typo, or need to add new parameter randomization for new model")
  }
  rand_parms <- apply(parm_box, 1, function(x) exp(runif(1, log(x[1]), log(x[2]))))
  sub_parms(rand_parms, coef(pomp_obj))
}

mif2_single_run <- function(parms,
                            pomp_obj, 
                            settings) {
  ## Function to run single mif
  ## With the specified parameters
  
  num_particles <- settings$mif_nparticles
  num_mif_iterations <- settings$mif_niter
  
  mif2(pomp_obj,
       start=parms,
       Np=num_particles,
       Nmif=num_mif_iterations,
       cooling.type="geometric",
       cooling.fraction.50=0.6,
       transform=TRUE,
       rw.sd= settings$parms_sd)
}

mif2_multirun <- function(pomp_obj, 
                          settings,
                          refresh = FALSE,
                          n = 10){
  ## Function to run multiple mifs
  
  ## First setup the location for the .rda file
  dest <- paste0("data_produced/outbreak_rda/", settings$model_used,"_", settings$outbreak, "_mif.rda")
  
  ## Generate list of parms for running models
  parm_list <- n %>% rerun(randomize_parms(pomp_obj, settings$model_used))
  # parm_list <- n %>% rerun(coef(pomp_obj))
  if(file.exists(dest) & refresh == FALSE){
    load(dest)
  } else{
    mifs_global <- parm_list %>% map(mif2_single_run, 
                                     pomp_obj=pomp_obj, 
                                     settings=settings)
    mifs_global <- as(mifs_global, "mif2List")
    save(mifs_global, file = dest)
  }
  
  return(mifs_global)
}
  
find_max_ll_mif <- function(mif2_list){
  ## Returns the MLE mif run
  mif2_list[[which.max(map(mif2_list, logLik) %>% flatten_dbl())]]
}


get_parm_bounds <- function(parm, mif2_obj){
  lower <- mif2_obj@params[parm] / 5
  upper <- mif2_obj@params[parm] * 5
  # if(parm == "beta0"){
  #   c(0.1, 15)
  # } else if(parm == "k"){
  #   c(0, 5)
  # }else if(parm == "p0"){
  #   c(1e-4, 0.99)
  # } else if(parm == "tau1"){
  #   c(0, 200)
  # } else{
  #   stop("Parameter not recognized, please check est_parms for accuracy")
  # }
  c(lower, upper)
}

get_list_parm_bounds <- function(parms, mif2_obj){
  parms %>% map(get_parm_bounds, mif2_obj=mif2_obj) %>% setNames(parms)
}

get_parm_slice_input <- function(parm, parm_bounds, settings){
  rep(seq(from=parm_bounds[[parm]][1], to=parm_bounds[[parm]][[2]], length=settings$slice_length), each=settings$slice_reps)
}

get_lik_slice <- function(mif2_obj, settings){
  ## Gives the slice design used for the likelihood profile
  est_parms <- settings$est_parms
  model_used <-settings$model_used
  
  parm_bounds <- get_list_parm_bounds(est_parms, mif2_obj)

  if(model_used == "combo"){
    sliceDesign(
      center=c(coef(mif2_obj)["beta0"],coef(mif2_obj)["p0"], 
               coef(mif2_obj)["tau1"], coef(mif2_obj)["k"],
               sigma = 1/9.312799, gamma = 1/7.411374, ff = 49/69),
      beta0 = get_parm_slice_input("beta0", parm_bounds, settings),
      p0 = get_parm_slice_input("p0", parm_bounds, settings),
      k = get_parm_slice_input("k", parm_bounds, settings),
      tau1 = get_parm_slice_input("tau1", parm_bounds, settings)
    )
  } else if (model_used == "ss"){
    sliceDesign(
      center=c(coef(mif2_obj)["beta0"],coef(mif2_obj)["p0"], 
               sigma = 1/9.312799, gamma = 1/7.411374, ff = 49/69),
      beta0 = get_parm_slice_input("beta0", parm_bounds, settings),
      p0 = get_parm_slice_input("p0", parm_bounds, settings)
    )
  } else if (model_used == "int"){
    sliceDesign(
      center=c(coef(mif2_obj)["beta0"], 
               coef(mif2_obj)["tau1"], coef(mif2_obj)["k"],
               sigma = 1/9.312799, gamma = 1/7.411374, ff = 49/69),
      beta0 = get_parm_slice_input("beta0", parm_bounds, settings),
      k = get_parm_slice_input("k", parm_bounds, settings),
      tau1 = get_parm_slice_input("tau1", parm_bounds, settings)
    )
  } else{
    stop("Do not recognize the model_used parameter. Either a typo, or need to add new parameter randomization for new model")
  }
}

get_single_lik <- function(pomp_obj, num_particles, ...){
  ## Completes the particle filtering for a certain set of parameters
  ## The parameters are specified by the ... and called through the prof_like_run function
  
  pfilter(pomp_obj, params=unlist(list(...)), Np = num_particles)
}

prof_lik_run <- function(mif2_obj, est_parms, settings, outbreak, model_used, refresh=FALSE){
  ## This function calculates the likelihood profile for the MLE parameters, and
  ## Returns a data_frame with the results (and saves it for later)
  
  dest <- paste0("data_produced/outbreak_rda/", settings$model_used, "_", settings$outbreak, "_prof.rda")
  if(file.exists(dest) & refresh == FALSE){
    load(dest)
  } else{
    parms <- get_lik_slice(mif2_obj, settings) %>% as_data_frame
    slice_runs <- parms %>% 
      pmap(.f = get_single_lik, mif2_obj, settings$prof_lik_nparticles) 
    
    prof_lik <- parms %>% mutate(ll = map(slice_runs, logLik) %>% unlist)
    save(prof_lik, file = dest)
  }
  prof_lik
}

plot_prof_lik <- function(df, est_parms){
  df %>% gather(key, value, est_parms) %>%
    filter(key==slice) %>%
    ggplot(aes(x=value, y=ll))+
    geom_point()+
    facet_wrap(~key, scales = "free_x") +
    geom_smooth()
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
  ref_parms 
}



calc_rnot <- function(fit_parms){
  unname(fit_parms["beta0"] * fit_parms["p0"] / fit_parms["gamma"])
}

calc_cv <- function(fit_parms){
  p <- unname(fit_parms["p0"])
  rnot <- calc_rnot(fit_parms)
  sqrt( (rnot * (2 / p - 1) + 1) / rnot)
}









# Old code ----------------------------------------------------------------
# prof_lik_run <- function(pomp_obj, outbreak, pars, settings) {
#   
#   prof_file <- paste0(outbreak, "_prof_lik.rda")
#   
#   parallel_vars <- new.env()
#   assign("pomp_obj", pomp_obj, envir= parallel_vars)
#   assign("num_particle1", settings[5], envir= parallel_vars)
#   
#   p <- pars[1]
#   b <- pars[2]
#   
#   p_lims <- calc_interval(p)
#   b_lims <- calc_interval(b)
#   
#   if (file.exists(prof_file)) {
#     load(file = prof_file)
#   } else {
#     sliceDesign(
#       center=c(beta0 = b , p0=p, sigma = 1/9.312799, 
#                gamma = 1/7.411374, ff = 49/69),
#       beta0=rep(seq(from=b_lims[1],to=b_lims[2],length=settings[6]),each=settings[7]),
#       p0=rep(seq(from=p_lims[1],to=p_lims[2],length=settings[6]),each=settings[7])
#     ) -> p
#     
#     #Length Used for current SS and SEIR prof lik 
#     # length = 250
#     # each = 100
#     # issue: prof lik is very jagged for some outbreaks
#     
#     cl <- makeCluster(settings[1])
#     clusterExport(cl,c("pomp_obj", "num_particle1"), envir = parallel_vars)
#     registerDoParallel(cl)
#     
#     set.seed(108028909,kind="L'Ecuyer")
#     
#     foreach (theta=iter(p,"row"),.packages = 'pomp',.combine = rbind,
#              .inorder = FALSE,
#              .options.multicore=list(set.seed=TRUE)
#     ) %dopar% {
#       pfilter(pomp_obj,params=unlist(theta),Np=num_particle1) -> pf
#       theta$loglik <- logLik(pf)
#       theta
#     } -> prof_lik
#     prof_lik
#     save(prof_lik, file = prof_file)
#     stopCluster(cl)
#     stopImplicitCluster()
#     closeAllConnections()
#   }
#   
#   
#   prof_lik %>% group_by(beta0, p0, slice)  %>%
#     summarize(loglik = mean(loglik)) %>%
#     ungroup() %>%
#     mutate(loglik = loglik - max(loglik)) %>%
#     gather(key, value, beta0:p0) %>%
#     filter(slice == key, loglik > -10) %>%
#     filter(slice == "beta0") -> spline_dat_beta0
#   
#   beta_results <- conf_int_find(spline_dat_beta0$value, spline_dat_beta0$loglik,-1.92)
#   
#   prof_lik %>% group_by(beta0, p0, slice)  %>%
#     summarize(loglik = mean(loglik)) %>%
#     ungroup() %>%
#     mutate(loglik = loglik - max(loglik)) %>%
#     gather(key, value, beta0:p0) %>%
#     filter(slice == key, loglik > -10) %>%
#     filter(slice == "p0") -> spline_dat_p0
#   
#   p_results <- conf_int_find(spline_dat_p0$value,spline_dat_p0$loglik,-1.92)
#   
#   prof_out <- c(round(beta_results[1],digits = 2),round(beta_results[2], digits = 2),
#                 round(p_results[1], digits = 2),round(p_results[2], digits = 2))
#   
#   rm(prof_file)
#   return(prof_out)
#   
# }
# mif2_run_old <- function(pomp_obj, outbreak, settings) {
#   ## Old version kept for archival purposes
#   parallel_vars <- new.env()
#   assign("pomp_obj", pomp_obj, envir= parallel_vars)
#   assign("num_particle1", settings[2], envir= parallel_vars)
#   assign("num_mif", settings[3], envir= parallel_vars)
# 
# 
#   traj.match(pomp_obj,
#              method = c("Nelder-Mead"),
#              est=c("beta0","p0"),
#              transform=TRUE) -> t_match
#   assign("t_match", t_match, envir= parallel_vars)
# 
#   seir_parm <- sub_parms(ref_parms = t_match@params)
#   seir_parm <- get_parms(seir_parm)
#   assign("seir_parm", seir_parm, envir= parallel_vars)
# 
#   dest <- paste0("/", outbreak, "_mif_seir.rda")
# 
# 
#   # Np and Nmif Suggestion
#   # Np=10000
#   # Nmif=500
# 
#   cl <- makeCluster(settings[1])
#   clusterExport(cl,c("pomp_obj","t_match","seir_parm","num_particle1","num_mif"), envir = parallel_vars)
#   registerDoParallel(cl)
# 
#   stew(file=paste0("data_produced/outbreak_rda", dest),{
#     t_local <- system.time({
#       mifs_global <- foreach(i=1:10,.packages='pomp', .combine=c, .options.multicore=list(set.seed=TRUE)) %dopar% {
#         mif2(
#           pomp_obj,
#           start=seir_parm,
#           Np=num_particle1,
#           Nmif=num_mif,
#           cooling.type="geometric",
#           cooling.fraction.50=0.6,
#           transform=TRUE,
#           rw.sd=rw.sd(
#             beta0=.1,
#             p0=0.02)
#         )
#       }
#     })
#   },seed=80101,kind="L'Ecuyer")
#   browser()
#   stopCluster(cl)
#   stopImplicitCluster()
# 
#   dest <- paste0("/", outbreak, "_mif_pf.rda")
# 
#   parallel_vars <- new.env()
#   assign("mifs_global", mifs_global, envir = parallel_vars)
#   assign("pomp_obj", pomp_obj, envir= parallel_vars)
#   assign("num_particle2", settings[4], envir= parallel_vars)
# 
#   cl <- makeCluster(settings[1])
#   clusterExport(cl,c("pomp_obj","mifs_global", "num_particle2"), envir = parallel_vars)
#   registerDoParallel(cl)
# 
#   # Np Suggestion
#   # Np = 2000
# 
#   stew(file=paste0("data_produced/outbreak_rda", dest),{
#     t_local_eval <- system.time({
#       liks_global <- foreach(i=1:10,.packages='pomp',.combine=rbind) %dopar% {
#         evals <- replicate(10, logLik(pfilter(pomp_obj,params=coef(mifs_global[[i]]),Np=num_particle2)))
#         logmeanexp(evals, se=TRUE)
#       }
#     })
#   },seed=900242057,kind="L'Ecuyer")
# 
#   stopCluster(cl)
#   stopImplicitCluster()
#   closeAllConnections()
# 
# 
#   #plot(mifs_global)
# 
#   mif2_best_match <- mifs_global[[which.max(map(mifs_global, logLik) %>% flatten_dbl())]]
# 
#   LL <- mif2_best_match$loglik
# 
#   cv <- calc_cv(mif2_best_match$params)
#   r_0 <- calc_rnot(mif2_best_match$params)
#   k <- calc_k(mif2_best_match$params)
#   out_par <- ((mif2_best_match$params))
#   out_par <- c(unname(out_par['p0']),unname(out_par['beta0']), r_0, cv, k, LL)
#   return(out_par)
# 
# }
# 
# get_parms <- function(ref_parms){
#   parm_box <- rbind(
#     beta0=c(1e-4,10),
#     p0=c(1e-8, 1)
#   )
#   rand_parms <- apply(parm_box, 1, function(x) exp(runif(1, log(x[1]), log(x[2]))))
#   sub_parms(rand_parms, ref_parms)
# }
# 
# calc_interval <- function(val) {
#   lower_lim <- val / 5
#   upper_lim <- val * 5
#   lims <- c(lower_lim, upper_lim)
#   return(lims)
# }
# 
# # None of this works well... I think maybe we should do it manually -------
# 
# 
# 
# ## Plot the likelihood profile
# plot_prof_lik(prof_lik, est_parms)
# 
# 
# find_scam_root <- function(df, lower = T){
#   solver_fxn <- function(x, mod, mod_val_at_max){
#     ## Add in the 1.96, because likelihoods are negative, and we need to find -1.96
#     unname(predict(mod, data_frame(var = x))) - mod_val_at_max + 1.96
#   }
#   print(lower)
#   if(lower){
#     df <- df %>% filter(var <= df$var[which.max(df$rel_ll)] & rel_ll >= -100)  
#     browser()
#     mod <- scam(rel_ll ~ s(var, k = nrow(df), bs="mpd"), data = df)
#     root <- try(uniroot(f = solver_fxn, 
#                         interval = c(min(df$var), max(df$var)), 
#                         mod = mod, 
#                         mod_val_at_max = unname(predict(mod, data_frame(var = min(df$var))))))
#     
#   } else{
#     df <- df %>% filter(var >= df$var[which.max(df$rel_ll)] & rel_ll >= -100)  
#     mod <- scam(rel_ll ~ s(var, k = nrow(df), bs="mpd"), data = df)
#     root <- try(uniroot(f = solver_fxn, 
#                         interval = c(min(df$var), max(df$var)), 
#                         mod = mod, 
#                         mod_val_at_max = unname(predict(mod, data_frame(var = min(df$var))))))
#   }
#   
#   if(class(root) == "try-error"){
#     warning("May need to increase the range of parameters investigated, check plots")
#     if_else(lower, 0, Inf)
#   } else{
#     root$root  
#   }
# }
# 
# calc_single_parm_ci <- function(parm, prof_lik){
#   library(scam)
#   ## First average across pfilter likelihoods, and find the relative ll
#   df <- prof_lik %>% filter(slice==parm) %>% 
#     group_by_at(vars(-slice, -ll)) %>% 
#     summarize(avg_ll = mean(ll)) %>% 
#     ungroup() %>% 
#     mutate(rel_ll = avg_ll - max(avg_ll, na.rm=T)) %>% 
#     select(var = parm, rel_ll)
#   
#   ## Now find the lowerbound
#   lower <- find_scam_root(df, lower=TRUE)
#   
#   ## Now find upperbound
#   upper <- find_scam_root(df, lower=FALSE)
#   
#   return(data_frame(parm=parm, lower = lower, upper = upper))
# }
# 
# calc_all_parm_ci <- function(prof_lik, est_parms, outbreak, model_used){
#   ## Calculates each parameters confidence interval from the prof_lik
#   est_parms %>% map(calc_single_parm_ci, prof_lik) %>% bind_rows()
# }
