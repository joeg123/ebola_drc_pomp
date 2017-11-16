############################
## Slice Likelihoods
## 
############################



#####
# Log lik profile


source("R/althaus_traj_match.R")

t_match

sliceDesign(
  center=c(beta0 = .824, tau1 = 12.82, k=.133,sigma = 1/9.312799, 
           gamma = 1/7.411374, ff = plogis(49/69), beta1 = 1),
  beta0=rep(seq(from=0.01,to=1.5,length=50),each=10),
  tau1=rep(seq(from=0,to=30,length=50),each=10),
  k=rep(seq(from=.01,to=.30,length=50),each=10)
) -> p

registerDoParallel()
set.seed(108028909,kind="L'Ecuyer")

foreach (theta=iter(p,"row"),.combine = rbind,
         .inorder = FALSE,
         .options.multicore=list(set.seed=TRUE)
) %dopar% {
  pfilter(althaus_seir_pomp,params=unlist(theta),Np=1000) -> pf
  theta$loglik <- logLik(pf)
  theta
} -> p


p %>% 
  gather(parm, value, beta0, tau1, k) %>%
  filter(slice == parm) %>%
  ggplot(aes(x=value,y=loglik, color=parm))+
  geom_point()+
  facet_wrap(~parm, scales = "free_x") +
  geom_smooth()


p %>% group_by(beta0, tau1, k, slice)  %>%
  summarize(loglik = mean(loglik)) %>%
  ungroup() %>%
  mutate(loglik = loglik - max(loglik)) %>%
  gather(key, value, beta0:k) %>%
  filter(slice == key, loglik > -10) %>%
  ggplot(aes(x=value,y=loglik))+
    geom_point()+
    facet_wrap(~key, scales = "free_x") +
    geom_hline(yintercept = -1.92, lty=2)


#########
#End log lik profile
########

########
# Testing pomp mode
########


#Trajectory 
SEIR.traj <- trajectory(seir_pomp, params = c(t_match$params, seir.init.state), as.data.frame=TRUE)

# Create synthetic data

synth_data <- (SEIR.traj[c("time","C")])
#dC <- c(diff(SEIR.traj$C))
dC <- round(SEIR.traj$C)
#dC <- c(0, dC)
synth_data$C <- c(dC)
names(synth_data) <- c("times","cases")


#Deterministic

pomp(data = synth_data,
     times="times",
     t0=0,
     skeleton = vectorfield(seir_skel), 
     params = seir_parm,
     initializer=init,
     statenames= c("S","E","I","R", "D", "C"), 
     paramnames=c("tau1", "beta0","k", "sigma", "gamma", "ff","beta1"),
     zeronames = c("C"),
     fromEstimationScale = untrans,
     toEstimationScale = trans,
     dmeasure = seir_dmeasure,
     rmeasure = seir_rmeasure,
     rprocess = euler.sim(step.fun = seir_rprocess, delta.t = 1)) -> seir_pomp_test

traj.match(seir_pomp_test, 
           start=seir_parm,
           method = c("Nelder-Mead"),
           est=c("tau1", "beta0","k"),
           transform=TRUE) -> t_match_test

t_match_test$params
trajectory(seir_pomp_test, params=t_match_test$params, times=0:120, as.data.frame=TRUE) %>% 
  ggplot(aes(x=time, y=C)) + geom_line() +
  geom_point(data=synth_data, aes(x=times, y = cases), color = "red")

# Simulate using new rmeasure

sim <- simulate(seir_pomp, params = c(t_match_test$params, seir.init.state), nsim = 10, as.data.frame=TRUE)

sim %>% 
  ggplot(aes(x=time,y=C, group=sim))+
  geom_line()


# no clue how to graph this



#       sigma        gamma           ff         tau1        beta0        beta1            k 
#1.073791e-01 1.349277e-01 6.704332e-01 3.535753e+01 3.131401e-01 1.000000e+00 3.212365e-09 
#

#log lik of simulated data  

sliceDesign(
  center=c(beta0 = .7454, tau1 = 12.8366, k=.1249,sigma = 1/9.312799, 
           gamma = 1/7.411374, ff = plogis(49/69), beta1 = 1),
  beta0=rep(seq(from=0,to=2.0,length=20),each=10),
  tau1=rep(seq(from=0,to=50,length=20),each=10),
  k=rep(seq(from=.01,to=.30,length=20),each=10)
) -> p

registerDoParallel()
set.seed(108028909,kind="L'Ecuyer")

foreach (theta=iter(p,"row"),.combine = rbind,
         .inorder = FALSE,
         .options.multicore=list(set.seed=TRUE)
) %dopar% {
  library (pomp)
  pfilter(seir_pomp_test,params=unlist(theta),Np=500) -> pf_test
  theta$loglik <- logLik(pf_test)
  theta
} -> p_t


p_t %>% 
  gather(parm, value, beta0, tau1, k) %>%
  filter(slice == parm) %>%
  ggplot(aes(x=value,y=loglik, color=parm))+
  geom_point()+
  facet_wrap(~parm, scales = "free_x") +
  geom_smooth()


## Maximizing likelihood using pfilter
# pfilter give stochastic est of liklihood
# must use derivative free algorithm to 
# optimizing over unbounded parms (must transform)

# fixing a reference point in parameter space
neg.ll <- function (par, est) {
  allpars <- coef(flu,transform=TRUE)
  allpars[est] <- par
  try(
    freeze(
      pfilter(flu,params=partrans(flu,allpars,dir="fromEst"),
              Np=2000),
      seed=915909831
    )
  ) -> pf
  if (inherits(pf,"try-error")) 1e10 else -logLik(pf)
}

# Nelder mead optimizes function (neg.ll)
fit <- optim(
  par=c(),
  est=c(),
  fn=negloglik,
  method="Nelder-Mead",
  control=list(maxit=400,trace=0)
)

mle <- flu
coef(mle,c("Beta","mu_I","rho"),transform=TRUE) <- fit$par
coef(mle)

