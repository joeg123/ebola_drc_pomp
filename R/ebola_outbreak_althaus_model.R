########################
##### POMP SKEL
########################

# init
rm(list = ls()) 

# Necessary libraries
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

# Initialize random number generator
set.seed(529342)

# Read the data
ebola <- read.csv("data/Ebola_outbreak_DRC_data.csv")
ebola$Date <- chron(as.character(ebola$Date), format=c(dates = "day mon year"))     
ebola$Cases[1] <- 0
data <- na.omit(ebola[c("Date","Cases")])
names(data) <- c("times","cases")
begin <- chron("26 Jul 2014", format=c(dates = "day mon year")) 
data$times <- data$times - data$times[1]
data$times <- as.numeric(data$times)
#POMP Skel

seir_skel <- Csnippet('
                      double beta, x, N;
                      if (t < tau1) {
                      beta = beta0;
                      } else{
                      x = -k*(t-tau1);
                      beta = (beta0)*(exp(x));
                      }
                      N = S+E+I+R+D;
                      // Balance the equations
                      DS = -beta/N*S*I;
                      DE = beta/N*S*I - sigma*E;
                      DI = sigma*E - gamma*I;
                      DR = (1-ff)*gamma*I;
                      DD = ff*gamma*I;
                      DC = sigma*E;
                      ')

# dt = 1 -> daily data

seir_rprocess <- Csnippet('double beta, x, N;
                          if (t < tau1) {
                            beta = beta0;
                          } else{
                            x = -k*(t-tau1);
                            beta = (beta0)*(exp(x));
                          }
                          N = S+E+I+R+D;
                          
                          double new_exposed = rbinom(S, 1 - exp(-beta/N*I*dt));
                          double new_infected = rbinom(E, 1 - exp(-sigma*dt));
                          double new_not_infectious = rbinom(I, 1 - exp(-gamma*dt));       
                          double new_recovered = rbinom(new_not_infectious, (1 - ff) );
                          double new_dead = new_not_infectious - new_recovered;

                          S -= new_exposed;
                          E += new_exposed - new_infected;
                          I += new_infected - new_recovered - new_dead;
                          R += new_recovered;
                          D += new_dead;
                          C += new_infected;
                          ')

seir_dmeasure <- Csnippet('double f;
                          f = dpois(nearbyint(cases), C,1);
                          lik = (give_log) ? f : exp(f);
                          ')

seir_rmeasure <- Csnippet('
                          cases = rpois(C);
                          ')

trans <- Csnippet('
                    Tbeta0 = log(beta0);
                    Tbeta1 = log(beta1);
                    Tk = log(k);
                    Tff = plogis(ff,0,1,1,0);
                    //Ttau1 = log(tau1);
                  ')

untrans <- Csnippet('
                    Tbeta0 = exp(beta0);
                    Tbeta1 = exp(beta1);
                    Tk = exp(k);
                    Tff = qlogis(ff,0,1,1,0);
                    //Ttau1 = exp(tau1);
                    ')



pop <- 1e6      
#init <- c(S = N - 1, E = 0, I = 1, R = 0, D = 0, C = 0)

init <- Csnippet("
          S= 999999; 
          E = 0.0;
          I = 1.0;
          R = 0.0;
          D = 0.0;
          C = 0.0;
                 ")

seir.init.state <- c(S.0= 999999, 
                     E.0 = 0.0,
                     I.0 = 1.0,
                     R.0 = 0.0,
                     D.0 = 0.0,
                     C.0 = 0.0)

seir_parm <- c(sigma = 1/9.312799, 
               gamma = 1/7.411374, 
               ff = plogis(49/69),
               tau1 = 14.0,
               beta0 = 0.4,
               beta1 = 1,
               k = 1e-4)


pomp(data = data,
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
  rprocess = euler.sim(step.fun = seir_rprocess, delta.t = 1)) -> seir_pomp

#trajectory matching

traj.match(seir_pomp, 
           start=seir_parm,
           method = c("Nelder-Mead"),
           est=c("tau1", "beta0","k"),
           transform=TRUE) -> t_match

t_match$params


# Should match Althaus results
#!beta0!      beta1          !k!          f       tau0       !tau1!      sigma      gamma 
#-0.3634729       -Inf -2.0883390  0.8960880       -Inf 14.2890158  0.1073791  0.1349277

#####
# Log lik profile
#####

sliceDesign(
  center=c(beta0 = .7454, tau1 = 12.8366, k=.1249,sigma = 1/9.312799, 
           gamma = 1/7.411374, ff = plogis(49/69), beta1 = 1),
  beta0=rep(seq(from=0.01,to=2.0,length=20),each=10),
  tau1=rep(seq(from=0,to=30,length=20),each=10),
  k=rep(seq(from=.01,to=.30,length=20),each=10)
) -> p

registerDoParallel()
set.seed(108028909,kind="L'Ecuyer")

foreach (theta=iter(p,"row"),.combine = rbind,
         .inorder = FALSE,
         .options.multicore=list(set.seed=TRUE)
) %dopar% {
           library (pomp)
  pfilter(seir_pomp,params=unlist(theta),Np=500) -> pf
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



###############
# End of testing
################







































times <- 1:100
rnot <- if_else(times<t_match$params["tau1"], 
                t_match$params["beta0"]/t_match$params["gamma"], 
                exp(-t_match$params["k"]*(times - t_match$params["tau1"]))*t_match$params["beta0"]/t_match$params["gamma"])
plot(times, rnot, type = "l") # Plot matches althaus fig 2

## Do particle filtering of parameters
seir_parm["k"] = t_match$params["k"]
seir_parm["tau1"] = t_match$params["tau1"]
seir_parm["beta0"] = t_match$params["beta0"]

trajectory(pomp(seir_pomp, params = seir_parm), as.data.frame=TRUE)


pfilter(seir_pomp, params=seir_parm, Np=1000, filter.mean=TRUE) -> test
logLik(test)

#################################################
## Pfilter tester
## Can delete later
#################################################
stew(file="data_produced/like_seir.rda",{
  t_test_eval <- system.time({
    liks_test <- foreach(i=1:10,.packages='pomp',.combine=rbind) %dopar% {
      evals <- replicate(10, logLik(pfilter(seir_pomp,params=seir_parm,Np=1000)))
      test_mean<-logmeanexp(evals, se=TRUE)
    }
  })
},seed=39048,kind="L'Ecuyer")

results_test <- data.frame(logLik=liks_test[,1],logLik_se=liks_test[,2])
summary(results_test$logLik,digits=8)

#################################################
## mif2 Local search
#################################################
registerDoMC(cores=2) 

mcopts <- list(set.seed=TRUE)
rw_sd <- 0.002

stew(file="data_produced/local_mif_seir2.rda",{
  t_local <- system.time({
    mifs_local <- foreach(i=1:100,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  {
      mif2(
        seir_pomp,
        start=seir_parm,
        Np=1000,
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


stew(file="likes_mif_loca2l.rda",{
  t_local_eval <- system.time({
    liks_local <- foreach(i=1:100,.packages='pomp',.combine=rbind) %dopar% {
      evals <- replicate(10, logLik(pfilter(seir_pomp,params=coef(mifs_local[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=900242057,kind="L'Ecuyer")
results_local <- data.frame(logLik=liks_local[,1],logLik_se=liks_local[,2],t(sapply(mifs_local,coef)))
summary(results_local$logLik,digits=7)
#pairs(~logLik+beta_I+beta_W, data=results_local)
pairs(~logLik+tau1+k+beta0, data=results_local)


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

stew(file="data_produced/global_mif_seir3.rda",{
  t_local <- system.time({
    mifs_local <- foreach(i=1:1000,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  {
      mif2(
        seir_pomp,
        start=get_parms(seir_parm),
        Np=1000,
        Nmif=100,
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


stew(file="likes_mif_global3.rda",{
  t_local_eval <- system.time({
    liks_local <- foreach(i=1:1000,.packages='pomp',.combine=rbind) %dopar% {
      evals <- replicate(100, logLik(pfilter(seir_pomp,params=coef(mifs_local[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=900242057,kind="L'Ecuyer")

results_local <- data.frame(logLik=liks_local[,1],logLik_se=liks_local[,2],t(sapply(mifs_local,coef)))
summary(results_local$logLik,digits=7)
#pairs(~logLik+beta_I+beta_W, data=results_local)
pairs(~logLik+tau1+k+beta0, data=results_local)



