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
                          
                          double new_exposed = rbinom(S, 1 - exp(-beta/N*S*I*dt));
                          double new_infected = rbinom(E, 1 - exp(-sigma*E*dt));
                          double new_not_infectious = rbinom(I, 1 - exp(-gamma*I*dt));       
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
seir_parm <- c(S.0= pop-1, 
               E.0 = 0.0,
               I.0 = 1.0,
               R.0 = 0.0,
               D.0 = 0.0,
               C.0 = 0.0,
               sigma = 1/9.312799, 
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
  statenames= c("S","E","I","R", "D", "C"), 
  paramnames=c("tau1", "beta0","k", "sigma", "gamma", "ff","beta1"),
  zeronames = c("C"),
  fromEstimationScale = untrans,
  toEstimationScale = trans,
  dmeasure = seir_dmeasure,
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
times <- 1:100
rnot <- if_else(times<t_match$params["tau1"], 
                t_match$params["beta0"]/t_match$params["gamma"], 
                exp(-t_match$params["k"]*(times - t_match$params["tau1"]))*t_match$params["beta0"]/t_match$params["gamma"])
plot(times, rnot, type = "l") # Plot matches althaus fig 2

## Do particle filtering of parameters
seir_parm["k"] = t_match$params["k"]
seir_parm["tau1"] = t_match$params["tau1"]
seir_parm["beta0"] = t_match$params["beta0"]

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



