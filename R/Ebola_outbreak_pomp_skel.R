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
                      //beta = beta1 + (beta0-beta1)*(exp(x));
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

seir_dmeasure <- Csnippet('double f;
                          //if (C == ){
                            //Rprintf(\" %lg \", C);
                            //f = dpois(nearbyint(cases), I,1);
                            //lik = (give_log) ? f : exp(f);
                          //} else {
                          
                          //Rprintf(\" I = %lg \", I);
                          f = dpois(nearbyint(cases), C,1);
                          //if (f < 0) f = NA_REAL;
                          //Rprintf(\" F = %lg \", f);
                          lik = (give_log) ? f : exp(f);
                          //}
                          Rprintf(\" C = %lg \", C);
                          ')
trans <- Csnippet('
                    Tbeta0 = exp(beta0);
                    Tbeta1 = exp(beta1);
                    Tk = exp(k);
                    Tff = qlogis(ff,0,1,1,0);
                    //Ttau0 = exp(tau0);
                  ')

untrans <- Csnippet('
                    Tbeta0 = log(beta0);
                    Tbeta1 = log(beta1);
                    Tk = log(k);
                    Tff = plogis(ff,0,1,1,0);
                    //Ttau0 = log(tau0);
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
               ff = qlogis(49/69),
               tau1 = 14.0,
               beta0 = 0.4,
               beta1 = 1,
               k = log(1))

#Althouse Results
#!beta0!      beta1          !k!          f       tau0       !tau1!      sigma      gamma 
#-0.3634729       -Inf -2.0883390  0.8960880       -Inf 14.2890158  0.1073791  0.1349277

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
  dmeasure = seir_dmeasure) -> seir_pomp

# zeronames = c("C"),
traj <- trajectory(seir_pomp, as.data.frame=TRUE)

#traj %>% gather(state, value, C) %>%
#  ggplot(aes(time, value, color=state)) + geom_line()

#trajectory matching

#est will probably not use all of those args
traj.match(seir_pomp, 
           start=seir_parm,
           method = c("Nelder-Mead"),
           est=c("tau1", "beta0","k"),
           transform=TRUE) -> t_match

## Its the untransform?
## Transforms twice and then it untransforms and stops due to initial params, when untrans is first it immediately stops
## The else {beta = ...} is the problem!! after 14 iters of program t becomes greater than tau0 and then terrible things happen
plot(t_match)
summary(t_match)
rm(t_match)

# I
#tau1        beta0            k 
#1.512e+01   3.706847e-01     2.103075e-01 

