######################
#### POMP Stochastic Model of DRC Ebola Outbreak
######################

# init
rm(list = ls())

# Necessary libraries
library(deSolve)
library(chron)
library(bbmle)
library(mvtnorm)
library(pomp)
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

#START OF POMP

# Define the SEIR Step
seir_step <- Csnippet("
                      double dN_SE = rbinom(S,1-exp(Beta*I/N*dt));
                      double dN_EI = rbinom(E,1-exp(sigma*dt));
                      double dN_IR = rbinom(I,1-exp(gamma*dt));
                      S -= dN_SE;
                      E += dN_SE - dN_EI;
                      I += dN_EI - dN_IR;
                      R += dN_IR;
                      H += dN_EI;
                      ")
# User guide H += IR
#D += dN_RD;
#C += dN_EI;
#double dN_RD = rbinom(R,1-exp(f*gamma*I));


# Initialize
N <- 1e6
seir_init <- Csnippet("
                      S = (1e6)-1;
                      E = 0;
                      I = 1;
                      R = 0;
                      H = 0;
                      ")

seir <- pomp(data = data,
          time="times",
          t0=0,
          rprocess=euler.sim(seir_step,delta.t=1/365),
          initializer=seir_init,paramnames=c("sigma","Beta","gamma","N"),
          statenames=c("S","E", "I","R","H"))

# rho is the probability that new cases result in confinement
# H tracks the transition from E to I, and is set to zero after each iteration

pomp(seir,zeronames="H") -> seir
rmeas <- Csnippet("B = rbinom(H,rho);")
seir <- pomp(seir,rmeasure=rmeas,statenames="H",paramnames="rho")

sims <- simulate(seir,params=c(Beta=0.4,gamma=1/7.411374,sigma = 1/9.312799,rho=0.9,N=1e6),
                 nsim=20,as.data.frame=TRUE,include.data=TRUE)

summary(seir)
plot(seir)

ggplot(sims,mapping=aes(x=time,y=B,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)

