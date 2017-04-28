######################
#### POMP 
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


#START OF POMP

#Cummulative incidence

# Define the SEIR Step
seir_step <- Csnippet("
                      double dN_SE = rbinom(S,1-exp(beta/N*S*I));
                      double dN_EI = rbinom(E,1-exp(sigma*E));
                      double dN_IR = rbinom(I,1-exp(gamma*I));
                      double dN_RD = rbinom(R,1-exp(f*gamma*I));
                      S -= dN_SE;
                      E += dN_SE - dN_EI;
                      I += dN_EI - dN_IR;
                      R += dN_IR - dN_RD;
                      D += dN_RD;
                      C += dN_EI;
                      ")
# Initialize
seir_init <- Csnippet("
                      S = N-1;
                      E = 0;
                      I = 1;
                      R = 0;
                      D = 0;
                      C = 0;
                      ")

# Create the Pomp object
# May have to add another state to keep track of the transition in question
pomp_ebola <- pomp(data = ebola,time= "00",t0=00, 
                   rprocess=euler.sim(step.fun = seir_step, delta.t = 00), initializer = sir_init,
                   statenames = c("S","E","I","R"), paramnames = c("beta","gamma", "sigma"))

# Run the simulation
simStates <- simulate(pomp_ebola, nsim = 10,params = c(beta = 00, gamma=00, sigma=00), states=TRUE)

#Finding the negative log likelihood
#simulation of the observation process given the states and parameters (rmeasure)
#evaluation of the likelihood of a set of observations given the states and parameters (dmeasure)

#00 is the state var that keeps track of the transition in question
dmeas <- Csnippet("lik = dbinom(B,00,rho,give_log);")
rmeas <- Csnippet("B = rbinom(00,rho);")

pomp_ebola <- pomp(pomp_ebola,rmeasure=rmeas,dmeasure=dmeas,statenames="00",paramnames="rho")

sim <- simulate(pomp_ebola,params=c())

#END OF POMP