######################
#### POMP Stochastic Model of DRC Ebola Outbreak 
#### Stochastic Model made by following Section 5 of
### http://kingaa.github.io/pomp/vignettes/pompjss.pdf
######################

# init
rm(list = ls())

#Libraries
library(deSolve)
library(chron)
library(bbmle)
library(mvtnorm)
library(pomp)
library(lubridate)
library(tidyverse)
library(cowplot)

#random number generator
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
cases <- data$cases

#START OF POMP

rmeas <- "
  cases = rnbinom_mu(theta,rho*H);
"

dmeas <- "
  lik = dnbinom_mu(cases, theta, rho*H, give_log);
"
### H is the true incidence vs the cases which was the measured incidence
### rnbinom_mu: negative binomial simulator (R function)
### dnbinom_mu: density function



# Define the SEIR Step

### Using the rate list may be an easier way to transition to the
### different beta value after 14 iterations

### reulermultinom and deulermultinom parameterize the Euler-multinomial 
### distributions by the size, rates, and time interval.


### SEIR step still under construction, I last left off trying to adapt the example from page 33
### of the guide http://kingaa.github.io/pomp/vignettes/pompjss.pdf
### In the guide dN[2] and dN[4] are referred to in the model, but are not set to equal anything.
### I believe that they are death rates, but I was not certain.

seir.step <- "
    double rate[5];
    double dN[5];
    double P;
    P = S + E + I + R;
    rate[0] = beta/P*S*I; // Exposure
    rate[1] = sigma*E; // Transition from exposed to infected
    rate[2] = gamma*I; // recovery from I
    rate[3] = (1-mu)*gamma*I; // recovery from I minus death from I
    rate[4] = mu*gamma*I; // death from I

    reulermultinom(2, S, &rate[0], dt, &dN[0]);
    reulermultinom(2, E, &rate[1], dt, &dN[1]);
    reulermultinom(2, I, &rate[2], dt, &dN[2]);
    reulermultinom(1, R, &rate[4], dt, &dN[4]);
    S -= dN[0];
    E += dN[0] - dN[1];
    I += dN[1] - dN[2];
    R += dN[3];
    D += dN[4];
    H += dN[1];
    "

### I believe reported times were daily, if not then delta.t needs to be reset to 1/52 for weeks
### Theta and rho have been set to zero until I figure out what to set them to
### Beta was set to zero while I figure out what to do about the changing beta value
### I believe that rho is the reporting efficiency, but I am not sure about that
### mu was used instead of the Althouse variable f, because POMP uses a variable called 'f'



seir <- pomp(data = data,
                times = "times", t0 = 0, 
                dmeasure = Csnippet(dmeas),
                rmeasure = Csnippet(rmeas), 
                rprocess = euler.sim(step.fun = Csnippet(seir.step), 
                delta.t = 1/365),
                statenames = c("S", "E", "I", "R", "H", "D"),
                paramnames = c("gamma", "mu", "theta", "sigma","beta", "rho", "S.0", "E.0" ,"I.0", "R.0"), 
                zeronames = "H", 
                initializer = function(params, t0, ...) {
                fracs <- params[c("S.0", "E.0", "I.0", "R.0")]
                setNames(c(round(params["N"] * fracs/sum(fracs)), 0),
                c("S","E", "I", "R", "H", "D"))}, 
                params = c(N = 1e6, beta = 0, sigma = 1/9.312799, gamma = 1/7.411374, mu = qlogis(49/69), rho = 0.1, 
                           theta = 100, S.0 = (1e6)-1, E.0 = 0, I.0 = 1.0, R.0 = 0,D.0 = 0))

### King: seed = 1914679908L
### Althouse: seed = 529342

seir <- simulate(seir, seed = 529342)


