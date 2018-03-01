#archive

#Cummulative incidence

# Define the SEIR Step
seir_step <- Csnippet("
                      double dN_SE = rbinom(S,1-exp(beta1/N*S*I));
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
N <- 1e6
seir_init <- Csnippet("
                      S = N-1;
                      E = 0;
                      I = 1;
                      R = 0;
                      D = 0;
                      C = 0;
                      ")

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


# Create the Pomp object
# May have to add another state to keep track of the transition in question
#Working from http://kingaa.github.io/short-course/parest/parest.html

pomp_ebola <- pomp(data = data,
                   time="times",
                   t0=0, 
                   skeleton = vectorfield(
                     Csnippet("
                              DS = -Beta*S*I/N;
                              DE = Beta/N*S*I - sigma*E;
                              DI = sigma*E - gamma*I;
                              DR = (1-ff)*gamma*I;
                              DD = ff*gamma*I;
                              DC = sigma*E;
                              ")
                     ),
                   initializer = Csnippet("
                                          S = S_0;
                                          E = E_0;
                                          I = I_0;
                                          R = N-S_0-I_0;
                                          "),
                   statenames = c("S","E","I","R", "C", "D"), 
                   paramnames = c("Beta", "sigma", "gamma", "S_0", "E_0","I_0","R_0","ff", "N"))

# Sum of the squared errors 
sse <- function (params) {
  x <- trajectory(pomp_ebola,params=params)
  discrep <- x["I",,]-obs(pomp_ebola)
  sum(discrep^2)
}

f1 <- function (beta) {
  params <- c(Beta=beta,gamma=1/7.411374,sigma=1/9.312799,ff=qlogis(49/69),N=1e6,S_0=N-1,I_0=1,E_0=0,R_0=0)
  sse(params)
}
beta <- seq(from=0,to=3,by=0.5)
SSE <- sapply(beta,f1)
beta.hat <- beta[which.min(SSE)]
plot(beta,SSE,type='l')
abline(v=beta.hat,lty=2)


# rprocess=euler.sim(step.fun = seir_step, delta.t = 1/365),
# Run the simulation
simStates <- simulate(pomp_ebola, nsim = 10,params = c(beta = 0.3706847, gamma=1/7.411374, sigma=1/9.312799), states=TRUE)

# Stochastic Simulations
# MIF2: calculating likelihoods at 
# Iterated filtering 2 Pseudo


#Finding the negative log likelihood
#simulation of the observation process given the states and parameters (rmeasure)
#evaluation of the likelihood of a set of observations given the states and parameters (dmeasure)

#00 is the state var that keeps track of the transition in question
dmeas <- Csnippet("lik = dbinom(B,00,rho,give_log);")
rmeas <- Csnippet("B = rbinom(00,rho);")

pomp_ebola <- pomp(pomp_ebola,rmeasure=rmeas,dmeasure=dmeas,statenames="00",paramnames="rho")

sim <- simulate(pomp_ebola,params=c())

#END OF POMP