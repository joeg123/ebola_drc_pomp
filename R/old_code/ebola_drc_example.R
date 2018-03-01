#############################
## Testing SIR model
#############################
rm(list=ls())

library(pomp)
library(tidyverse)
library(cowplot)
library(lubridate)

ebola_data <- read_csv(file = "~/Desktop/elife-09015-supp1-v1-download.csv")
boende_cases <- ebola_data %>% filter(Outbreak == "Boende") %>%
  mutate(Date_of_onset_symp = parse_date(Date_of_onset_symp, format = "%m/%d/%y")) %>%
  group_by(Date_of_onset_symp) %>%
  summarize(cases = n()) %>% rename(date = Date_of_onset_symp) %>% ungroup() %>%
  mutate(day = as.numeric(date - min(date) + 1))

############################################################################################

seir_rprocess <- Csnippet('
                 double lambda, beta;
                 beta = R0 * gamma; // Transmission rate
                 lambda = beta * I / N; // Force of infection
                 
                 // Transitions
                 // From class S
                 double transS = rbinom(S, 1.0 - exp(- lambda * dt)); // No of infections
                 // From class E
                 double transE = rbinom(E, 1.0 - exp(-E * alpha * dt));
                 // From class I
                 double transI = rbinom(I, 1.0 - exp(-gamma * dt)); // No of transitions I->R
                 
                 // Balance the equations
                 S -= transS;
                 E += transS - transE;

                 I += transE - transI;
                 R += transI;
                 N_EI += transE; // No of transitions from E to I
                 ')

seir_dmeasure <- Csnippet('
                 double f;
                 f = dpois(nearbyint(cases), N_EI,1);
                 lik = (give_log) ? f : exp(f);
                 ')

seir_rmeasure <- Csnippet('
                  cases = rpois(N_EI);
                 ')


seir_skel <- Csnippet('
                      double lambda, beta;
                      beta = R0 * gamma; // Transmission rate
                      lambda = beta * I / N; // Force of infection
                      
                      // Balance the equations
                      DS = - lambda * S;
                      DE = lambda * S - alpha * E;
                      DI = alpha * E - gamma * I;
                      DR = gamma * I;
                      DN_EI = alpha * E;
                      ')

boende_cases %>% 
  select(cases:day) %>%
  pomp(
    times="day",
    t0=0,
    skeleton = vectorfield(seir_skel),
    rprocess = discrete.time.sim(seir_rprocess, delta.t = 0.1),
    rmeasure = seir_rmeasure,
    dmeasure = seir_dmeasure,
    statenames = c("S", "E", "I", "R", "N_EI"),
    paramnames = c("N","R0","alpha","gamma"),
    zeronames = c("N_EI")) -> mod


parms <- parmat(c(N=100000, R0 = 2, alpha = 1/9.312799, gamma = 1/7.411374, S.0 = 99999, E.0 = 0, I.0 = 1, R.0 =0, N_EI.0=0))

temp <- simulate(mod, params=parms, nsim = 10, states=TRUE, as.data.frame=TRUE, times=seq(1,500,by=0.1))
temp %>% ggplot(aes(time, N_EI, color=sim))+ geom_line()

temp <- trajectory(mod,params=parms,times=seq(1,500,by=0.1), as.data.frame=TRUE)

temp %>% gather(state, value, S:R) %>%
  ggplot(aes(time, value, color=state)) + geom_line()


fit <- traj.match(mod, start = c(N=100000, R0=1.1, alpha = 1/9.312799, gamma = 1/7.411374, S.0 = 99999, E.0 = 0, I.0 = 1, R.0 =0, N_EI.0=0),
           est = c("R0"))
summary(fit)

simulate(fit, nsim = 100, as.data.frame=TRUE) %>%
  ggplot(aes(time, N_EI, group = as.factor(sim)))+
  geom_line() 








set.seed(594709947L)
library(ggplot2)
theme_set(theme_bw())
library(plyr)
library(reshape2)
library(magrittr)
library(pomp)
stopifnot(packageVersion("pomp")>="1.6")


base_url <- "http://kingaa.github.io/sbied/"
read.csv(paste0(base_url,"ebola/ebola_data.csv"),stringsAsFactors=FALSE,
         colClasses=c(date="Date")) -> dat
sapply(dat,class)
head(dat)

populations <- c(Guinea=10628972,Liberia=4092310,SierraLeone=6190280)

rSim <- Csnippet('
                 double lambda, beta;
                 double *E = &E1;
                 beta = R0 * gamma; // Transmission rate
                 lambda = beta * I / N; // Force of infection
                 int i;
                 
                 // Transitions
                 // From class S
                 double transS = rbinom(S, 1.0 - exp(- lambda * dt)); // No of infections
                 // From class E
                 double transE[nstageE]; // No of transitions between classes E
                 for(i = 0; i < nstageE; i++){
                 transE[i] = rbinom(E[i], 1.0 - exp(-nstageE * alpha * dt));
                 }
                 // From class I
                 double transI = rbinom(I, 1.0 - exp(-gamma * dt)); // No of transitions I->R
                 
                 // Balance the equations
                 S -= transS;
                 E[0] += transS - transE[0];
                 for(i=1; i < nstageE; i++) {
                 E[i] += transE[i-1] - transE[i];
                 }
                 I += transE[nstageE-1] - transI;
                 R += transI;
                 N_EI += transE[nstageE-1]; // No of transitions from E to I
                 N_IR += transI; // No of transitions from I to R
                 ')

rInit <- Csnippet("
                  double m = N/(S_0+E_0+I_0+R_0);
                  double *E = &E1;
                  int j;
                  S = nearbyint(m*S_0);
                  for (j = 0; j < nstageE; j++) E[j] = nearbyint(m*E_0/nstageE);
                  I = nearbyint(m*I_0);
                  R = nearbyint(m*R_0);
                  N_EI = 0;
                  N_IR = 0;
                  ")

#' 
## ----skel,include=FALSE--------------------------------------------------
skel <- Csnippet('
                 double lambda, beta;
                 const double *E = &E1;
                 double *DE = &DE1;
                 beta = R0 * gamma; // Transmission rate
                 lambda = beta * I / N; // Force of infection
                 int i;
                 
                 // Balance the equations
                 DS = - lambda * S;
                 DE[0] = lambda * S - nstageE * alpha * E[0];
                 for (i=1; i < nstageE; i++)
                 DE[i] = nstageE * alpha * (E[i-1]-E[i]);
                 DI = nstageE * alpha * E[nstageE-1] - gamma * I;
                 DR = gamma * I;
                 DN_EI = nstageE * alpha * E[nstageE-1];
                 DN_IR = gamma * I;
                 ')

#' 
## ----measmodel,include=FALSE---------------------------------------------
dObs <- Csnippet('
                 double f;
                 if (k > 0.0)
                 f = dnbinom_mu(nearbyint(cases),1.0/k,rho*N_EI,1);
                 else
                 f = dpois(nearbyint(cases),rho*N_EI,1);
                 lik = (give_log) ? f : exp(f);
                 ')

rObs <- Csnippet('
                 if (k > 0) {
                 cases = rnbinom_mu(1.0/k,rho*N_EI);
                 } else {
                 cases = rpois(rho*N_EI);
                 }')

#' 
## ----partrans,include=FALSE----------------------------------------------
toEst <- Csnippet('
                  const double *IC = &S_0;
                  double *TIC = &TS_0;
                  TR0 = log(R0);
                  Trho = logit(rho);
                  Tk = log(k);
                  to_log_barycentric(TIC,IC,4);
                  ')

fromEst <- Csnippet('
                    const double *IC = &S_0;
                    double *TIC = &TS_0;
                    TR0 = exp(R0);
                    Trho = expit(rho);
                    Tk = exp(k);
                    from_log_barycentric(TIC,IC,4);
                    ')

#' 
## ----pomp-construction,include=FALSE-------------------------------------
ebolaModel <- function (country=c("Guinea", "SierraLeone", "Liberia"),
                        timestep = 0.1, nstageE = 3) {
  
  ctry <- match.arg(country)
  pop <- unname(populations[ctry])
  nstageE <- as.integer(nstageE)
  
  globs <- paste0("static int nstageE = ",nstageE,";")
  
  dat <- subset(dat,country==ctry,select=-country)
  
  incubation_period <- 11.4/7
  infectious_period <- 7/7
  index_case <- 10/pop
  dt <- timestep
  theta <- c(N=pop,R0=1.4,
             alpha=-1/(nstageE*dt)*log(1-nstageE*dt/incubation_period),
             gamma=-log(1-dt/infectious_period)/dt,
             rho=0.2,
             k=0,
             S_0=1-index_case,E_0=index_case/2-5e-9,
             I_0=index_case/2-5e-9,R_0=1e-8)
  
  ## Create the pomp object
  dat %>% 
    extract(c("week","cases")) %>%
    pomp(
      times="week",
      t0=min(dat$week)-1,
      globals=globs,
      statenames=c("S",sprintf("E%1d",seq_len(nstageE)),
                   "I","R","N_EI","N_IR"),
      zeronames=c("N_EI","N_IR"),
      params=theta,
      paramnames=c("N","R0","alpha","gamma","rho","k",
                   "S_0","E_0","I_0","R_0"),
      dmeasure=dObs, rmeasure=rObs,
      rprocess=discrete.time.sim(step.fun=rSim, delta.t=timestep),
      skeleton=vectorfield(skel),
      toEstimationScale=toEst,
      fromEstimationScale=fromEst,
      initializer=rInit) -> po
}

ebolaModel("Guinea") -> gin

ebola_pf <- pfilter(gin, Np = 5000)

logLik(ebola_pf)

sims <- simulate(gin, nsim=20,as.data.frame=TRUE,include.data=TRUE)

ggplot(sims,mapping=aes(x=time,y=cases,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)



## Calculating the parameter estimates
library(foreach)
library(doParallel)
registerDoParallel()

library(doRNG)
registerDoRNG(625904618)
bake(file="guinea/pf.rds",{
  foreach(i=1:10,.packages='pomp',
          .export=c("gin")
  ) %dopar% {
    pfilter(gin,Np=10000)
  }
}) -> pf
(L_pf <- logmeanexp(sapply(pf,logLik),se=TRUE))


results <- as.data.frame(as.list(c(coef(pf[[1]]),loglik=L_pf[1],loglik=L_pf[2])))
write_csv(results,path="guinea/ebola_params.csv")

registerDoRNG(482947940)
bake(file="box_search_local.rds",{
  foreach(i=1:20,
          .packages='pomp',
          .combine=c, 
          .export=c("gin")
  ) %dopar%  
  {
    mif2(
      gin,
      Np=2000,
      Nmif=50,
      cooling.type="geometric",
      cooling.fraction.50=0.5,
      transform=TRUE,
      rw.sd=rw.sd(R0 = .02, alpha = .02, gamma = 0.02)
    ) -> mifs_local
  }
}) -> mifs_local

############################################################################################