############################
## superspreading POMP Model Code
## Creates various functions that will be put into the POMP Model
############################


seir_skel <- Csnippet('
                      double N;
                      
                      N = S+E+I+R+D;
                      // Balance the equations
                      DS = -beta0/N*S*I;
                      DE = beta0/N*S*I - sigma*E;
                      DI = sigma*E - gamma*I;
                      DR = (1-ff)*gamma*I;
                      DD = ff*gamma*I;
                      DC = sigma*E;
                      ')

seir_rprocess <- Csnippet('double beta, N;
            beta = beta0;
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
                          //Rprintf(\"%lg \\n\",f);
                          lik = (give_log) ? f : exp(f);
                          ')

seir_rmeasure <- Csnippet('
                          cases = rpois(C);
                          ')

trans <- Csnippet('
                  Tbeta0 = log(beta0);
                  ')

untrans <- Csnippet('
                    Tbeta0 = exp(beta0);
                    ')

ebola_seir_model <- function (outbreak=c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende"),
                            data = NULL, sim=FALSE) {
  
  # populations <- c(Yambuku=275000,Kikwit=200000,Mweka2007=170000,Mweka2008=170000,Isiro=700000,Boende=250000)
  outbrk <- match.arg(outbreak)
  # pop <- unname(populations[outbrk])
  
  # browser()
  data %>%
    ungroup() %>%
    filter(outbreak==outbrk) %>%
    select(times,cases) -> data
  
  data <- as.data.frame(data)
  
  init <- Csnippet("
                   S = 999999;
                   E = 0.0;
                   I = 1.0;
                   R = 0.0;
                   D = 0.0;
                   C = 1.0;
                   ")
  if (sim==FALSE) {
  seir_parm <- c(
    sigma = 1/9.312799, 
    gamma = 1/7.411374, 
    ff = plogis(49/69),
    beta0 = 2
    )} else {
  seir_parm <- c(
    sigma = 1/9.312799, 
    gamma = 1/7.411374, 
    ff = plogis(49/69),
    beta0 = .1281814
  )}
  
  names_seir <- c("S","E","I", "R", "D", "C")
  
  para_seir <- c("beta0", "sigma", "gamma", "ff")
  
  
  pomp(data = data,
       times="times",
       t0= -1,
       skeleton = vectorfield(seir_skel), 
       params = seir_parm,
       initializer=init,
       statenames= names_seir, 
       paramnames=para_seir,
       zeronames = c("C"),
       fromEstimationScale = untrans,
       toEstimationScale = trans,
       dmeasure = seir_dmeasure,
       rmeasure = seir_rmeasure,
       rprocess = euler.sim(step.fun = seir_rprocess, delta.t = 1)) -> ss_seir_pomp
  
  
  return(ss_seir_pomp)
}

