############################
## Althaus POMP Model Code
## Creates various functions that will be put into the POMP Model
############################

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
                  Tk = log(k);
                  Tff = qlogis(ff,0,1,1,0);
                  Ttau1 = log(tau1);
                  ')

untrans <- Csnippet('
                    Tbeta0 = exp(beta0);
                    Tk = exp(k);
                    Tff = plogis(ff,0,1,1,0);
                    Ttau1 = exp(tau1);
                    ')


generate_pomp_model <- function (outbreak=c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende", "Equator"),
                                 data = NULL) {
  pop <- 1e6   
  
  if(outbreak == "Equator"){
    init <- Csnippet("
                   S= 999998; 
                   E = 0.0;
                   I = 2.0;
                   R = 0.0;
                   D = 0.0;
                   C = 2.0;
                   ")
  } else{
    init <- Csnippet("
                   S= 999999; 
                   E = 0.0;
                   I = 1.0;
                   R = 0.0;
                   D = 0.0;
                   C = 1.0;
                   ")  
  }
  
  
  seir.init.state <- c(S.0= 999999, 
                       E.0 = 0.0,
                       I.0 = 1.0,
                       R.0 = 0.0,
                       D.0 = 0.0,
                       C.0 = 1.0)
  
  seir_parm <- c(sigma = 1/9.312799, 
                 gamma = 1/7.411374, 
                 ff = plogis(49/69),
                 tau1 = 14.0,
                 beta0 = 0.4,
                 k = 1e-4)
  
  pomp_statenames= c("S","E","I","R", "D", "C")
  
  pomp_paramnames=c("tau1", "beta0","k", "sigma", "gamma", "ff")
  
    
  data <- pull_outbreak_data(outbreak,data)
  
  pomp(data = data,
       times = "times",
       t0 = -1,
       skeleton = vectorfield(seir_skel), 
       params = seir_parm,
       initializer=init,
       statenames=pomp_statenames, 
       paramnames=pomp_paramnames,
       zeronames = c("C"),
       fromEstimationScale = untrans,
       toEstimationScale = trans,
       dmeasure = seir_dmeasure,
       rmeasure = seir_rmeasure,
       rprocess = euler.sim(step.fun = seir_rprocess, delta.t = 1))
}
