############################
## superspreading + intervention model code
## Creates various functions that will be put into the POMP Model
############################


pomp_skel <- Csnippet('
                      double N;
                      double beta, x;
                      if (t < tau1) {
                        beta = beta0;
                      } else{
                        x = -k*(t-tau1);
                        beta = (beta0)*(exp(x));
                      }
                      
                      N = S+E+Ih+Il+R+D;
                      // Balance the equations
                      DS = -beta/N*S*Ih;
                      DE = beta/N*S*Ih - sigma*E;
                      DIl = (1 - p0) * sigma*E - gamma*Il;
                      DIh = p0 * sigma*E - gamma*Ih;
                      DR = (1-ff)*gamma*Ih + (1-ff)*gamma*Il;
                      DD = ff*gamma*Ih + ff*gamma*Il;
                      DC = sigma*E;
                      ')

pomp_rprocess <- Csnippet('double beta, x, N;
                          if (t < tau1) {
                            beta = beta0;
                          } else{
                            x = -k*(t-tau1);
                            beta = (beta0)*(exp(x));
                          }
                          
                          N = S+E+Il+Ih+R+D;
                          //Rprintf(\"%lg %lg \\n\",p0, beta0);
                          double new_exposed = rbinom(S, 1 - exp(-beta/N*Ih*dt));
                          
                          double new_infected = rbinom(E, 1 - exp(-sigma*dt));
                          double new_infected_h = rbinom(new_infected, p0);
                          double new_infected_l = new_infected - new_infected_h;
                          
                          double new_not_infectious_h = rbinom(Ih, 1 - exp(-gamma*dt));
                          double new_not_infectious_l = rbinom(Il, 1 - exp(-gamma*dt));
                          
                          double new_recovered_h = rbinom(new_not_infectious_h, (1 - ff) );
                          double new_recovered_l = rbinom(new_not_infectious_l, (1 - ff) );
                          double new_dead = (new_not_infectious_h - new_recovered_h) + (new_not_infectious_l - new_recovered_l);
                          
                          S -= new_exposed;
                          E += new_exposed - new_infected;
                          Ih += new_infected_h - new_not_infectious_h;
                          Il += new_infected_l - new_not_infectious_l;
                          R += new_recovered_h + new_recovered_l;
                          D += new_dead;
                          C += new_infected;
                          ')


pomp_dmeasure <- Csnippet('double f;
                          f = dpois(nearbyint(cases), C,1);
                          lik = (give_log) ? f : exp(f);
                          ')

pomp_rmeasure <- Csnippet('
                          cases = rpois(C);
                          ')

pomp_trans <- Csnippet('
                  Tk = log(k);
                  Tbeta0 = log(beta0);
                  Tp0 = qlogis(p0, 0, 1, 1, 0);
                  ')

pomp_untrans <- Csnippet('
                    Tbeta0 = exp(beta0);
                    Tp0 = plogis(p0,0,1,1,0);
                    Tk = exp(k);
                    ')
pomp_init <- Csnippet("
                 S = 999998;
                 E = 0.0;
                 Il = 1.0;
                 Ih = 1.0;
                 R = 0.0;
                 D = 0.0;
                 C = 1.0;
                 ")
pomp_parms <- c(
  sigma = 1/9.312799, 
  gamma = 1/7.411374, 
  ff = plogis(49/69),
  beta0 = 2,
  k = .2,
  p0 = 0.05,
  tau1 = 15)

pomp_state_names <- c("S","E","Il", "Ih", "R", "D", "C")
pomp_param_names <- c("beta0","p0", "sigma", "gamma", "ff", "tau1", "k")

generate_pomp_model <- function (outbreak=c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende"),
                            data = NULL) {
  
  data <- pull_outbreak_data(outbreak, data)
  
  pomp(data = data,
       times="times",
       t0=-1,
       skeleton = vectorfield(pomp_skel), 
       params = pomp_parms,
       initializer = pomp_init,
       statenames = pomp_state_names, 
       paramnames = pomp_param_names,
       zeronames = c("C"),
       fromEstimationScale = pomp_untrans,
       toEstimationScale = pomp_trans,
       dmeasure = pomp_dmeasure,
       rmeasure = pomp_rmeasure,
       rprocess = euler.sim(step.fun = pomp_rprocess, delta.t = 1))
}

