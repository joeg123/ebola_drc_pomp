############################
## superspreading POMP Model Code
## Creates various functions that will be put into the POMP Model
############################


seir_skel <- Csnippet('
                      double N;
                      
                      N = S+E+Ih+Il+R+D;
                      // Balance the equations
                      DS = -beta0/N*S*Ih;
                      DE = beta0/N*S*Ih - sigma*E;
                      DIl = (1 - p0) * sigma*E - gamma*Il;
                      DIh = p0 * sigma*E - gamma*Ih;
                      DR = (1-ff)*gamma*Ih + (1-ff)*gamma*Il;
                      DD = ff*gamma*Ih + ff*gamma*Il;
                      DC = sigma*E;
                      ')

seir_rprocess <- Csnippet('double N;
                          N = S+E+Il+Ih+R+D;
                          //Rprintf(\"%lg %lg \\n\",p0, beta0);
                          double new_exposed = rbinom(S, 1 - exp(-beta0/N*Ih*dt));

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


seir_dmeasure <- Csnippet('double f;
                          f = dpois(nearbyint(cases), C,1);
                          lik = (give_log) ? f : exp(f);
                          ')

seir_rmeasure <- Csnippet('
                          cases = rpois(C);
                          ')

trans <- Csnippet('
                  Tbeta0 = log(beta0);
                  Tp0 = qlogis(p0, 0, 1, 1, 0);
                  ')

untrans <- Csnippet('
                    Tbeta0 = exp(beta0);
                    Tp0 = plogis(p0,0,1,1,0);
                    ')

generate_pomp_model <- function (outbreak=c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende", "Equator"),
                        data = NULL, sim=FALSE) {

  data <- pull_outbreak_data(outbreak,data)

  if(outbreak == "Equator"){
    init <- Csnippet("
                 S = 999997;
                 E = 0.0;
                 Il = 1.0;
                 Ih = 2.0;
                 R = 0.0;
                 D = 0.0;
                 C = 3.0;
                 ")
  } else{
    init <- Csnippet("
                 S = 999998;
                 E = 0.0;
                 Il = 1.0;
                 Ih = 1.0;
                 R = 0.0;
                 D = 0.0;
                 C = 2.0;
                 ") 
  }
  


  if (sim==FALSE) {
    seir_parm <- c(
      sigma = 1/9.312799, 
      gamma = 1/7.411374, 
      ff = plogis(49/69),
      beta0 = 2,
      p0 = 0.05
  )} else {
    seir_parm <- c(
      sigma = 1/9.312799, 
      gamma = 1/7.411374, 
      ff = plogis(49/69),
      beta0 = 1.281814,
      p0 = .1
    )}

  names_seir <- c("S","E","Il", "Ih", "R", "D", "C")
  
  para_seir <- c("beta0","p0", "sigma", "gamma", "ff")
  
  
  pomp(data = data,
       times="times",
       t0=-1,
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

