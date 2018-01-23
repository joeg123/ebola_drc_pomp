


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