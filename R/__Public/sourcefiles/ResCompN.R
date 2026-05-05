## To do ...
# This is where we have to stagger the starts ...
# Megan says, smart thing to do would be to:
    # (1) calculate R until min(\tau_{g,i} days (evapotranspiration period) 
    # (2) start ODE with first germinating species at \tau{g,i} for first germinating species with remaining R
    # (3) then re-start ODE \tau{g,j} (second germinating species) with R remaining and species 1 biomass as initial conditions
# Megan shoud do this as it seems triky and annoying to Lizzie

#define function for ODE solver
#R is resource
#B is vector of consumers
ResCompN <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)),{
    dR = -eps*R - sum(B*a*R^theta/(1+a*u*R^theta))
    dB = (c*a*R^theta/(1+a*u*R^theta)-m)*B
    return(list(c(dR,dB)))
  })
}
Pars <- c(c=c,a=a,u=u,m=m,theta=theta,eps=eps)

##add rootfunction so that the integration stops once (opt1) R is numerically 0 (threshold hh)
rootfun <- function(Time, State, Pars) State[1] -min(Rstar)

