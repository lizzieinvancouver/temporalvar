#define function for ODE solver
#R is resource
#B is vector of nsp consumers

RstarComp <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)),{
    dR = -eps*R - sum(B*a*R^theta/(1+a*u*R^theta))
    dB = (c*a*R^theta/(1+a*u*R^theta)-m)*B
    return(list(c(dR,dB)))
  })
}

##add root function so that the integration stops once (opt1) R is numerically 0 (threshold hh)
rootfun <- function(Time, State, Pars) State[1] - min(Rstar)


