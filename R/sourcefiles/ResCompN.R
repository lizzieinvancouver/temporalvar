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
rootfun <- function(Time, State, Pars) R - min(Rstar)
eventfun <- function(Time, State, Pars) State
events <- list(func = eventfun, root = TRUE, terminalroot = 1)

# rootfun<-function(Time, State, Pars) {
#  dstate <- unlist(ResCompN(Time, State, Pars))
  #opt 1. this works, makes threshold when R is 0
#    hh <- abs(dstate[1]) - 1e-35  
#  return(hh)
# }
