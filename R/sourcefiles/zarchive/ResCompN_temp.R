#define function for ODE solver
#R is resource
#B is vector of consumers
ResCompN <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)),{
    Evap = eps*R
    Uptake =  B*a*R^theta/(1+a*u*R^theta)
    dR = -Evap - sum(Uptake)
    dB = (c*Uptake-m)*B
    return(list(c(dR,dB,Uptake)))
  })
}
Pars <- c(c=c,a=a,u=u,m=m,theta=theta,eps=eps)

##add rootfunction so that the integration stops once (opt1) R is numerically 0 (threshold hh) 
rootfun<-function(Time, State, Pars) {
  dstate <- unlist(ResCompN(Time, State, Pars))
  #opt 1. this works, makes threshold when R is 0
    hh <- abs(dstate[1]) - 1e-35  
  return(hh)
}
