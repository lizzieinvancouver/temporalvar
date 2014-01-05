#define function for ODE solver
#R is resource
#B is vector of consumers
NoCompN <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)),{
    dR = -eps*R 
    dB = (c*a*R^theta/(1+a*u*R^theta)-m)*B
    return(list(c(dR,dB)))
  })
}
Pars <- c(c=c,a=a,u=u,m=m,theta=theta,eps=eps)
