#define function for ODE solver
#R is resource
#B is vector of nsp consumers

RstarComp <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)),{
    
    uptake1 <- a[1] * R^theta[1] / (1 + a[1] * u[1] * R^theta[1])
    uptake2 <- a[2] * R^theta[2] / (1 + a[2] * u[2] * R^theta[2])
    
    dR  <- -eps * R - B1 * uptake1 - B2 * uptake2 
    dB1 <- (c[1] * uptake1 - m[1]) * B1
    dB2 <- (c[2] * uptake2 - m[2]) * B2
    
    return(list(c(dR,dB1,dB2)))
    # dR = -eps*R - sum(B*a*R^theta/(1+a*u*R^theta))
    # dB = (c*a*R^theta/(1+a*u*R^theta)-m)*B
    # return(list(c(dR,dB)))
  })
}

##add root function so that the integration stops once (opt1) R is numerically 0 (threshold hh)
rootfun <- function(Time, State, Pars) State[1] - min(Rstar)

#Define all parameters in the model
Pars <- c(c=c,a=a,u=u,m=m,theta=theta,eps=eps)

#Names for events
Bnames <- c("B1","B2")
