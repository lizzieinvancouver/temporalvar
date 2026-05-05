#define function for ODE solver
#R is resource
#B is vector of nsp consumers
# Megan found this promising stackoverflow 
# https://stackoverflow.com/questions/69843063/how-to-solve-a-system-of-ode-with-time-dependent-parameters-in-r
# https://tpetzoldt.github.io/deSolve-forcing/deSolve-forcing.html


RstarComp <- function(Time,State,Pars, Germinate1, Germinate2) {
  g1 <- Germinate1(Time)
  g2 <- Germinate2(Time)
  with(as.list(c(State,Pars)),{
    uptake1 <- a[1] * R^theta[1] / (1 + a[1] * u[1] * R^theta[1])
    uptake2 <- a[2] * R^theta[2] / (1 + a[2] * u[2] * R^theta[2])
    
    dR  <- -eps * R - B1 * uptake1 - B2 * uptake2 
    dB1 <- (c[1] * uptake1 - m[1]) * B1 + (N0[1] * s[1] * b[1])*g1
    dB2 <- (c[2] * uptake2 - m[2]) * B2 + (N0[2] * s[2] * b[2])*g2
    
    return(list(c(dR,dB1,dB2),g = c(g1,g2)))
    # dR = -eps*R - sum(B*a*R^theta/(1+a*u*R^theta))
    # dB = (c*a*R^theta/(1+a*u*R^theta)-m)*B
    # return(list(c(dR,dB)))
  })
}

#Germination is a forcing function with new values each year
#  g1,g2 spread the germination over the growing season
#  version 1 puts all germination on day tau_g50, gmax[yr,sp]
#  version 2 spreads germination by Hill Eqn with g50 at tau_g50

