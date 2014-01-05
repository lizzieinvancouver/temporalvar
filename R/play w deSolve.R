#deSolve vignette

library(deSolve)

LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx <- x*(alpha - beta*y)
    dy <- -y*(gamma - delta*x)
    return(list(c(dx, dy)))
  })
}

Pars <- c(alpha = 2, beta = .5, gamma = .2, delta = .6)
State <- c(x = 10, y = 10)
Time <- seq(0, 100, by = 1)

out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))

matplot(out[,-1], type = "l", xlab = "time", ylab = "population")
legend("topright", c("Cute bunnies", "Rabid foxes"), lty = c(1,2), col = c(1,2), box.lwd = 0)

#2 consumers, 1 resource

ResComp2 <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)),{
    dN1 = (c1*a1*R^theta1/(1+a1*u1*R^theta1)-m1)*N1
    dN2 = (c2*a2*R^theta2/(1+a2*u2*R^theta2)-m2)*N2
    dR  = -eps*R - N1*a1*R^theta1/(1+a1*u1*R^theta1) - N2*a2*R^theta2/(1+a2*u2*R^theta2)
    return(list(c(dN1,dN2,dR)))
  })
}
Pars <- c(c1=12,c2=12,a1=10,a2=20,u1=1, u2=1,m1=.05,m2=.05,theta1=1,theta2=1,eps=1)
State <- c(N1=10,N2=20,R=log(2))
Time <- seq(0,1000,by=.1)

outResComp2 <- as.data.frame(ode(func = ResComp2, y = State, parms = Pars, times = Time))

matplot(outResComp2[,-1], type = "l", xlab = "time", ylab = "population")


#try N consumers, 1 resource
nsp = 20
ndays=3
dt=.0001
eps=1
source("getSpecies.R")
source("ResCompN.R")
# ResCompN <- function(Time,State,Pars) {
#   with(as.list(c(State,Pars)),{
#     dR  = -eps*R - sum(B*a*R^theta/(1+a*u*R^theta))
#     dB = (c*a*R^theta/(1+a*u*R^theta)-m)*B
#     return(list(c(dR,dB)))
#   })
# }
#Pars <- c(c=c,a=a,u=u,m=m,theta=theta,eps=eps)
R0=20
B0=c(rep(0,10), seq(5,10,length.out=10))
B=B0
State <- c(R=R0,B=B)
#Time <- seq(0,3,by=.1)
outResCompN <- as.data.frame(ode(func = ResCompN, y = State, parms = Pars, times = Time))
matplot(Time, outResCompN[,-1],type = "l", xlab = "time", ylab = "population")

#No competition model
NoCompN <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)),{
    dN = (c*a*R^theta/(1+a*u*R^theta)-m)*N
    dR  = -eps*R
    return(list(c(dN,dR)))
  })
}
source("getSpecies.R")
eps=1
R0=2
tstar <- (1/eps)*(log(R0)-(1/theta)*log(m/(a*c-a*u*m)))
Pars <- c(c=c,a=a,u=u,m=m,theta=theta,eps=eps)
N=seq(1,20,1)
State <- c(N=N,R=R0)
Time <- seq(0,max(tstar)+.01,by=.01)
outNoCompN <- as.data.frame(ode(func = NoCompN, y = State, parms = Pars, times = Time))
matplot(Time, out[,-1],type = "l", xlab = "time", ylab = "population")
##ooooh! it works!