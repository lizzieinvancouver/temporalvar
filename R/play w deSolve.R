#deSolve vignette

library(deSolve)

LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*(alpha - beta*y)
    dy = -y*(gamma - delta*x)
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
Pars <- c(c1=12,c2=12,a1=20,a2=20,u1=1, u2=1,m1=.05,m2=.05,theta1=1,theta2=1,eps=1)
State <- c(N1=10,N2=20,R=log(2))
Time <- seq(0,1000,by=.1)

out <- as.data.frame(ode(func = ResComp2, y = State, parms = Pars, times = Time))

matplot(out[,-1], type = "l", xlab = "time", ylab = "population")


#try N consumers, 1 resource
nsp = 12
ResCompN <- function(Time,Inits,Pars) {
  with(as.list(c(State,Pars)),{
    dN = (c*a*R^theta/(1+a*u*R^theta)-m)*N
    dR  = -eps*R - sum(a*R^theta/(1+a*u*R^theta))
    return(list(c(dN,dR)))
  })
}
Pars <- c(c=rep(12,nsp),a=rep(20,nsp),u=rep(1,nsp),m=rep(0.05,nsp),theta=rep(1,nsp),eps=1)
Inits <- c(N=seq(1,12,by=1),R=100)
Time <- seq(0,25,by=.1)

out <- as.data.frame(ode(func = ResCompN, y = Inits, parms = Pars, times = Time))

matplot(out[,-1], type = "l", xlab = "time", ylab = "population")

##ooooh! it works!