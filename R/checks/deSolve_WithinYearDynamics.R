##September 22, 2014
##The purpose of this file is to create a "map" of the ratio of BiomassIn/BiomassOut for each species
##as a function of BiomassIn_Sp1, BiomassIn_Sp2,R0
##This should give us a better grasp of the within season dynamics, so that we can ignore them and look at the 
##between season dynamics where the action is

library(deSolve)
nonsta=0
nsp = 4  
set.seed(1)
nyrs <- 1

source("sourcefiles/getRunParms.R") #define runtime parameters
source("sourcefiles/getGraphParms.R")  #define graphics parameters
source("sourcefiles/getEnvt.R")  #get constant and time-varying envt parms
source("sourcefiles/getSpecies.R")  #get species characteristics and Rstar

N0 <- 10
B0<- b*g[1,]*N0
R0 <- 2
R<-R0
B<-B0
State <- c(R=R0,B=B0)
Time <- seq(0,2,by=dt)  
Pars <- c(c=c,a=a,u=u,m=m,theta=theta,eps=eps)


## Within-season dynamics set-up- ODE
#source("sourcefiles/ResCompN.R") # define within-season ode solver AND defines Pars
#**************
ResCompN <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)),{
    dR = -eps*R - sum(B*a*R^theta/(1+a*u*R^theta))
    dB = (c*a*R^theta/(1+a*u*R^theta)-m)*B
    return(list(c(dR,dB)))
  })
}

##add rootfunction so that the integration stops once (opt1) R is numerically 0 (threshold hh)
rootfun <- function(Time, State, Pars) State[1] -min(Rstar)
#   dstate <- unlist(ResCompN(Time, State,Pars))
#                    dstate[2] - 1
# }
eventfun <- function(Time, State, Pars) State
events <- list(func = eventfun, root = TRUE, terminalroot = 1)
#***********

Brun <- as.data.frame(ode(func = ResCompN, y = State, parms = Pars, times = Time, rootfun=rootfun))


#Plot dynamics of a single year
# plot(Brun$B1~Brun$time,type="l",col="red",xlim = c(0,2),ylim = c(0,max(Brun)),ylab="Biomass",xlab="Time")
# points(Brun$B2~Brun$time,type="l", col="Blue")
# par(new=TRUE)
# plot(Brun$R~Brun$time,type="l",axes=FALSE,bty="n",xlab="",ylab="")
# axis(side=4,at=pretty(range(Brun$R)))
# mtext("Resource",side=4,line=3)
# 
