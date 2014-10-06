#This program is to compare the solution to the within season model using an ODE solver 
# versus stepping through the season in discrete steps.  Unsatifyingly, it does not give the same outcome.
# This is a concern, and we could check solutions for R, then R and one species of B, and then the full community.  
#  However, I don't have time to deal with this now, so we will go eith the solver.



library(deSolve)
nonsta=0
nsp = 2  #when nsp=2, tauI is assigned known values from chesson 2004
set.seed(1)

source("sourcefiles/getRunParms.R") #define runtime parameters
source("sourcefiles/getGraphParms.R")  #define graphics parameters
source("sourcefiles/getEnvt.R")  #get constant and time-varying envt parms
source("sourcefiles/getSpecies.R")  #get species characteristics and Rstar

N0 <- rep(10,nsp)          # initial number of seeds (per meter square?)
B0<- b*g[1,]*N0 
## Within-season dynamics set-up- ODE
source("sourcefiles/ResCompN_temp.R") # define within-season ode solver
source("sourcefiles/NoCompN.R")  # define within-season ode solver for no competition

## Within Season set up for STEP
ndays = 2.0
dt = 0.001
tss <- seq(1,ndays/dt)
Rss <- array(rep(0),length(tss)) # R is the resource level through the growing season
Bss <- matrix(rep(0), nrow=nsp,ncol=length(tss)) # where B is spp biomass within the year
#Bss <- array(rep(0), dim=c(ndays/dt)) # where B is ONE SP biomass within the year
Uss <- matrix(rep(0), nrow=nsp,ncol=length(tss)) # where USs is the uptake by each species at each timesetp

#step step
k<-1
Rss[1] <- R0[1]
Bss[,1] <- B0
#Bss[1]<- B0[1]
#while (Rss[k]>min(Rstar)){
for (k in seq(1,length(tss)-1)) {
  f <- (a*Rss[k]^theta)/(1+a*u*Rss[k]^theta)
  Bss[,k+1] <- Bss[,k]+(c*f-m)*Bss[,k]*dt
  Uss[,k+1] <- dt*f*Bss[,k]
  Rss[k+1] <- Rss[k] -dt*(t(Bss[,k]) %*% f + eps*Rss[k])
#  k <- k+1
}
#Rss[k+1] <- Rss[k] -dt*(eps*Rss[k])  #R only
#Rss[k+1] <- Rss[k+1]*(Rss[k+1]>0)

#Test for R only
#use deSolve for ResCompN
OnlyR <- function(Time,State,Pars) {
  with(as.list(c(State,Pars)),{
    dR = -eps*R
    return(list(dR))
  })
}


R<-R0[1]
B<-B0
U<-B0*a*R^theta/(1+a*u*R^theta)
#State<-c(R=R)  #Test for R only
#State<-c(R=R,B=B[1]) #Test for R and one B
#Pars <- c(c=c[1],a=a[1],u=u[1],m=m[1],theta=theta[1],eps=eps[1])
State<-c(R=R,B=B,U=U) 
Time <- seq(0,4,by=dt)
Bout <- as.data.frame(ode(func = ResCompN, y = State, parms = Pars, times = Time))
#Bout <- as.data.frame(ode(func = OnlyR, y = State, parms = Pars, times = Time))


#COmpare output
# plot(Bout$R~Bout$time,main="ODE", type="l",ylim = c(0,3),xlim = c(0,4))
# plot(Bout$B1~Bout$time,main="ODE",type="l",xlim = c(0,4))
# points(Bout$B2~Bout$time,type="l",main="ODE", col="Blue")
png("../figures/Biomass_R_time_ODE.png",width=6, height=4,units="in",res=400)
par(mar = c(5,4,4,4+.3))
plot(Bout$B1~Bout$time,main="ODE",type="l",col="red",xlim = c(0,4),ylab="Biomass",xlab="Time")
points(Bout$B2~Bout$time,type="l", col="Blue")
par(new=TRUE)
plot(Bout$R~Bout$time,type="l",axes=FALSE,bty="n",xlab="",ylab="")
axis(side=4,at=pretty(range(Bout$R)))
mtext("Resource",side=4,line=3)
mtext(Sys.time(),side=3,line=0)
dev.off()

png("../figures/Biomass_R_time_Step.png",width=6, height=4,units="in",res=400)
par(mar = c(5,4,4,4+.3))
plot(Bss[1,]~tss,main="Step-Step",type="l",col="red",xlim = c(0,40000),ylab="Biomass", xlab="Time")
points(Bss[2,]~tss,type="l", col="Blue")
par(new=TRUE)
plot(Rss~tss,type="l",axes=FALSE,bty="n",xlab="",ylab="")
axis(side=4,at=pretty(range(Rss)))
mtext("Resource",side=4,line=3)
mtext(Sys.time(),side=3,line=0)
dev.off()
