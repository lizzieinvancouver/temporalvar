### Started 10 July 2011 ###
### By Lizzie & Megan ###

## VarEnvironments & Coexistence ##

# define all parameters

nsp <- 2    # number of spp
nyrs <- 100  # number of yrs
ndays <- 10  # number of days in a growing season
dt <- 0.001 # within yr timestep
y <- c(1:nyrs)
tsteps <- ndays/dt
t <-c(1:tsteps)


##
## species characteristics
##
b <-  c(1,1)          # biomass of seedling
s <-  c(0.8,0.8)      # seedbank survival overwinter
a <-  c(20,20)        # slope of species uptake rate with increasing R
d <-  c(1,1)          # inverse of the max uptake rate
c <-  c(12,12)        # conversion of resource to biomass
m <-  c(0.05,.05)     # mortality
G <-  c(0.5, 0.5)     # max germination fraction
h <-  100             # max rate of germination decrease following pulse
phi <- c(0.05,0.05)     # conversion of end-of-season plant biomass to seeds
tauI <- c(0.3, 0.4)    # time of max germ for sp i
theta <- c(1,1)         # shape of species i uptake curve
N0 <- c(10,10)          # initial number of seeds
Rstar <- (m/(a*(c-m*d)))^(1/theta)

##
## time-varying env variables
##
R0 <- rlnorm(nyrs, ln(2), 0.2) # intial R in a season
eps <- 1              # evaporative stress
#tauP <- 0.3           # timing of pulse
p <- 2  #first parameter for beta distribution of tau
q <- 2  #second parameter for beta distribution of tau
##
## Within-growing season dynamics set-up
##
R <- matrix(rep(0), nyrs, ndays/dt) # R is the resource level through the growing season (each yr)
B <- array(rep(0), dim=c(nyrs,nsp,ndays/dt)) # where B is an array with yr (nyrs), spp biomass
    # through growing season (ndays)
N <- matrix(rep(0), nyrs, nsp) # number of seeds by yr and spp
Bfin <- matrix(rep(0),nyrs,nsp) # biomass at end of year y
##
## change to mapply someday?
## for now, a loop
## Better to use ODE solver within each year?
##

for (y in c(1:nyrs)){
  tauP <- rbeta(1,p,q)
  g <- G*exp(-h*(tauP-tauI)^2)
  t <- 1
  R[y,t] <- R0[y]
  print(y)
  if(y==1) N[y,] <- N0
  else N[y,] <- s*(N[y-1,]*(1-g)+phi*Bfin[y-1,])
  B[y,, t] <- b*g*N[y,]
  f <- (a*R[y,1]^theta)/(1+a*d*R[y,1]^theta)
  print(f)
  while (R[y,t]>min(Rstar)){
      t <- t+1
      f <- (a*R[y,t-1]^theta)/(1+a*d*R[y,t-1]^theta)
      if((sum(c*f<m)>0)) print(paste("c*f-m=", c*f-m,", t=",t)) #what is this for?  
      B[y,,t] <- B[y,,t-1]+(c*f-m)*B[y,,t-1]*dt
      R[y,t] <- R[y, t-1]-dt*(t(B[y,,t-1]) %*% f - eps*R[y,t-1])
      #J <- which(R[y, t-1]-Rstar < 0)
      #  B[y,,t][J] <- 0
    }
  print(y)
  Bfin[y,] <- apply(B[y,,], 1, max)
}
plot(Bfin[,1]~c(1:nyrs),xlab="year",ylab="Abundance",type="l")
points(Bfin[,2]~c(1:nyrs),col='blue',type="l")

