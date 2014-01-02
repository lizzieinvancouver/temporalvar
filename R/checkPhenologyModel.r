### Started 10 December 2013 ###
### By Lizzie & Megan ###

## Checking assumptions of PhenologyModel.R ##
## Also some code below to check coding etc.. ##

options(stringsAsFactors=FALSE)

## First! Does species with highest Rstar win when E is constant? ##
## Yes she does, go you, you low Rstar species! ##
## How checked? I changed the below ##

nyrs <- 5000 # Run it for a while
tauI <- rep(0.6, nsp) # Make sure all spp parameters are the same, except ...
theta[5] <- 0.8 # for one, make one spp the winner in an Rstar model

# Next, make the environment the same:
R0 <- rep(R0[1], length(R0)) 
tauP <- rep(0.5, nyrs)

# couple things to look at:
N[5000,]
Rstar
colMeans(Bfin) #means

##
## Check that with R0 constant and just tauP varying
## species with tauI closest to tauP wins

nyrs <- 1000 
ndays <- 20 # I know, crazy
# tauI <- runif(nsp,0.1, 0.9)    # time of max germ for sp i
tauI <- rep(0.7, nsp)
tauI[5] <- 0.5
R0 <- rep(R0[1], length(R0))

# appears to be headed in right direction, though we did not make it to extinction
N[1000,]
E[999,]

##
## Check that E covaries with Rstar ##
##

ndays <- 5  # number of days in a growing season

# tauI <- runif(nsp,0.1, 0.9)    # time of max germ for sp i
tauI <- rep(0.7, nsp)
theta <- runif(nsp,0.1, 0.9)

tauP <- 0.3           # timing of pulse
p <- 2  #first parameter for beta distribution of tau
q <- 2  #second parameter for beta distribution of tau
# tauP <- rbeta(nyrs,p,q) 
R0 <- rlnorm(nyrs, mu, sigma) # intial R in a season
R0 <- rep(R0[1], length(R0)) 

# check code
order(colMeans(E))
order(Rstar) # Yay!

##
## model run here (currently this is the one that checks Rstar and E)
##

# define all parameters
nsp <- 20    # number of spp
nyrs <- 100  # number of yrs 
ndays <- 5  # number of days in a growing season
dt <- 0.001 # within yr timestep
y <- c(1:nyrs)
tsteps <- ndays/dt

## Extinction Threshold:
# 1 seed per hectare (assuming that initial density is 10 seeds per meter)
ext <- 1/10000

##
## species characteristics
##
b <-  rep(1,nsp)          # biomass of seedling
s <-  rep(0.8,nsp)      # seedbank survival overwinter
a <-  rep(20,nsp)        # slope of species uptake rate with increasing R
d <-  rep(1,nsp)          # inverse of the max uptake rate
c <-  rep(12,nsp)        # conversion of resource to biomass
m <-  rep(0.05,nsp)     # mortality
G <-  rep(0.5,nsp)     # max germination fraction
h <-  rep(100,nsp)             # max rate of germination decrease following pulse
phi <- rep(0.05,nsp)     # conversion of end-of-season plant biomass to seeds
# tauI <- runif(nsp,0.1, 0.9)    # time of max germ for sp i
tauI <- rep(0.7, nsp)

theta <- runif(nsp,0.1, 0.9)  # shape of species i uptake curve
N0 <- rep(10,nsp)          # initial number of seeds (per meter square?)
Rstar <- (m/(a*(c-m*d)))^(1/theta)

##
## time-varying env variables
##
mu <- log(2)  #mean of resource distribution
sigma <- 0.2  #sd of resource distribution
R0 <- rlnorm(nyrs, mu, sigma) # intial R in a season
R0 <- rep(R0[1], length(R0)) 
eps <- 1              # evaporative stress
# tauP <- 0.3           # timing of pulse
p <- 2  #first parameter for beta distribution of tau
q <- 2  #second parameter for beta distribution of tau
# tauP <- rbeta(nyrs,p,q)
tauP <- rep(0.3, nyrs)

##
## Within-growing season dynamics set-up
##
R <- matrix(rep(0), nyrs, ndays/dt) # R is the resource level through the growing season (each yr)
B <- array(rep(0), dim=c(nyrs,nsp,tsteps)) # where B is an array with yr (nyrs), spp biomass
    # through growing season (ndays)
N <- matrix(rep(0), nyrs, nsp) # number of seeds by yr and spp
N[1,] <- N0  #initialize
B0 <- matrix(rep(0),nyrs,nsp) # biomass at beginning of year y
Bfin <- matrix(rep(0),nyrs,nsp) # biomass at end of year y
BnoC <- matrix(rep(0),nyrs,nsp) # biomass at end of year y
rcrt <- matrix(rep(0),nyrs,nsp) # recruitment in year y
E <- matrix(rep(0),nyrs,nsp)

##
## change to mapply someday?
## for now, a loop
## Better to use ODE solver within each year?
##

for (y in c(1:(nyrs-1))){
  g <- G*exp(-h*(tauP[y]-tauI)^2)  #germination fraction in year y
  k<-1  #k is the within year step
  R[y,] <- R0[y]
  B0[y,] <- b*g*N[y,]
  B[y,,k] <- B0[y,]
  while (R[y,k]>min(Rstar)){
      f <- (a*R[y,k]^theta)/(1+a*d*R[y,k]^theta)
      B[y,,k+1] <- B[y,,k]+(c*f-m)*B[y,,k]*dt
      R[y,k+1] <- R[y, k] -dt*(t(B[y,,k]) %*% f + eps*R[y,k])
      R[y,k+1] <- R[y,k+1]*(R[y,k+1]>0)
      k <- k+1
    }
  #final biomass: max biomass for each species; assumes seed set for sp j at the point that R<R*(sp=j)
  Bfin[y,] <- apply(B[y,,], 1, max)  
  tstar <- (1/eps)*(log(R0[y]-(1/theta)*log(m/(a*c-a*d*m))))
  e1 <- 1+a*d*R0[y]^theta
  e2 <- 1+a*d*R0[y]^theta*exp(-eps*tstar*theta)
  e3 <- c/d/eps/theta
  BnoC[y,] <- B0[y,] * e1^(-e3) * e2^(-e3) * exp(-m*tstar)
  E[y,] <- s*g*(phi*BnoC[y,]-1)      #the -1 accounts for the loss of adults to germination
  rcrt[y,] <- s*g*(phi*Bfin[y,]-1)   #to recruit, convert biomass to seeds and overwinter
  C[y,] <- E[y,]/rcrt[y,]
  N[y+1,] <- s*(N[y,]) + rcrt[y,]    #survival + recruitment
  N[y+1,] <- N[y+1,]*(N[y+1,]>ext)   #if density does not exceed ext, set to zero
}




stop(print("Shifting focus here..."))

## Currently, some double coding to check work and a simulation to work
# on some issues ##


# double-coding to try to avoid errors
#tauIstar <- (log(R0[y])/eps)-1/(theta*eps)*log(m/(a*c-a*u*m))
#Bnocomp <- B[y,,1]*((1+a*u*R0[y]^theta)/(1+a*u*(R0[y]^theta)*exp(-eps*tauIstar*theta)))*exp((-c/(u#*eps*theta))-m*tauIstar)
tstar <- (1/eps)*(log(R0[y]-(1/theta)*log(m/(a*c-a*d*m))))
e1 <- 1+a*d*R0[y]^theta
e2 <- 1+a*d*R0[y]^theta*exp(-eps*tstar*theta)
e3 <- c/d/eps/theta
BnoC[y,] <- B0[y,] * e1^(-e3) * e2^(-e3) * exp(-m*tstar)

## simulation without competition

ndays <- 10
tsteps <- ndays/dt

R <- matrix(rep(0), nyrs, ndays/dt) # R is the resource level through the growing season (each yr)
B <- array(rep(0), dim=c(nyrs,nsp,tsteps)) # where B is an array with yr (nyrs), spp biomass
    # through growing season (ndays)
N <- matrix(rep(0), nyrs, nsp) # number of seeds by yr and spp
N[1,] <- N0  #initialize
Bfin <- matrix(rep(0),nyrs,nsp) # biomass at end of year y
BnoC <- matrix(rep(0),nyrs,nsp) # B without competition at end of year y
##

for (y in c(1:(nyrs-1))){
  g <- gmax*exp(-h*(tauP[y]-tauI)^2)  #germination fraction in year y
  k<-1
  R[y,k] <- R0[y]
  B[y,,k] <- b*g*N[y,]
  while (R[y,k]>min(Rstar)){
      f <- (a*R[y,k]^theta)/(1+a*u*R[y,k]^theta)
      B[y,,k+1] <- B[y,,k]+(c*f-m)*B[y,,k]*dt
      R[y,k+1] <- R[y, k] -dt*(eps*R[y,k])
      R[y,k+1] <- R[y,k+1]*(R[y,k+1]>0)
      k <- k+1
    }
  Bfin[y,] <- apply(B[y,,], 1, max)  #final biomass
  # add some internal calculations to make other calculations easier
#  ts <- (log(R0[y])/eps)-(1/(theta*eps))*log(m/(a*c-a*u*m))
#  Bnocomp[y,] <- B[y,,1]*((1+a*u*R0[y]^theta)/(1+a*u*R0[y]^theta*
#      exp(-eps*tauIstar*theta)))*exp((-c/(u*eps*theta))-m*tauIstar) # ack! Why so small?
  tstar <- (1/eps)*(log(R0[y]-(1/theta)*log(m/(a*c-a*d*m))))
  e1 <- 1+a*d*R0[y]^theta
  e2 <- 1+a*d*R0[y]^theta*exp(-eps*tstar*theta)
  e3 <- c/d/eps/theta
  BnoC[y,] <- B0[y,] * e1^(-e3) * e2^(-e3) * exp(-m*tstar)
  E[y,] <- s*g*(phi*BnoC[y,]-1)      #the -1 accounts for the loss of adults to germination
  rcrt[y,] <- s*g*(phi*Bfin[y,]-1)   #to recruit, convert biomass to seeds and overwinter
  C[y,] <- E[y,]/rcrt[y,]
  N[y+1,] <- s*(N[y,]) + rcrt[y,]    #survival + recruitment
  #E[y,] <- log(s*g*(phi*Bnocomp[y,]-1))
  #C[y,] <- log((phi*Bnocomp[y,]-1)/(phi*Bfin[y,]-1))
  #N[y+1,] <- s*(N[y,]*(1-g)+phi*Bfin[y,])  #convert biomass to seeds and overwinter
  N[y+1,] <- N[y+1,]*(N[y+1,]>ext)  #if density does not exceed ext, set to zero
}
