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
R0 <- rlnorm(nyrs, dlnorm(2), 0.2) # intial R in a season
## alert, lizzie changed ln above because my computer was confused by it, log-normal command okay?
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
## set up graphics parameters
##
colerz <- topo.colors(nsp)
rcol <- "darkslateblue"
linez <- rep(1:6, 100) # enough for 600 species for now
lspbyrs <- 1
lresbyrs <- 2
lwd=2

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
      if((sum(c*f<m)>0)) print(paste("c*f-m=", c*f-m,", t=",t)) #c*f<m is asking whether any species is below its R*
      B[y,,t] <- B[y,,t-1]+(c*f-m)*B[y,,t-1]*dt
      R[y,t] <- R[y, t-1]-dt*(t(B[y,,t-1]) %*% f - eps*R[y,t-1])
      #J <- which(R[y, t-1]-Rstar < 0)
      #  B[y,,t][J] <- 0
    }
  #print(y)
  Bfin[y,] <- apply(B[y,,], 1, max)
} 
  #Megan's lame plotting.  It would be nice to call a plot function that showed within year increase in biomass of all species on the left axis
  # and within year decrease in resource on the right axis
  if (y%%10 == 0) {
    dev.new(width=14, height=10)
    plot(B[y,1,]~c(1:tsteps), ylim=c(min(B[y,,]), max(B[y,,])),
        xlab="year", ylab="Abundance", type="n")
    lines(B[y,1,]~c(1:tsteps), col=colerz[1], lty=lspbyrs, lwd=lwd)
    lines(B[y,2,]~c(1:tsteps), col=colerz[2], lty=lspbyrs, lwd=lwd)
    lines(R[y,]~c(1:tsteps), col=rcol, lty=lresbyrs, lwd=lwd)
    legend(nyrs, (max(B[y,,])-(0.1*max(B[y,,]))), "resource",
       col=c(rcol), lty=c(lresbyrs), lwd=1, bty="n")
  }
}

# between years plot
dev.new(width=14, height=10)
plot(Bfin[,1]~c(1:nyrs), ylim=c(min(Bfin), max(Bfin)),
    xlab="year", ylab="Abundance", type="n")
lines(Bfin[,1]~c(1:nyrs), col=colerz[1], lty=lspbyrs, lwd=lwd)
lines(Bfin[,2]~c(1:nyrs), col=colerz[2], lty=lspbyrs, lwd=lwd)

# within years plots
# B[is.na(B)] <- 0 
dev.new(width=14, height=10)
par(mfrow=c(3,3))
par(mar=c(5, 4, 4, 8) + 0.1)
selectyrs <- seq(1, nyrs, by=floor(nyrs/8))
yrlength <- ndays/dt
for (yr in seq_along(selectyrs)){
    allsp <- B[selectyrs[yr], ,]
    plot(B[selectyrs[yr], 1,]~c(1:yrlength), ylim=c(min(allsp), max(allsp)), 
        xlab="step, step, step",  ylab="Abundance", type="n",
        main=paste("year: ", selectyrs[yr], sep=""))
    par(new=TRUE)
    plot(R[selectyrs[yr],]~c(1:yrlength), axes=FALSE, xlab="", ylab="",
        ylim=c(min(R[selectyrs[yr],]), max(R[selectyrs[yr],])), type="l",
        col=rcol, lty=lresbyrs, lwd=lwd)
    raxis <- seq(0, max(R[selectyrs[yr],]), by=max(R[selectyrs[yr],])/10)
    axis(4, at=raxis, labels=round(raxis, digits=2))
    mtext("Resource", side=4, line=3, cex=0.75)
    for (sp in c(1:nsp)){
        lines(B[selectyrs[yr], sp ,]~c(1:yrlength), ylim=c(min(allsp), max(allsp)),
            col=colerz[sp])
  }
}
    

## notes for lizzie (by lizzie):
# length of vectors is nsp
# while loop runs through season until minimum R* is met
# while loop runs until SECOND species is below its R*
# but we should run it so at R* each species converts to seeds (maybe?)
# for now we convert max biomass to seeds

# y is years
# t is within years

# %% is mod (integer division: it's the remainder once a even division)
