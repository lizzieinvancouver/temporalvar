### Started 10 July 2011 ###
### By Lizzie & Megan ###

## VarEnvironments & Coexistence ##

# safety feature(s)
options(stringsAsFactors=FALSE)

# packages
library(ggplot2)

#set.seed(2)

# define all parameters
nsp <- 20    # number of spp
nyrs1 <- 100  # number of yrs to run first stationary period
ndays <- 1  # number of days in a growing season
dt <- 0.0005 # within yr timestep
dtnc <- 10*dt
tsteps <- ndays/dt
# params for adding on a nonstationary run after stationary run
nyrs2 <- 50 # number of yrs for second run (nonstationary for now)
nyrs <- nyrs1+nyrs2 # ALERT! change below once not doing stationary+nonstationary run
y <- c(1:nyrs)

## Extinction Threshold:  1 seed per hectare (assuming that initial density is 10 seeds per meter)
ext <- 1/10000

# set up graphics parameters
colerz <- topo.colors(nsp)
rcol <- "darkslateblue"
linez <- rep(1:6, 100) # enough for 600 species for now
lspbyrs <- 1
lresbyrs <- 2
lwd=2

## Setting up loop for multiple model runs
# To do still:
# (1) add in everything crossyrsvars need (was lazy about this)
# (2) add in resource stuff to withinyrs
# (3) decide on list for each run, versus some other format

modelruns <- list() # place to store output of runs
nruns <- 2 # number of model runs to do
for (j in c(1:nruns)){ # assuming, we will vary species characteristics between yrs ... 
  
  ##
  ## species characteristics
  ##
  b <-  rep(1,nsp)          # biomass of seedling
  s <-  rep(0.8,nsp)      # seedbank survival overwinter
  a <-  rep(20,nsp)        # slope of species uptake rate with increasing R
  u <-  rep(1,nsp)          # inverse of the max uptake rate
  c <-  rep(12,nsp)        # conversion of resource to biomass
  m <-  rep(0.05,nsp)     # mortality
  gmax <-  rep(0.5,nsp)     # max germination fraction
  h <-  rep(100,nsp)             # max rate of germination decrease following pulse
  phi <- rep(0.05,nsp)     # conversion of end-of-season plant biomass to seeds
  tauI <- runif(nsp,0.1, 0.9)    # time of max germ for sp i
  # set up tracking (to do for Lizzie, there should be code under here someday)
  
  
  theta <- rep(1,nsp)         # shape of species i uptake curve, remember it should be an integer!
  N0 <- rep(10,nsp)          # initial number of seeds (per meter square?)
  Rstar <- (m/(a*(c-m*u)))^(1/theta)
  
  crossyrsvars <- as.data.frame(cbind(b, s, a, u, c, m, gmax, h, phi, theta, tauI, Rstar))
  
  ##
  ## time-varying env variables
  ##
  mu <- log(2)  #mean of resource distribution
  sigma <- 0.2  #sd of resource distribution
  R0 <- rlnorm(nyrs, mu, sigma) # intial R in a season
  eps <- 1              # evaporative stress
  #tauP <- 0.3           # timing of pulse
  p <- 2  #first parameter for beta distribution of tau
  q <- 2  #second parameter for beta distribution of tau
  tauPs <- rbeta(nyrs1, p, q) # change once not doing stationary+nonstationary run
  
  # nonstationary tauP, change once not doing stationary+nonstationary run
  qns <- seq(2, 20, length.out=nyrs2)
  tauPns <- rbeta(nyrs2, p, qns) # yes, it takes a vector! Yay!
  plot(tauPns~c(1:50))
  tauP <- c(tauPs, tauPns)
  
  ##
  ## Within-growing season dynamics set-up
  ##
  R <- matrix(rep(0), nyrs, ndays/dt) # R is the resource level through the growing season (each yr)
  B <- array(rep(0), dim=c(nyrs,nsp,tsteps)) # where B is an array with yr (nyrs), spp biomass
  # through growing season (ndays)
  N <- matrix(rep(0), nyrs, nsp) # number of seeds by yr and spp
  N[1,] <- N0  #initialize
  Bfin <- matrix(rep(0),nyrs,nsp) # biomass at end of year y
  
  ##
  ## set-up for different coexistence mechanisms
  ##
  BnoC <- array(rep(0), dim=c(nyrs,nsp,tsteps)) # where B is an array with yr (nyrs), spp biomass
  RnoC <- array(rep(0), dim=c(nyrs,nsp,tsteps)) # where B is an array with yr (nyrs), spp biomass
  BnoCfin <- matrix(rep(0),nyrs,nsp) # biomass at end of year y
  tauIstar <- matrix(rep(0),nyrs,nsp)
#  Bnocomp <- matrix(rep(0),nyrs,nsp) # B without competition at end of year y
  E <- matrix(rep(0),nyrs,nsp)
  C <- matrix(rep(0),nyrs,nsp)
  
  ##
  ## change to mapply someday?
  ## for now, a loop
  ## Better to use ODE solver within each year?
  ##
  
  for (y in c(1:(nyrs-1))){
    g <- gmax*exp(-h*(tauP[y]-tauI)^2)  #germination fraction in year y
    k<-1
    R[y,k] <- R0[y]
    B[y,,k] <- b*g*N[y,]
    while (R[y,k]>min(Rstar)){
      f <- (a*R[y,k]^theta)/(1+a*u*R[y,k]^theta)
      B[y,,k+1] <- B[y,,k]+(c*f-m)*B[y,,k]*dt
      R[y,k+1] <- R[y, k] -dt*(t(B[y,,k]) %*% f + eps*R[y,k])
      R[y,k+1] <- R[y,k+1]*(R[y,k+1]>0)
      k <- k+1
    }
    Bfin[y,] <- apply(B[y,,], 1, max)  #final biomass
    # add some internal calculations to make other calculations easier
    tauIstar[y,] <- (log(R0[y])/eps)-(1/(theta*eps))*log(m/(a*c-a*u*m))
    Bnocomp[y,] <- B[y,,1]*((1+a*u*R0[y]^theta)/(1+a*u*R0[y]^theta*exp(-eps*tauIstar[y,]*theta)))*
      exp((-c/(u*eps*theta))-m*tauIstar[y,])
    # E[y,] <- log(s*g*(phi*Bnocomp[y,]-1))
    # C[y,] <- log((phi*Bnocomp[y,]-1)/(phi*Bfin[y,]-1))
    N[y+1,] <- s*(N[y,]*(1-g)+phi*Bfin[y,])  #convert biomass to seeds and overwinter
    N[y+1,] <- N[y+1,]*(N[y+1,]>ext)  #if density does not exceed ext, set to zero
  }
 
  for (y in c(1:(nyrs-1))){
    for (s in c(1:nsp)){
      g <- gmax[s]*exp(-h[s]*(tauP[y]-tauI[s])^2)  #germination fraction in year y
      k<-1
      RnoC[y,s,k] <- R0[y]
      BnoC[y,s,k] <- b[s]*g*N[y,s]
      while (RnoC[y,s,k]>min(Rstar[s])){
        f <- (a[s]*RnoC[y,s,k]^theta[s])/(1+a[s]*u[s]*RnoC[y,s,k]^theta[s])
        BnoC[y,s,k+1] <- BnoC[y,s,k]+(c[s]*f-m[s])*BnoC[y,s,k]*dtnc
        RnoC[y,s,k+1] <- RnoC[y,s,k] -dtnc*(t(BnoC[y,s,k]) %*% f + eps*RnoC[y,s,k])
        RnoC[y,s,k+1] <- RnoC[y,s,k+1]*(RnoC[y,s,k+1]>0)
        k <- k+1
      }

    BnoCfin[y,] <- max(B[y,s,])  #max biomass
    }
 }
  
  
  modelruns[[j]] <- list(crossyrsvars, Bfin, E)
  # could also make each run a multi-part dataframe with common names
  # so something like:
  # modelruns[[paste("crossyrs", j, sep="")]] <- crossyrsvars
  # modelruns[[paste("withinyrs", j, sep="")]] <- Bfin
}


# between years plot
dev.new(width=14, height=10)
plot(Bfin[,1]~c(1:nyrs), ylim=c(min(Bfin), max(Bfin)),
     xlab="year", ylab="Abundance", type="n")
for (i in 1:nsp) {
  lines(Bfin[,i]~c(1:nyrs), col=colerz[i], lty=lspbyrs, lwd=lwd)
}

dev.new(width=14, height=10)
range=c(1:(tsteps/3))
plot(R[1,range]~range, ylim=c(min(R), max(R)),
     xlab="step, step, step", ylab="Resource", type="n")
for (i in 2:nyrs) {
  lines(R[i,range]~range, col=colerz[i], lty=lspbyrs, lwd=lwd)
}


# within years plots
# tweaked a little to make biomass differences among years clearer
dev.new(width=14, height=10)
par(mfrow=c(3,3))
par(mar=c(5, 4, 4, 8) + 0.1)
selectyrs <- seq(1, nyrs, by=floor(nyrs/8))
yrlength <- ndays/dt
for (yr in seq_along(selectyrs)){
  allsp <- B[selectyrs[yr], ,]
  plot(B[selectyrs[yr], 1,]~c(1:yrlength), ylim=c(min(allsp), max(B)), 
       xlab="step, step, step",  ylab="Abundance", type="n",
       main=paste("year: ", selectyrs[yr], sep=""))
  for (sp in c(1:nsp)){
    lines(B[selectyrs[yr], sp ,]~c(1:yrlength), ylim=c(min(allsp), max(B)),
          col=colerz[sp])
  }
  par(new=TRUE)
  plot(R[selectyrs[yr],]~c(1:yrlength), axes=FALSE, xlab="", ylab="",
       ylim=c(min(R[selectyrs[yr],]), max(R[selectyrs[yr],])), type="l",
       col=rcol, lty=lresbyrs, lwd=lwd)
  raxis <- seq(0, max(R[selectyrs[yr],]), by=max(R[selectyrs[yr],])/10)
  axis(4, at=raxis, labels=round(raxis, digits=2))
  mtext("Resource", side=4, line=3, cex=0.75)
}



# two options for overlay histograms
par(mfrow=c(1,1))

dev.new(width=5, height=6)
# (1) using the base package (I hate this, let's rm it)
maxhist <- max(tauI, tauP)*1.1
tauIfin <- tauI[which(Bfin[max(y),]>0)]
tauIlosers <- tauI[which(Bfin[max(y),]<=0)]
hist(tauIfin, col=rgb(1, 0, 0, 0.5), xlim=c(0, maxhist), ylim=c(0, nsp),
     main=paste(sum(Bfin[max(y),]>0), "out of", nsp, "coexisted", sep=" "), xlab="taus")
hist(tauP, col=rgb(0, 0, 1, 0.5), add=TRUE)
hist(tauIlosers, col=rgb(0, 0, 1, 0.5), add=TRUE)

dev.new(width=7, height=6)
# (2) using ggplot, which really is good for this sort of thing
tau.df <- data.frame(coexisted=Bfin[max(y),]>0, tauI=tauI)
tau.df$coexisted[tau.df$coexisted==TRUE] <- "coexisted"
tau.df$coexisted[tau.df$coexisted==FALSE] <- "doomed"
tauP.df <- data.frame(coexisted=rep("tauP"), tauP=tauP)

ggplot(tau.df, aes(tauI, fill = coexisted)) + geom_histogram(alpha=0.5) +
  geom_density(data=tauP.df[1:100,], aes(tauP),  alpha = 0.2) +
  geom_density(data=tauP.df[101:150,], aes(tauP),  alpha = 0.4) +
  labs(title=paste(sum(Bfin[max(y),]>0), "out of", nsp, "coexisted", sep=" "))


## notes for lizzie (by lizzie):
# length of vectors is nsp
# while loop runs through season until minimum R* is met
# while loop runs until SECOND species is below its R*
# but we should run it so at R* each species converts to seeds (maybe?)
# for now we convert max biomass to seeds

# y is years
# t is within years

# %% is mod (integer division: it's the remainder once a even division)

print(c("The number of coexisting species are",sum(Bfin[max(y),]>0),"out of",nsp))
