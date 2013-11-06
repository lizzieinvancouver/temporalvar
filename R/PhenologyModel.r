### Started 10 July 2011 ###
### By Lizzie & Megan ###

## VarEnvironments & Coexistence ##

#set.seed(2)

# define all parameters
nsp <- 20    # number of spp
nyrs <- 100  # number of yrs
ndays <- 2  # number of days in a growing season
dt <- 0.001 # within yr timestep
y <- c(1:nyrs)
tsteps <- ndays/dt


## Extinction Threshold:  1 seed per hectare (assuming that initial density is 10 seeds per meter)
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
tauI <- runif(nsp,0.1, 0.9)    # time of max germ for sp i

theta <- rep(1,nsp)         # shape of species i uptake curve
N0 <- rep(10,nsp)          # initial number of seeds (per meter square?)
Rstar <- (m/(a*(c-m*d)))^(1/theta)

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
tauP <- rbeta(nyrs,p,q)
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

for (y in c(1:(nyrs-1))){
  g <- G*exp(-h*(tauP[y]-tauI)^2)  #germination fraction in year y
  k<-1
  R[y,k] <- R0[y]
  B[y,,k] <- b*g*N[y,]
  while (R[y,k]>min(Rstar)){
      f <- (a*R[y,k]^theta)/(1+a*d*R[y,k]^theta)
      B[y,,k+1] <- B[y,,k]+(c*f-m)*B[y,,k]*dt
      R[y,k+1] <- R[y, k] -dt*(t(B[y,,k]) %*% f + eps*R[y,k])
      R[y,k+1] <- R[y,k+1]*(R[y,k+1]>0)
      k <- k+1
    }
  Bfin[y,] <- apply(B[y,,], 1, max)  #final biomass
  N[y+1,] <- s*(N[y,]*(1-g)+phi*Bfin[y,])  #convert biomass to seeds and overwinter
  N[y+1,] <- N[y+1,]*(N[y+1,]>ext)  #if density does not exceed ext, set to zero
} 

  #Megan's lame plotting.  It would be nice to call a plot function that showed within year increase in biomass of all species on the left axis
  # and within year decrease in resource on the right axis
#   if (y%%10 == 0) {
#     dev.new(width=14, height=10)
#     plot(B[y,1,]~c(1:tsteps), ylim=c(min(B[y,,]), max(B[y,,])),
#         xlab="year", ylab="Abundance", type="n")
#     #lines(B[y,1,]~c(1:tsteps), col=colerz[1], lty=lspbyrs, lwd=lwd)
#     #lines(B[y,2,]~c(1:tsteps), col=colerz[2], lty=lspbyrs, lwd=lwd)
#     lines(R[y,]~c(1:tsteps), col=rcol, lty=lresbyrs, lwd=lwd)
#     for (sp in c(1:nsp)){
#       lines(B[y,sp,]~c(1:tsteps), col=colerz[i], lty=lspbyrs, lwd=lwd)
#     }
#     legend(nyrs, (max(B[y,,])-(0.1*max(B[y,,]))), "resource",
#        col=c(rcol), lty=c(lresbyrs), lwd=1, bty="n")
#   }
# }

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

print(c("The number of coexisting species are",sum(Bfin[max(y),]>0),"out of",nsp))
