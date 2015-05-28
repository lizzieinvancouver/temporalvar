### Started 28 Jan 2015 ###
### By Lizzie & Megan ###

## Simpler model for within year dynamics (L-V within year), with variable between year germination ##

# safety feature(s)
setwd(getwd()) # Lizzie: setwd("~/Documents/git/temporalvar/R")
options(stringsAsFactors=FALSE)

# packages
library(ggplot2)
library(deSolve)

set.seed(2)

#Number of species
nsp = 2  #when nsp=2, tauI is assigned known values from chesson 2004

# multiple runs
nruns <- 10

# Set up dataframes to hold data from multiple runs
# Assumes lambda, gmax, s, etc. (see below) are constant across species
sppabund <- data.frame(matrix(ncol = nsp, nrow = nruns))
colnames(sppabund) <- c(1:nsp)
modelparams <- data.frame(nsp=numeric(), lambda=numeric(), gmax=numeric(),
    s=numeric(), h=numeric(), p=numeric(), q=numeric())
coexist <- data.frame(meanFit=numeric(), storEff=numeric(), relN=numeric())
   # need to finish coexist below once finalized (format here may not be correct)
alphas <- array(dim=c(nsp, nsp, nruns))

for (arun in (c(1:nruns))){

source("sourcefiles/simple/getVaryingTauAlpha.R")  #get species characteristics
ext <- 1/100000  #Extinction Threshold:  1 seed/ha (assuming that initial density is 10 seeds per meter)

#Number of years to run
nyrs <- 1000

#Define arrays
#interannual dynamics set-up (R0 is in getEnvt.R)
N0 <- rep(10,nsp)          # initial number of seeds (per meter square?)
N <- matrix(rep(0), nyrs, nsp) # number of seeds by yr and spp
N[1,] <- N0  #initialize
#rcrt <- matrix(rep(0),nyrs,nsp) # recruitment in year y
#rcrt0 <- matrix(rep(0),nyrs,nsp) # recruitment WO competition in year y

#Define environment parameters & germination
p <- 2  #first parameter for beta distribution of tau
q <- 2  #second parameter for beta distribution of tau
tauP <- rbeta(nyrs, p, q) # change once not doing stationary+nonstationary run
h <-  rep(100,nsp)      # max rate of germination decrease following pulseif (nsp==2) tauI <- c(0.35, 0.4) else tauI <-runif(nsp,0.1, 0.9)  # time of max germ for sp i
g <- gmax*exp(-h*(matrix(rep(tauP,nsp),nrow=length(tauP),ncol=nsp)-
    matrix(rep(tauI,nyrs),ncol=nsp,nrow=nyrs,2))^2)  #germination fraction in year y


for (y in c(1:(nyrs-1))){
  C <- rep(NA,nsp)
  for (i in c(1:nsp)){
    a <- sum(alpha[i,]*g[y,]*N[y,])
  }
  N[y+1,] <- N[y,]*s + g[y,]*N[y,]*(lambda/(1+a) - s)
  N[y+1,] <- N[y+1,]*(N[y+1,]>ext)  #if density does not exceed ext, set to zero
}

source("sourcefiles/simple/Envt & Comp Calcs.R")
source("sourcefiles/simple/Fitness_Components.R")

sppabund[arun,] <- colMeans(N[1:5,], na.rm=TRUE)
modelparams[arun,] <- data.frame(nsp=nsp, lambda=lambda[1], gmax=gmax[1],
    s=s[1], h=h[1], p=p, q=q)
alphas[,,arun] <- alpha

}
     
#most basic plot
plot(c(1:nyrs),N[,1],type="l",ylim=c(0,max(N)), main="Abundance")
lines(c(1:nyrs),N[,2],col="red")
#plot(c(1:nyrs),g[,1],type="l",col="black",main = "Germination")
#lines(c(1:nyrs),g[,2],col="red")

#Fitness Components plot by species
bars<- matrix(c(meanFit + storEff - relNL,meanFit,storEff,-relNL),nrow=4,ncol = nsp, byrow = TRUE)
barplot(bars,main= "Fitness Components by species",beside = TRUE, xlab = "Species", col = c("black","gray","red","blue"),legend=c("p.c. Growth Rate","Mean Fitness","Storage Effect", "Rel Nonlin"),args.legend=c(cex=0.7))

