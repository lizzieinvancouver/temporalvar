### Started 10 July 2011 ###
### By Lizzie & Megan ###

## VarEnvironments & Coexistence ##

## 8 Jan 2016 we deleted all calculations of E and C (and relations)
# the last commit with those calcs was 4af34aa, the 131st commit ##

# safety feature(s)
rm(list=ls()) 
options(stringsAsFactors=FALSE)

# set the working directory
# print("Hey! you need to set the correct working directory here")
setwd("/n/regal/wolkovich_lab/tempvar")

# packages
library(ggplot2)
library(deSolve)

set.seed(2)

## Setting up loop for multiple model runs
# To do still:
# (1) add in everything sppvars need (was lazy about this)
# (2) add in resource stuff to withinyrs
# (3) decide on list for each run, versus some other format
# (4) Figure out how to write output here and on RC 

# set up model runs
modelruns <- list() # place to store output of runs 
nruns <- 4 # number of model runs to do
nonsta = c(100,100,100)   #number of [1] initial stationary,[2]nonstationary,[3]final nonstationary years
tracking = 1   #tracking in these runs?
varRstar = 1   #flag for variation in Rstar; if 1, then c is drawn randomly, and R* varies
nsp = 2        #Number of species to start in these runs?
source("sourcefiles/getRunParms.R") #define runtime parameters

numcoexist <- matrix(NA,nruns,length(nonsta))

for (j in c(1:nruns)){ # assuming, we will vary species characteristics between yrs ... 
  
source("sourcefiles/getGraphParms.R")  #define graphics parameters
source("sourcefiles/getEnvt.R")  #get constant and time-varying envt parms
source("sourcefiles/getSpecies.R")  #get species characteristics and Rstar

#Define arrays
#interannual dynamics set-up (R0 is in getEnvt.R)
N0 <- rep(100,nsp)          # initial number of seeds (per meter square?)
N <- matrix(rep(0), nyrs, nsp) # number of seeds by yr and spp
N[1,] <- N0  #initialize
rcrt <-  matrix(rep(0),nyrs,nsp) # recruitment in year y

## Within-season dynamics set-up
Bfin <- matrix(rep(0),nyrs,nsp) # biomass at end of year y
B0  <- matrix(rep(0),nyrs,nsp) # biomass at beginning of year y
Bout <- list() #each year has a dataframe with time,R(t),Bi(t) of dims(2+nsp,tsteps)
source("sourcefiles/ResCompN.R") # define within-season ode solver

## And away we go!
for (y in c(1:(nyrs-1))){
  #include a flag for runs with initial and final stationary periods, but no nonstationary period
  if (y==(nonsta[1]+1) && nonsta[2]==0) N[y,] <- N0 #if no ns period, reset for final stationary period
  #get initial biomass for year y
  B0[y,] <- b*g[y,]*N[y,] 
  #use deSolve for ResCompN
  R<-R0[y]
  B<-B0[y,]
  State<-c(R=R,B=B)
  Time <- seq(0,ndays,by=dt)
  Bout[[y]] <- as.data.frame(ode(func = ResCompN, y = State, parms = Pars, times = Time,
     rootfun=rootfun))
  #Bout[[y]] <- as.data.frame(lsodar(func = ResCompN, y = State, parms = Pars, times = Time,rootfun=rootfun))
  Bfin[y,] <-  apply(Bout[[y]][3:(2+nsp)],2,FUN=max)  #final biomass
  N[y+1,] <- N[y,]*s*(1-g[y,]) + phi*Bfin[y,]    #note Bfin already includes N(t) as init cond; USES g here!
  N[y+1,] <- N[y+1,]*(N[y+1,]>ext)  #if density does not exceed ext, set to zero
  if (!(sum(N[y+1,]>0))) break    #if all species have gone extinct, go to next run
  rcrt[y,] <- g[y,]*(phi*Bfin[y,]-s)   #to recruit, convert biomass to seeds and overwinter
  #set Rstarmin threshold & update rootfun; this accounts for case where the spp with min R* goes extinct
  Rstarmin <- min(Rstar[(N[y+1,]!=0)])
  rootfun <- function(Time, State, Pars) State[1] -Rstarmin
}

## write out the model run output
# numcoexist records the number of species above 0 close to end of the period of each run)
numcoexist[j,1] <- sum(Bfin[(nonsta[1]-2),]>0)
if (nonsta[2]>0) numcoexist[j,2] <- sum(Bfin[(nonsta[1]+nonsta[2]-2),]>0)
if (nonsta[3]>0) numcoexist[j,3] <- sum(Bfin[(nyrs-2),]>0) 

## modelruns includes the variables that are constant across years in one dataframe...
# then tauI, tauP and Bfin for each year
modelruns[[j]] <- list("sppvars"=sppvars, "tauI"=tauI, "tauP"=tauP, "Bfin"=Bfin)

## could also make each run a multi-part dataframe with common names
# so something like:
# modelruns[[paste("crossyrs", j, sep="")]] <- sppvars
# modelruns[[paste("withinyrs", j, sep="")]] <- Bfin
## This is how I did it in the other file
# sppabund <- data.frame(matrix(ncol = nsp, nrow = nruns))
# colnames(sppabund) <- c(1:nsp)
# modelparams <- data.frame(nsp=numeric(), lambda=numeric(), gmax=numeric(),
#     s=numeric(), h=numeric(), p=numeric(), q=numeric())

} # close the model run loop

save(modelruns, file="modelruns.RData")
