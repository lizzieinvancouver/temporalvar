### Started 10 July 2011 ###
### By Lizzie & Megan ###

## VarEnvironments & Coexistence ##

## 8 Jan 2016 we deleted all calculations of E and C (and relations)
# the last commit with those calcs was 4af34aa, the 131st commit ##

# safety feature(s)
rm(list=ls()) 
options(stringsAsFactors=FALSE)

# set the working directory
print("Hey! you need to set the correct working directory here")
setwd(getwd()) # Lizzie: setwd("~/Documents/git/projects/temporalvar/R")

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
nruns <- 40 # number of model runs to do
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

## Helpful reminders on indexing within lists
# how to call something (given that you have at least 2 runs) ....
modelruns[[2]][1][[1]]$gmax
# how to call internal list by name
modelruns[[2]]["sppvars"]
# when slicing out bits of vector within a list you have to get beyond the $ 
# this is often something I forget so here are two examples
modelruns[[2]]["sppvars"][[1]]$gmax
modelruns[[2]]["tauP"][[1]][1:10]

#source("sourcefiles/plotNyears.R")  #plots dynamics of seedbank abundance over years
#source("sourcefiles/plotBinSeason.R")  #plot within season dynamics of biomass & R for a subset of years
#source("sourcefiles/plotBinSeason_Lizzie.R")

## Did program ever jump out of loop? YES!
for (i in c(1:100)){
   print(length(Bout[[i]]$time))
}
# looking at R for last run (compare to min(Rstar))
Bout[[99]]$R

### Megan stopped tweaking plots here, but they will need to be adjusted for new within-year output structure from ode

# within years plots
# tweaked a little to make biomass differences among years clearer
# dev.new(width=14, height=10)
# par(mfrow=c(3,3))
# par(mar=c(5, 4, 4, 8) + 0.1)
# selectyrs <- seq(1, nyrs, by=floor(nyrs/8))
# yrlength <- ndays/dt
# for (yr in seq_along(selectyrs)){
#     allsp <- B[selectyrs[yr], ,]
#     plot(B[selectyrs[yr], 1,]~c(1:yrlength), ylim=c(min(allsp), max(B)), 
#         xlab="step, step, step",  ylab="Abundance", type="n",
#         main=paste("year: ", selectyrs[yr], sep=""))
#     for (sp in c(1:nsp)){
#         lines(B[selectyrs[yr], sp ,]~c(1:yrlength), ylim=c(min(allsp), max(B)),
#             col=colerz[sp])
#     }
#     par(new=TRUE)
#     plot(R[selectyrs[yr],]~c(1:yrlength), axes=FALSE, xlab="", ylab="",
#         ylim=c(min(R[selectyrs[yr],]), max(R[selectyrs[yr],])), type="l",
#         col=rcol, lty=lresbyrs, lwd=lwd)
#     raxis <- seq(0, max(R[selectyrs[yr],]), by=max(R[selectyrs[yr],])/10)
#     axis(4, at=raxis, labels=round(raxis, digits=2))
#     mtext("Resource", side=4, line=3, cex=0.75)
# }


dev.new(width=4, height=3)

# (2) using ggplot, which really is good for this sort of thing
# need to update this now that tauI varies by year
tau.df <- data.frame(coexisted=Bfin[max(y),]>0, tauI=colMeans(tauI[,]),
    tauIhat=colMeans(tauIhat[,]), alpha=colMeans(alpha[,]))
tau.df$coexisted[tau.df$coexisted==TRUE] <- "coexisted"
tau.df$coexisted[tau.df$coexisted==FALSE] <- "doomed"
tauP.df <- data.frame(coexisted=rep("tauP"), tauP=tauP)


if (nonsta[1]>0){
  ggplot(tau.df, aes(tauI, fill = coexisted)) + geom_histogram(alpha=0.5) +
    geom_density(data=tauP.df[1:(nyrs-sum(nonsta)),], aes(tauP),  alpha = 0.2) +
    geom_density(data=tauP.df[(sum(nonsta)+1):nyrs,], aes(tauP),  alpha = 0.4) +
    labs(title=paste(sum(Bfin[max(y),]>0), "out of", nsp, "coexisted", sep=" "))
}else {
  ggplot(tau.df, aes(tauI, fill = coexisted)) + geom_histogram(alpha=0.5) +
    geom_density(data=tauP.df[1:nyrs,], aes(tauP),  alpha = 0.2) +
    labs(title=paste(sum(Bfin[max(y),]>0), "out of", nsp, "coexisted", sep=" "))
}


# plot(tauP, type="l")
# lines(tauI[,34], col="red", pch=16, cex=0.5) 
# lines(tauI[,49], col="blue", pch=16, cex=0.5)

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
