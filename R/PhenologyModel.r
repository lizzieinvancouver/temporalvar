### Started 10 July 2011 ###
### By Lizzie & Megan ###

## VarEnvironments & Coexistence ##

# safety feature(s)
options(stringsAsFactors=FALSE)

# packages
library(ggplot2)
library(deSolve)

set.seed(2)

## Setting up loop for multiple model runs
# To do still:
# (1) add in everything crossyrsvars need (was lazy about this)
# (2) add in resource stuff to withinyrs
# (3) decide on list for each run, versus some other format

#Temporarily disabling the multiple runs
#modelruns <- list() # place to store output of runs
#nruns <- 2 # number of model runs to do
#for (j in c(1:nruns)){ # assuming, we will vary species characteristics between yrs ... 
  
#Stationarity in this run?
nonsta = 0  #flag for stationary (0) vs nonstationary (=num yrs nonstationary)

#Number of species to start?
nsp = 5  #when nsp=2, tauI is assigned known values from chesson 2004

source("getRunParms.R") #define runtime parameters

source("getGraphParms.R")  #define graphics parameters

source("getEnvt.R")  #get constant and time-varying envt parms

source("getSpecies.R")  #get species characteristics and Rstar

#Define arrays
#interannual dynamics set-up (R0 is in getEnvt.R)
N0 <- rep(100,nsp)          # initial number of seeds (per meter square?)
N <- matrix(rep(0), nyrs, nsp) # number of seeds by yr and spp
N[1,] <- N0  #initialize
rcrt <- matrix(rep(0),nyrs,nsp) # recruitment in year y
rcrt0 <- matrix(rep(0),nyrs,nsp) # recruitment WO competition in year y

## Within-season dynamics set-up
Bfin <- matrix(rep(0),nyrs,nsp) # biomass at end of year y
B0  <- matrix(rep(0),nyrs,nsp) # biomass at beginning of year y
Bout <- list() #each year has a dataframe with time,R(t),Bi(t) of dims(2+nsp,tsteps)
source("ResCompN.R") # define within-season ode solver
source("NoCompN.R")  # define within-season ode solver for no competition

## set-up for different coexistence mechanisms
#I have considered 3 defns for E and C, but only one isn't problematic numerically
  #defn1: E= ln(rcrt0), C=ln(rcrt0/rcrt)  ->compare recruitment w and wo comp
  #defn2: E= ln(g), C=-ln(phi*Bfin-s) -> this the easiest division of the eqn for rcrt
  #defn3: E =ln(g), C=-ln(phi*Bfin)  -> in this case, (1-g) is included in the seedback lifespan
  #defn4: E = ln(g*phi*BnoC), C=ln(BnoC/Bfin)
BnoCout <- list()
BnoC <- matrix(rep(0),nyrs,nsp) # B without competition at end of year y
E <- matrix(rep(0),nyrs,nsp)  
C <- matrix(rep(0),nyrs,nsp)  

for (y in c(1:(nyrs-1))){
  #get initial biomass for year y
  B0[y,] <- b*g[y,]*N[y,]

  #use deSolve for ResCompN
  R<-R0[y]
  B<-B0[y,]
  State<-c(R=R,B=B)
  Time <- seq(0,ndays,by=dt)
  #Bout[[y]] <- as.data.frame(ode(func = ResCompN, y = State, parms = Pars, times = Time))
  Bout[[y]] <- as.data.frame(lsodar(func = ResCompN, y = State, parms = Pars, times = Time,rootfun=rootfun))
  Bfin[y,] <-  apply(Bout[[y]][3:(2+nsp)],2,FUN=max)  #final biomass

  #use deSolve for NoCompN to solve for noCompetition condition
  #would be faster to used solved equation, but calculations were coming out wrong
  #tstar is when species cross their Rstar threshold and we stop the season under the no competition condition; checked against ODE solver, it is when the  biomass starts decreasing
  tstar <- (1/eps)*(log(R0[y]) - (1/theta)*log(m/(a*c-a*u*m)))  
  TimeNC <- seq(0,max(tstar)+5*dt,by=dt)
  BnoCout[[y]] <- as.data.frame(ode(func = NoCompN, y = State, parms = Pars, times = TimeNC))
  BnoC[y,] <- apply(BnoCout[[y]][3:(2+nsp)],2,FUN=max)
  
  rcrt[y,] <- g[y,]*(phi*Bfin[y,]-s)   #to recruit, convert biomass to seeds and overwinter
  #rcrt0[y,] <- g[y,]*(phi*BnoC[y,]-s)   #no-competition recruitment
  
  #calculate E and C
  E[y,] <- log(g[y,]*phi*BnoC[y,])         #defn 4
  C[y,] <- log(BnoC[y,]/Bfin[y,])         #defn 4   
  
  N[y+1,] <- N[y,]*(s + rcrt[y,])    #N(t+1) = N(t)* (survival + recruitment)
  N[y+1,] <- N[y+1,]*(N[y+1,]>ext)  #if density does not exceed ext, set to zero
}
                                                                     

#modelruns[[j]] <- list(crossyrsvars, Bfin, E)

# could also make each run a multi-part dataframe with common names
# so something like:
# modelruns[[paste("crossyrs", j, sep="")]] <- crossyrsvars
# modelruns[[paste("withinyrs", j, sep="")]] <- Bfin

#}

source("plotNyears.R")  #plots dynamics of seedbank abundance over years

source("plotBinSeason.R")  #plot within season dynamics of biomass & R for a subset of years

source("plotBwCnoC.R")  #plot within season biomass resource dynamics w and wo competition


###Megan stopped tweaking plots here, but they will need to be adjusted for new within-year output structure from ode

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
