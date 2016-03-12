### Started 10 July 2011 ###
### By Lizzie & Megan ###
### This executes a single run of the model and writes out

#define run conditions
source("R/sourcefiles/getRunParms.R") #define runtime parameters

#define parameters and functions for this run
source("R/sourcefiles/getEnvt.R")  #get constant and time-varying envt parms
source("R/sourcefiles/getSpecies.R")  #get species characteristics and Rstar
source("R/sourcefiles/ResCompN.R") # define within-season ode solver

for (j in c(1:nruns)){
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
  
  ## modelruns includes the variables that are constant across years in one dataframe...
  # then tauI, tauP and Bfin for each year
  save(sppvars, tauI, tauP,Bfin,file=paste("R/output/",filename,"_", jobID[1],"-",jobid[2],"-run",j,".Rdata",sep="")) 
  if (writeBout>0) {
    save(Bout,file=paste("R/output/",filename,"_Bout_",jobID[1],"-",jobid[2],j,".Rdata",sep="")) #("out_",i,".Rdata"))
  }
}
