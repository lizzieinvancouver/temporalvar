# safety feature(s)
rm(list=ls()) 
options(stringsAsFactors=FALSE)

# packages
library(deSolve)

# set the working directory
setwd(getwd()) # Lizzie: setwd("~/Documents/git/projects/temporalvar/R")
#SCRATCH:  setwd("~/n/regal/wolkovich_lab/tempvar")

#random seed
set.seed(2)

#Runtime Parameters
batch <-0  #flag if running as batch file w SLURM
nruns <- 5
writeBout <- 0  #flag indicating whether Bout should be written
nonsta = c(100,100,100)   #number of [1] initial stationary,[2]nonstationary,[3]final nonstationary years
tracking = 1   #tracking in these runs?
varRstar = 1   #flag for variation in Rstar; if 1, then c is drawn randomly, and R* varies
nsp = 2        #Number of species to start in these runs?
jobID <- ""
if (batch==1) {jobID <- Sys.getenv("SLURM_JOB_ID")}


#between year
nyrs <- sum(nonsta)  # number of yrs to run if nonsta=0 or for initial period if nonsta>0

y <- c(1:nyrs)
ext <- 1/10000  #Extinction Threshold:  1 seed/ha (assuming that initial density is 10 seeds per meter)

#within year
ndays <- 5  # number of days in a growing season
dt <- 0.005 # within yr timestep
tsteps <- ndays/dt
