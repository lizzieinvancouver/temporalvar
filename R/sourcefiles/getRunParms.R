# safety feature(s)
rm(list=ls()) 
options(stringsAsFactors=FALSE)

# packages
library(deSolve)

# set the working directory
#setwd(getwd()) # Lizzie: setwd("~/Documents/git/projects/temporalvar")
#setwd(getwd()) # Megan: setwd("~/Documents/GitHub/temporalvar")
#setwd("~/n/regal/wolkovich_lab/temporalvar")

#random seed
#set.seed(2)

#Runtime Parameters
runname <- "Track_varR_2spp"
batch <-1  #flag if running as batch array w SLURM
nruns <- 100
writeBout <- 1  #flag indicating how often Bout should be written (0=never, n = every n runs)
nonsta = c(200,0,0)   #number of [1] initial stationary,[2]nonstationary,[3]final nonstationary years
tracking = 1   #tracking in these runs?
varRstar = 1   #flag for variation in Rstar; if 1, then c is drawn randomly, and R* varies
nsp = 2        #Number of species to start in these runs?
jobID <- ""
if (batch==1) {jobID <- Sys.getenv(c("SLURM_ARRAY_JOB_ID","SLURM_ARRAY_TASK_ID"))}


#between year
nyrs <- sum(nonsta)  # number of yrs to run if nonsta=0 or for initial period if nonsta>0

y <- c(1:nyrs)
ext <- 1/10000  #Extinction Threshold:  1 seed/ha (assuming that initial density is 10 seeds per meter)

#within year
ndays <- 5  # number of days in a growing season
dt <- 0.005 # within yr timestep
tsteps <- ndays/dt
