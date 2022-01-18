#Priority Effects - main code

#libraries
library(deSolve)
require(MultiRNG)

#define the run:  read in a csv file with the following headers:
# runname, nyrs,
run_parms <- read.csv("PriEff_runparms.csv",header = TRUE)

#define the environment for this run
source("/sourcefiles/PriEff_Envt.R")

#define the species in this run
source("/sourcefiles/PriEff_Species.R")

#run the model for nyrs
source("/sourcefiles/PriEff_Comp.R")
source("/sourcefiles/PriEff_Model.R")



