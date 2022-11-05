#Priority Effects - main code
rm(list=ls()) 
#libraries
library(deSolve)
require(MultiRNG)
library(here)
here()

#define the run - Consider creating a dataframe with combinations of parms to test

for (n in c(1,nruns)) {

  nyrs <- 300
  
  #define the environment for this run
  source(here("R","sourcefiles","PriEff_Envt.R"))
  
  #define the species in this run
  source(here("R","sourcefiles","PriEff_Species.R"))
  
  #run the model for nyrs
  source(here("R","sourcefiles","PriEff_Comp.R"))
  source(here("R","sourcefiles","PriEff_Model.R"))

  # #write out the results for this run
  # source(here("R","sourcefiles","PriEff_Output.R"))
  
}

# #plot for testing
# par(mfrow=c(1,1))
# plot(seq(1,nyrs,dt),Bfin[,2],col="blue")
# points(seq(1,nyrs,dt),Bfin[,1],col="red")
# 


