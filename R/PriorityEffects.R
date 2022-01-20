#Priority Effects - main code

#libraries
library(deSolve)
require(MultiRNG)
library(here)
here()

#define the run - Consider creating a dataframe with combinations of parms to test
nyrs <- 1

#define the environment for this run
source(here("R","sourcefiles","PriEff_Envt.R"))

#define the species in this run
source(here("R","sourcefiles","PriEff_Species.R"))

#run the model for nyrs
source(here("R","sourcefiles","PriEff_Comp.R"))
source(here("R","sourcefiles","PriEff_Model.R"))



