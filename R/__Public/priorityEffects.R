#Priority Effects - main code for review and public
rm(list=ls()) 
#libraries and working directories, oh my!
library(deSolve)
library(scales)
library(dplyr)
require(MultiRNG)
library(ggpubr)
if(length(grep("lizzie", getwd())>0)) {
  setwd("~/Documents/git/projects/temporalvar/R")
}

if(length(grep("dbuona", getwd())>0)) {
  setwd("~/Documents/git/temporalvar/")
}
library(here)
here()

#define the run 
nruns<-10 ## number of iterations
g_notxi <- 1 ##this determines whether or not germination fraction varies with climate. 1=no.
c_warm<- 0.5 ## frequencey of drawing from current and future cliamte distrubtions
makeplots <- TRUE # whether or not plots are made
outputy<-data.frame() #this is the main results file that can be written out for plotting.

for (j in c(1:nruns)) {
  
  nyrs <- 300 # how many years is each run
  #define the environment for this run
  source(here("R","sourcefiles","PriEff_Envt.R"))
  
  #define the species in this run
  source(here("R","sourcefiles","PriEff_Species.R"))
  
  #run the model for nyrs
  source(here("R","sourcefiles","PriEff_Comp.R"))
  source(here("R","sourcefiles","PriEff_Model.R"))
  
  # #write out the results for this run
  if (g_notxi==1) {
    source(here("R","sourcefiles","PriEff_Output.R"))
    
  }else {
    source(here("R","sourcefiles","PriEff_OutputwFract.R"))
  }
}

