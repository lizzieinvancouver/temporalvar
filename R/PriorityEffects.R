#Priority Effects - main code
rm(list=ls()) 
#libraries and working directories, oh my!
library(deSolve)
library(scales)
require(MultiRNG)
if(length(grep("lizzie", getwd())>0)) {
    setwd("~/Documents/git/projects/temporalvar/R")
    }
library(here)
here()


#define the run - Consider creating a dataframe with combinations of parms to test
nruns<-1:1
outputy<-data.frame()
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
   #source(here("R","sourcefiles","PriEff_Output.R"))
  mydat<-data.frame(sp1_Rstar=Rstar[1],sp2_Rstar=Rstar[2],sp1_sense=xi_tau[1],sp2_sense=xi_tau[2],sp1.Bfin=Bfin[300,1],sp2.Bfin=Bfin[300,2])
  outputy<-rbind(mydat,outputy)
}

outputy$coexist<-ifelse(outputy$sp1.Bfin!=0 & outputy$sp2.Bfin!=0,"yes","no")
coexist<-dplyr::filter(coexist,sp2.Bfin!=0)

library(ggplot2)
ggplot(outputy,aes(sp1_Rstar))+geom_histogram(aes(fill=coexist),bins=100)
ggplot(outputy,aes(sp2_Rstar))+geom_histogram(aes(fill=coexist),bins=100)


outputy2<-dplyr::filter(outputy,sp1_Rstar<1)
outputy2<-dplyr::filter(outputy2,sp2_Rstar<1)
outputy2$Rratio<-outputy2$sp1_Rstar/outputy2$sp2_Rstar

outputy2$senratio<-outputy2$sp1_sense/outputy2$sp2_sense

ggplot(outputy2,aes(Rratio))+geom_histogram(aes(fill=coexist),bins=100)
ggplot(outputy2,aes(senratio))+geom_histogram(aes(fill=coexist),bins=100)

# #plot for testing
# par(mfrow=c(1,1))
# plot(seq(1,nyrs,dt),Bfin[,2],col="blue")
# points(seq(1,nyrs,dt),Bfin[,1],col="red")
# 


