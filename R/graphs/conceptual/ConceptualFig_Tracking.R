#Conceptual Figure to define key tracking concepts
#left panel show three environment timeseries:  temperature, daylength, incident PAR for two different locations (2 lines per panel)
#right panel shows three types of filters:  
#     fitness filter (fitness given time of event t): make one envt give clear peak, and second envt give wide peak
#     two cue filter panels (prob of event given envt) for two envts, showing ENVT'L TRACKING = how the same cue function results in differnet timing
#     measurement filter (what researchers measure)
# inset somewhere showing correlation between ideal day of event v predicted day of event = cue quality
# 
rm(list=ls())
library(tidyverse)
library(dplyr)
library(lubridate)
setwd("C:/Users/Megan/Documents/GitHub/temporalvar/R")

# get environmental data from NEON
# load packages
library(neonUtilities)
library(geoNEON)
library(raster)
library(rhdf5)
# Set global option to NOT convert all character variables to factors
options(stringsAsFactors=F)

#dataproducts:  IR biological temp, PAR, precip, soil water content, airtemp 
dpIDs <- c("DP1.00005.001","DP1.00024.001","DP1.00006.001","DP1.00094.001","DP1.00002.001") 
# IRBT <- loadByProduct(dpIDs[1],site=c("SRER", "JORN", "MOAB"), 
#                       startdate="2017-07", enddate="2020-02",
#                       nCores=2)
# precip <- loadByProduct(dpIDs[3],site=c("SRER", "JORN", "MOAB"), 
#                         startdate="2017-07", enddate="2020-02",
#                         nCores=2)
# par <- loadByProduct(dpIDs[2],site=c("SRER", "JORN", "MOAB"), 
#                      startdate="2017-07", enddate="2020-02",
#                      nCores=2)
# swc <- loadByProduct(dpIDs[4],site=c("SRER", "JORN", "MOAB"), 
#                      startdate="2017-07", enddate="2020-02", avg=30,
#                      nCores=2)
#  airt <- loadByProduct(dpIDs[5],site=c("SRER", "JORN", "MOAB"), 
#                       startdate="2017-07", enddate="2020-02", 
#                       nCores=2)
# save(irbt, precip,par,swc, airt, file="neondata/neon.Rdata")
load("neondata/neon.Rdata")

irbt %>% 
  group_by(siteID) %>%
  group_by()
  summarise((d.min =min(bioTempMinimum)))


#for fitness calc, create threshold values
#growth increments with temp
rT.aa <- 2
rT.min <- 15
rT.max <- 35
rT[d] <- rT.aa * temp[d] * (temp[d] - rT.min) * (rT.max - temp[d])^(1 / 2)
surv.T <- mintemp[d]< (-5)

rSM[d] <- 



ndays <- 180
fit <- rep(1,ndays)
for (d in seq(1,ndays,1)){
  for (dd in seq(1,d,1)){
  fit[dd]<- fit[dd] + rT[d]
  }
}
