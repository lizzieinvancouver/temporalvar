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
library(zoo)
setwd("C:/Users/Megan/Documents/GitHub/temporalvar/R")

# get environmental data from NEON
# load packages
library(neonUtilities)
#library(geoNEON)
library(raster)
library(rhdf5)
# Set global option to NOT convert all character variables to factors
options(stringsAsFactors=F)


# get Neon Data -----------------------------------------------------------

#dataproducts:  IR biological temp, PAR, precip, soil water content, airtemp, soil temp 
sites <- c("ABBY","WREF","JORN","MOAB")

if(file.exists("neondata/neon.Rdata")) {
  load("neondata/neon.Rdata")  
} else {
  
  startD <- "2017-04"
  endD <- "2020-04"
  dpIDs <- c("DP1.00005.001",  #biotemp
             "DP1.00024.001",  #par
             "DP1.00006.001",  #precip
             "DP1.00094.001",  #soil water content
             "DP1.00002.001",  #airtemp
             "DP1.00041.001")  #soiltemp
  
  # irbt <- loadByProduct(dpIDs[1],site=sites,
  #                       startdate=startD, enddate=endD,
  #                       nCores=2)
  # par <- loadByProduct(dpIDs[2],site=sites,
  #                      startdate=startD, enddate=endD,
  #                      nCores=2)
  # precip <- loadByProduct(dpIDs[3],site=sites,
  #                         startdate=startD, enddate=endD,
  #                         nCores=2)
  # swc <- loadByProduct(dpIDs[4],site=sites,
  #                      startdate=startD, enddate=endD, avg=30,
  #                      nCores=2)
  # airt <- loadByProduct(dpIDs[5],site=sites,
  #                       startdate=startD, enddate=endD,
  #                       nCores=2)
  # soilt <- loadByProduct(dpIDs[6],site=sites,
  #                       startdate=startD, enddate=endD,
  #                       nCores=2)
  
  stackByTable("neondata/NEON_temp-bio.zip",nCores = 2)
  #stackByTable("neondata/NEON_temp-soil.zip",nCores = 2)
  stackByTable("neondata/NEON_temp-air-single.zip",nCores = 2)
  #stackByTable("neondata/NEON_conc-h2o-soil-salinity.zip",nCores = 2)
  stackByTable("neondata/NEON_precipitation.zip",nCores = 2)
  #soil files don't work with this function, so do by hand later
  
  #precip1 <- read_csv("neondata/NEON_precipitation/stackedFiles/PRIPRE_30min.csv")
  precip2 <- read_csv("neondata/NEON_precipitation/stackedFiles/SECPRE_30min.csv")
  bioT <- read_csv("neondata/NEON_temp-bio/stackedFiles/IRBT_30_minute.csv")
  airT <- read_csv("neondata/NEON_temp-air-single/stackedFiles/SAAT_30min.csv")
  
  save(bioT, precip1, precip2, airT, file="neondata/neon.Rdata")
}

# Make dataframe ----------------------------------------------------------
airT.thresh.lo <- 0
airT.thresh.hi <- 32
bioT.thresh.lo <- 0
bioT.thresh.hi <- 32

d <- 
  left_join(
    airT %>% 
      filter(verticalPosition=="010") %>% 
      mutate(lt.lo = if_else(tempSingleMean<airT.thresh.lo,1,0),
             gt.hi = if_else(tempSingleMean>airT.thresh.hi,1,0)) %>% 
      group_by(site = siteID, 
               date = date(startDateTime)) %>%
      summarise(airT.mean = mean(tempSingleMean, na.rm=TRUE),
                airT.min = ifelse(sum(!is.na(tempSingleMean))>0,min(tempSingleMean,na.rm=TRUE),NA),
                airT.max = ifelse(sum(!is.na(tempSingleMean))>0,max(tempSingleMean,na.rm=TRUE),NA),
                airT.lt.lo = ifelse(sum(!is.na(lt.lo))>0, sum(lt.lo, na.rm=TRUE)*0.5,NA),
                airT.gt.hi = ifelse(sum(!is.na(gt.hi))>0, sum(gt.hi, na.rm=TRUE)*0.5,NA),
                airT.dch = ifelse(sum(!is.na(lt.lo))>0, sum(lt.lo*0.5 * (airT.thresh.lo - tempSingleMean), na.rm=TRUE),NA),
                airT.dhh = ifelse(sum(!is.na(gt.hi))>0, sum(gt.hi*0.5 * (tempSingleMean - airT.thresh.hi), na.rm=TRUE),NA),
                airT.ct = sum(!is.na(tempSingleMean))) %>%  
      mutate(airT.dhh.5d = rollapplyr(airT.dhh, 5, sum, na.rm=TRUE, fill=NA)) %>% 
      mutate(airT.dch.90d = rollapplyr(airT.dch, 90, sum, na.rm=TRUE, fill=NA)),
    bioT %>%  
      filter(verticalPosition == "010") %>% 
      mutate(lt.lo = if_else(bioTempMean<bioT.thresh.lo, 1, 0),
             gt.hi = if_else(bioTempMean>bioT.thresh.hi, 1, 0)) %>% 
      group_by(site = siteID,
               date = date(startDateTime)) %>% 
      summarise(bioT.mean = mean(bioTempMean, na.rm=TRUE),
              bioT.min = ifelse(sum(!is.na(bioTempMean))>0,min(bioTempMean,na.rm=TRUE),NA),
              bioT.max = ifelse(sum(!is.na(bioTempMean))>0,max(bioTempMean,na.rm=TRUE),NA),
              bioT.lt.lo = ifelse(sum(!is.na(lt.lo))>0, sum(lt.lo, na.rm=TRUE)*0.5,NA),
              bioT.gt.hi = ifelse(sum(!is.na(gt.hi))>0, sum(gt.hi, na.rm=TRUE)*0.5,NA),
              bioT.dhh = ifelse(sum(!is.na(lt.lo))>0, sum(gt.hi*0.5*(bioTempMean-bioT.thresh.hi), na.rm=TRUE),NA),
              bioT.dch = ifelse(sum(!is.na(gt.hi))>0, sum(lt.lo*0.5 *(bioT.thresh.lo - bioTempMean), na.rm=TRUE),NA),
              bioT.ct = sum(!is.na(bioTempMean))) %>%
      mutate(bioT.dhh.5d = rollapplyr(bioT.dhh, 5, sum, na.rm=TRUE, fill=NA)) %>% 
      mutate(bioT.dch.90d = rollapplyr(bioT.dch, 90, sum, na.rm=TRUE, fill=NA)),
    by=c("date","site")) %>% 
  left_join(.,
    precip2 %>%
      group_by(site = siteID, 
               date = date(startDateTime)) %>% 
      summarise(precip.daily = ifelse(sum(!is.na(secPrecipBulk))>0,sum(secPrecipBulk, na.rm = TRUE),NA),
                precip.daily.any = any(secPrecipBulk>0)) %>% 
      mutate(precip.cum.15d = rollapplyr(precip.daily, 15, sum, na.rm=TRUE, fill=NA)) %>% 
      mutate(precip.cum.30d = rollapplyr(precip.daily, 30, sum, na.rm=TRUE, fill=NA)),
    by=c("date","site")) #%>% 
  # left_join(.,
  #           swc %>% 
  #             group_by(site = siteID, 
  #                      date = date(startDateTime)) %>% 
  #             summarise(swc.mean = mean(VSWCMean), 
  #                       swc.min = mean(VSWCMinimum),
  #                       swc.max = mean(VSWCMaximum)),
  #           by=c("date","site")) %>% 
  # left_join(.,
  #         soilt %>% 
  #           group_by(site = siteID, 
  #                    date = date(startDateTime)) %>% 
  #           summarise(soilT.mean = mean(soilTempMean), 
  #                     soilT.min = mean(soilTempMinimum),
  #                     soilT.max = mean(soilTempMaximum)),
  #         by=c("date","site"))  

rm(airT,bioT,precip1,precip2)

# basic plots -------------------------------------------------------------

#use grid.arrange
d %>% 
  group_by(site) %>% 
  ggplot() +
  facet_wrap(vars(site)) +
  geom_line(aes(x=date,y=precip.daily,color=site)) +
  labs(x = "Day",y="Daily Precipitation")
# d %>% 
#   ggplot() +
#   geom_line(aes(x=date,y=swc.mean, color=site)) +
#   labs(x = "Day",y="Soil Moisture")
d %>% 
  group_by(site) %>% 
  ggplot() +
  facet_wrap(vars(site)) +
  geom_line(aes(x=date,y=airT.mean)) +
  geom_line(aes(x=date,y=airT.min),color="light blue") +
  geom_line(aes(x=date,y=airT.max),color="pink") +
  labs(x = "Day",y="Air Temperature")
d %>% 
  group_by(site) %>% 
  ggplot() +
  facet_wrap(vars(site)) +
  geom_line(aes(x=date,y=bioT.mean)) +
  geom_line(aes(x=date,y=bioT.min),color="light blue") +
  geom_line(aes(x=date,y=bioT.max),color="pink") +
  labs(x = "Day",y="Bio Temperature")
d %>% 
  group_by(site) %>% 
  ggplot() +
  facet_wrap(vars(site)) +
  geom_point(aes(x=bioT.mean,y=airT.mean))


# Fitness -----------------------------------------------------------------

#for fitness calc:
# growth is a unimodal function of temperature
# P(success) is 0 if subsequent freezing exceeds 4 hours, and declines with degree heating hours in excess of 40deg 
#growth increments depend on daily temperature
grT.aa <- .2
grT.min <- 15
grT.max <- 35
t <- seq(0,40,.2)
growth.T <- ifelse(between(t,15,35), grT.aa * t * (t - grT.min) * (grT.max - t)^(1 / 2),0)
plot(t,growth.T)
#survival depends on cumulative cold or heat stress
cumFr <- 4  #excess of 4 hours below freezing after EVENT
cumHe <- 48  #48 degree heating hours that exceed 48
surv.T <- ifelse(airT.lt0>4,0,ifelse(airT.gt40>)
surv.T <- rT.aa * t * abs((rT.min-t) * (t-rT.max))^(1 / 2)
plot(t,surv.T)
rSM[d] <- 



ndays <- 180
fit <- rep(1,ndays)
for (d in seq(1,ndays,1)){
  for (dd in seq(1,d,1)){
  fit[dd]<- fit[dd] + rT[d]
  }
}
