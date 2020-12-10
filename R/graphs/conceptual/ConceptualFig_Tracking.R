#Conceptual Figure to define key tracking concepts
#left panel show three environment timeseries:  temperature, daylength, incident PAR for two different locations (2 lines per panel)
#right panel shows three types of filters:  
#     fitness filter (fitness given time of event t): make one envt give clear peak, and second envt give wide peak
#     two cue filter panels (prob of event given envt) for two envts, showing ENVT'L TRACKING = how the same cue function results in differnet timing
#     measurement filter (what researchers measure)
# inset somewhere showing correlation between ideal day of event v predicted day of event = cue quality
# 

# Libraries Etc -----------------------------------------------------------

rm(list=ls())
library(tidyverse)
library(dplyr)
library(lubridate)
library(zoo)
library(neonUtilities)
library(raster)
library(rhdf5)
library(gridExtra)
library(RColorBrewer)
library(insol)
library(svglite)
options(stringsAsFactors=F)
setwd("C:/Users/Megan/Documents/GitHub/temporalvar/R")

# get Neon Data -----------------------------------------------------------

#dataproducts:  IR biological temp, PAR, precip, soil water content, airtemp, soil temp 
sites <- c("ABBY","WREF","JORN","MOAB")
sitelist <- data.frame(name = sites,
                       latitude = c(45.76244,45.82049,32.59069,38.24828),
                       longitude = c(-122.33032,-121.95191,-106.84254,-109.38827))

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
  # par <- loadByProduct(dpIDs[2],site="ABBY",
  #                      startdate=startD, enddate=endD)
  # precip <- loadByProduct(dpIDs[3],site=sites,
  #                         startdate=startD, enddate=endD,
  #                         nCores=2)
  # swc <- loadByProduct(dpIDs[4],site=sites,
  #                        startdate=startD, enddate=endD, avg=30)
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
  
  save(bioT, precip1, precip2, airT, par, swc, file="neondata/neon.Rdata")
}

# Make dataframe ----------------------------------------------------------
airT.thresh.lo <- 0
airT.thresh.hi <- 32
bioT.thresh.lo <- 0
bioT.thresh.hi <- 32
swc.thresh.lo <- 0.5

d <- 
  left_join(
    airT %>% 
      filter(verticalPosition=="010") %>% 
      mutate(lt.lo = if_else(tempSingleMean<airT.thresh.lo,1,0),
             gt.lo = if_else(tempSingleMean>airT.thresh.lo,1,0),
             gt.hi = if_else(tempSingleMean>airT.thresh.hi,1,0)) %>% 
      group_by(site = siteID, 
               date = date(startDateTime)) %>%
      summarise(airT.mean = mean(tempSingleMean, na.rm=TRUE),
                airT.min = ifelse(sum(!is.na(tempSingleMean))>0,min(tempSingleMean,na.rm=TRUE),NA),
                airT.max = ifelse(sum(!is.na(tempSingleMean))>0,max(tempSingleMean,na.rm=TRUE),NA),
                airT.lt.lo = ifelse(sum(!is.na(lt.lo))>0, sum(lt.lo, na.rm=TRUE)*0.5,NA),
                airT.gt.lo = ifelse(sum(!is.na(gt.lo))>0, sum(gt.lo, na.rm=TRUE)*0.5,NA),
                airT.gt.hi = ifelse(sum(!is.na(gt.hi))>0, sum(gt.hi, na.rm=TRUE)*0.5,NA),
                airT.dch = ifelse(sum(!is.na(lt.lo))>0, sum(lt.lo*0.5 * (airT.thresh.lo - tempSingleMean), na.rm=TRUE),NA),
                airT.dhh.hi = ifelse(sum(!is.na(gt.hi))>0, sum(gt.hi*0.5 * (tempSingleMean - airT.thresh.hi), na.rm=TRUE),NA),
                airT.dhh.lo = ifelse(sum(!is.na(gt.lo))>0, sum(gt.lo*0.5 * (tempSingleMean - airT.thresh.lo), na.rm=TRUE),NA),
                airT.ct = sum(!is.na(tempSingleMean))) %>%  
      mutate(airT.dhh.5dH = rollapplyr(airT.dhh.hi, 5, sum, na.rm=TRUE, fill=NA),
             airT.dch.90d = rollapplyr(airT.dch, 90, sum, na.rm=TRUE, partial=TRUE),
             airT.dhh.30dL = rollapplyr(airT.dhh.lo,30,sum,na.rm=TRUE,partial=TRUE)),
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
      mutate(precip.cum.15d = rollapplyr(precip.daily, 15, sum, na.rm=TRUE, fill=NA),
             precip.cum.30d = rollapplyr(precip.daily, 30, sum, na.rm=TRUE, fill=NA)),
    by=c("date","site")) %>% 
  left_join(.,
    par %>%
      group_by(siteID, startDateTime) %>% 
      summarise(PARMeanV = mean(PARMean,na.rm=TRUE)) %>% 
      ungroup %>% 
      group_by(site = siteID,
               date = date(startDateTime)) %>%
      summarise(par.cum = ifelse(sum(!is.na(PARMeanV))>0, sum(PARMeanV,na.rm=TRUE),NA)) %>% 
      mutate(par.7d = rollapplyr(par.cum, 7,mean,na.rm=TRUE, partial=TRUE)), 
    by=c("date","site")) %>% 
  left_join(.,
    swc %>%
      filter(verticalPosition == "501") %>%
      group_by(siteID, startDateTime) %>% 
      summarise(VSWCMeanH = mean(VSWCMean,na.rm=TRUE),
                VSWCMinH = min(VSWCMean,na.rm=TRUE),
                VSWCMaxH= max(VSWCMean, na.rm=TRUE)) %>% 
      ungroup %>% 
      mutate(lt.lo = if_else(VSWCMeanH < swc.thresh.lo,1,0)) %>%
      group_by(site = siteID,
               date = date(startDateTime)) %>%
      summarise(swc.mean = mean(VSWCMeanH, na.rm=TRUE),
                swc.min = ifelse(sum(!is.na(VSWCMinH))>0,min(VSWCMinH,na.rm=TRUE),NA),
                swc.max = ifelse(sum(!is.na(VSWCMaxH))>0,max(VSWCMaxH,na.rm=TRUE),NA),
                swc.lt.lo = ifelse(sum(!is.na(lt.lo))>0, sum(lt.lo, na.rm=TRUE)*0.5,NA),
                swc.ct = sum(!is.na(VSWCMeanH))),
    by=c("date","site"))

d <- left_join(d,sitelist,by=c("site"="name"))
d <- d %>% mutate(airT.dch.90d.nAll = airT.dch.90d/max(airT.dch.90d,na.rm="TRUE"),
                  airT.dhh.30dL.nAll = airT.dhh.30dL/max(airT.dhh.30dL,na.rm="TRUE"))
d <- mutate(d,jdate = yday(date))
d <- d %>% rowwise() %>% mutate(dayL = daylength(latitude,longitude,jdate,1)[3]) #use rowwise bc functio is not vectorized

rm(airT,bioT,precip1,precip2,swc,par)

# basic plots -------------------------------------------------------------

xax <- c(mdy("11/1/2017"),
         mdy("12/1/2017"),
         mdy("1/1/2018"),
         mdy("2/1/2018"),
         mdy("3/1/2018"),
         mdy("4/1/2018"),
         mdy("5/1/2018"),
         mdy("6/1/2018"),
         mdy("7/1/2018"),
         mdy("8/1/2018"),
         mdy("9/1/2018"),
         mdy("10/1/2018"))

# xax <- c(mdy("11/1/2017"),
#          mdy("1/1/2018"),
#          mdy("3/1/2018"),
#          mdy("5/1/2018"),
#          mdy("7/1/2018"),
#          mdy("9/1/2018"))
xaxL <-c("N","D","J","F","M","A","M","J","J","A","S","O")


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
  geom_line(aes(x=date,y=par.cum,color=site)) +
  labs(x = "Day",y="Cum. Daily PAR")
d %>% 
  group_by(site) %>% 
  ggplot() +
  facet_wrap(vars(site)) +
  geom_line(aes(x=date,y=swc.mean)) +
  geom_line(aes(x=date,y=swc.min),color="light blue") +
  geom_line(aes(x=date,y=swc.max),color="pink") +
  labs(x = "Day",y="Daily Soil Water Content")
d %>% 
  group_by(site) %>% 
  ggplot() +
  facet_wrap(vars(site)) +
  geom_point(aes(x=bioT.mean,y=airT.mean))
d %>% 
  group_by(site) %>% 
  ggplot() +
  facet_wrap(vars(site)) +
  geom_point(aes(x=swc.mean,y=precip.cum.30d))
d %>% 
  group_by(site) %>% 
  ggplot() +
  facet_wrap(vars(site)) +
  geom_point(aes(x=date,y=dayL))

# Add Growth, Survival, Fitness -----------------------------------------------------------------
#gives fitness at each date as a function of:
#airT.lt.lo = hours of air temperature below the low threshold each day (e.g., freezing)
#airT.dhh.5dH = degree heating hours above upper threshold (e.g., 35C) in 5 days
#bioT.mean = mean daily "biological" temperature (based on IR)
#par.cum = cumulative daily par
#swc.mean = mean daily soil water content
gr.Taa <- .2
gr.Tmin <- 15
gr.Tmax <- 35
gr.PARb <- .0002  #rate of decline in growth multiplier w increasing PAR (growht increases with PAR to max = 1)
gr.swcA <- 500  #logistic curve for effect of swc
gr.swcB <- 120
#survival depends on cumulative cold or heat stress
s.LTa <- 0.45  #rate of decline of surv with increase hours below freezing
s.HTa <- 0.07  #rate of decline of surv with increase DHH>35 in 5d
season = 15

dd<- d %>%
#  filter(date%within%interval(ymd("2018-1-01"),ymd("2018-10-31")) , site =="ABBY"|site=="MOAB") %>%
  filter(date%within%interval(ymd("2018-1-01"),ymd("2018-07-31")) , site =="ABBY"|site=="MOAB") %>%
  mutate(growth.T = ifelse(between(bioT.mean,15,35),gr.Taa*bioT.mean*(bioT.mean-gr.Tmin)*(gr.Tmax-bioT.mean)^(1/2),0),
         growth.PAR = (1-exp(-gr.PARb*par.cum)),
         growth.SWC = (1/(1 + gr.swcA*exp(-gr.swcB*swc.mean))),
         growth.day = growth.T/max(growth.T,na.rm=TRUE) * growth.PAR/max(growth.PAR,na.rm=TRUE) * growth.SWC/max(growth.SWC,na.rm=TRUE)) %>% 
  mutate(surv =  1*exp(-s.LTa*airT.lt.lo)*exp(-s.HTa*airT.dhh.5dH),
         surv.n = surv/max(surv,na.rm=TRUE),
         fit.day = growth.day * surv) %>% 
  group_by(site) %>% 
  mutate(surv.start=rollapply(surv, season, prod, na.rm=TRUE, fill=NA, align="left"),
         growth.start = rollapply(growth.day,season,sum,na.rm=TRUE,fill=NA,align="left"),
         fit.start = ifelse(surv.start*growth.start<1e-30,0,surv.start*growth.start),
         fit.start.n = fit.start/max(fit.start,na.rm=TRUE))
# Add Cue Function -----------------------------------------------------

#Degree chilling hours over 90 days - airT.dch.90d
#Degree heating hours avobe freezing over 30 days - airT.dhh.30dL
# Need heating hours to come *after* chilling hour threshold is met

dd %>% 
  ggplot(aes(x=date,y=airT.dch.90d, color=site)) +
  geom_line(size=1.5)+
  labs(x="", y="Degree Chilling Hours") +
  scale_x_date(date_labels = "%b", breaks=xax)

dd %>% 
  ggplot(aes(x=date,y=airT.dhh.30dL, color=site)) +
  geom_line(size=1.5)+
  labs(x="", y="Degree Heating Hours") +
  scale_x_date(date_labels = "%b", breaks=xax)

#Create probability thresholds for degree heating * degree cooling
dh.a = 3000
dh.b = 1/1000
dc.a = 500
dc.b = 15 #scale by thousands of cooling hours

dd<-dd %>%
  group_by(site) %>% 
  mutate(airT.dch.90d.n = airT.dch.90d/max(airT.dch.90d,na.rm=TRUE),
         airT.dhh.30dL.n = airT.dhh.30dL/max(airT.dhh.30dL,na.rm=TRUE),
         pEvent.dc = 1/(1 + dc.a*exp(-dc.b*airT.dch.90d.n)),
         pEvent.dh = 1/(1 + dh.a*exp(-dh.b*airT.dhh.30dL)),
         pEvent.day = (pEvent.dc)*pEvent.dh,
         pNoEvent.day = 1-pEvent.day,
         pEvent.by = 1-rollapply(pNoEvent.day,30,prod, na.rm=TRUE,partial=TRUE),
         pEvent.by2 = cummax(pEvent.by))


# Measurement Filter ------------------------------------------------------

dd<-dd %>% 
  group_by(site) %>% 
  mutate(airT.30d = rollapply(airT.mean,30,mean,na.rm=TRUE, partial=TRUE,align="right"),
         airT.max.30d = cummax(airT.30d)) 

# Figure parameters

q <- brewer.pal(n=12,name="Paired")
site.color.12 <- q[3:4] #c("tan3","chocolate4") #c("sienna1","sienna4")
site.color.precip <- q[1:2] #c("skyblue1","royalblue4")
site.color.PAR <- q[7:8] #c("gold","orange")#c("darkgoldenrod3","darkgoldenrod1")
site.color.airT <- q[5:6] #c("pink1","red3")#c("lightpink","maroon")
site.names12 <- c("Site 1", "Site 2")

theme_fig3 <-   theme_minimal() +
                theme(axis.ticks.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.text.x = element_text(size=8),
#                  axis.text.x = element_text(margin=margin(t=-10)),
                  axis.title = element_text(size=9),
                  legend.position = c(0.15,0.8),
                  legend.text = element_text(size=8),
                  legend.title = element_blank(),
#                  legend.direction = "horizontal",
                  panel.grid = element_line(color="white"),
                  panel.grid.major = element_line(size=0.25),
                  panel.grid.minor = element_line(size=0.25),
                  panel.background = element_rect(fill = "transparent",colour = NA),
                  plot.background = element_rect(fill = "transparent",colour = NA))


# Figure - Fitness ----------------------------------------------------------


Fig.Fitness.main <- 
  dd %>% 
  ggplot(aes(x=date,color=site)) +
  geom_line(aes(y=fit.start, color=site),size=1.2) +
  # geom_line(aes(y=surv.start, color=site),size=1,linetype="dashed")+
  # geom_line(aes(y=growth.start, color=site),size=1,linetype="dotted")+
  labs(x="", y="Fitness") +
  scale_x_date(date_labels = xaxL, breaks=xax)+
  scale_color_manual(values=site.color.12, labels = site.names12)+
  theme_fig3
Fig.Fitness.main

#for inset
x = rnorm(60,mean=0,1)
y = x+c(rnorm(30,mean=0,sd=.5),rnorm(30,mean=0,sd=1.7))
nom <- c(rep(1,30),rep(2,30))
tt <- data.frame(x=x,y=y,nom=as.factor(nom))
ymax.inset = 0.9*max(dd$fit.start,na.rm="TRUE")
ymin.inset = 0.4*max(dd$fit.start,na.rm="TRUE")

Fig.Fitness.inset <-
  tt %>% 
  ggplot(aes(x=x, y=y)) +
  geom_point(aes(color=nom),size=1)+
  geom_abline(intercept=0,slope=1,color="dark gray",size=1,alpha=0.5)+
  labs(x="Date of Event", y="Fitness") +
  scale_color_manual(values=site.color.12, labels = site.names12) +
  theme_fig3 +
  theme(legend.position = "none",
        axis.text.x = element_blank())
Fig.Fitness.inset

Fig.Fitness <- Fig.Fitness.main + 
  annotation_custom(
    ggplotGrob(Fig.Fitness.inset),
    xmin=ymd("2018-01-01"),xmax=ymd("2018-03-31"),ymin=ymin.inset,ymax=ymax.inset)

Fig.Fitness

# Figure - Event Cue ------------------------------------------------------
# geom_line(data = mdt[name == "A"], col = "#ff5a32", size = 2) 
Fig.EventCue <-
  dd %>%
  ggplot(aes(x=date,color=site)) +
  geom_line(aes(y=pEvent.by2), size=1.2) +
  geom_line(aes(y=0.8*airT.dch.90d.nAll),size=1.2, alpha=0.35,linetype="dashed") +
  geom_line(aes(y=0.8 *airT.dhh.30dL.nAll),size=1.2, alpha=0.35) +
  geom_segment(aes(x=ymd("2018-03-24"),xend=ymd("2018-04-10"),y=0.5,yend=0.5),size=0.7,color="gray40",arrow=arrow(length = unit(7,"pt"),ends="both"))+
  #geom_line(size=1)+
  labs(x="", y="Probability of Event") +
  theme_fig3+
  scale_x_date(date_labels = xaxL, breaks=xax) +
  scale_color_manual(values=site.color.12, labels = site.names12)
Fig.EventCue

# Figure - Measurement Filter ---------------------------------------------

Fig.Measurement <-
  dd %>% 
  mutate(airT.exceed = 1*cummax(airT.30d)>7.5) %>% 
  ggplot(aes(x=date, color=site)) +
  geom_line(aes(y=airT.30d), size=1.2,alpha=0.35) +
  geom_line(aes(y=airT.exceed*25),size=1.2) +
  geom_segment(aes(x=ymd("2018-04-8"),xend=ymd("2018-04-25"),y=7.5,yend=7.5),size=0.7,color="gray40",arrow=arrow(length = unit(7,"pt"),ends="both"))+
  # geom_line(aes(y=dayL),size=1.2)+
#  scale_y_continuous(sec.axis = sec_axis(~./25)) +
  labs(x="", y="Air Temperature (30-day) > Threshold") +
  theme_fig3 +
  scale_x_date(date_labels = xaxL, breaks=xax) +
  scale_color_manual(values=site.color.12, labels = site.names12)
Fig.Measurement

# #Plots of Environment in focal year ---------------------------------------------------


Fig.DailyAirT <-
  dd %>% 
  ggplot(aes(x=date, y=airT.mean,color=site)) +
  geom_line(size=1)+
  labs(x="",y="Temperature") +
  theme_fig3 +
  scale_x_date(date_labels = xaxL, breaks=xax) +
  scale_color_manual(values=site.color.airT, labels = site.names12)
Fig.DailyAirT

Fig.PAR7d <-
  dd %>% 
  ggplot(aes(x=date,y=par.7d,color=site)) +
  geom_line(size=1) +
  labs(x="", y="PAR") +
  theme_fig3 +
  scale_x_date(date_labels = xaxL, breaks=xax) +
  scale_color_manual(values=site.color.PAR, labels = site.names12)
Fig.PAR7d

Fig.Precip <-
  dd %>% 
  ggplot(aes(x=date, y=precip.daily,color=site)) +
  geom_line(size=1)+
  labs(x="", y="Precipitation") +
  theme_fig3 +
  scale_x_date(date_labels = xaxL, breaks=xax) +
  scale_color_manual(values=site.color.precip, labels = site.names12)
Fig.Precip

Fig.Daylength <-
  dd %>% 
  ggplot(aes(x=date, y=dayL,color=site)) +
  geom_line(size=1)+
  labs(x="", y="Day Length") +
  theme_minimal() +
  theme(axis.text.y = element_blank()) +
  theme(legend.position = "none") +
  scale_x_date(date_labels = xaxL, breaks=xax) +
  scale_color_brewer(palette="Dark2")
Fig.Daylength

dd %>% 
  ggplot +
  geom_line(aes(x=date,y=swc.mean,color=site)) +
  labs(x="Day of Year",y="Mean Soil Moisture")

dd %>% 
  ggplot +
  geom_line(aes(x=date,y=precip.cum.15d, color=site))+
  labs(x="Day of Year",y="Cum. Precipitation (15d)")

dd %>% 
  ggplot +
  geom_line(aes(x=date,y=airT.lt.lo,color=site)) +
  labs(x="Day of Year",y="Hours below Freezing")


# Assemble Figure ---------------------------------------------------------


Fig.Concept <- grid.arrange(Fig.DailyAirT,
             Fig.PAR7d, 
             Fig.Precip,
             Fig.Measurement,
             Fig.EventCue,
             Fig.Fitness,
             layout_matrix=cbind(c(1,2,3),c(4,5,6)))

ggsave(filename="graphs/conceptual/Fig_ConceptTrack.pdf",plot=Fig.Concept, width = 7, height= 7, units="in")

Fig.Concept.RHS <- grid.arrange(Fig.Measurement,
                                Fig.EventCue,
                                Fig.Fitness,
                                layout_matrix=cbind(c(1,2,3)))
Fig.Concept.LHS <- grid.arrange(Fig.DailyAirT,
                                Fig.PAR7d, 
                                Fig.Precip,
                                layout_matrix=cbind(c(1,2,3)))
ggsave(filename="graphs/conceptual/Fig_ConceptTrack.RHS.pdf", plot=Fig.Concept.RHS,width = 4.25, height=7,units="in")
ggsave(filename="graphs/conceptual/Fig_ConceptTrack.LHS.pdf", plot=Fig.Concept.LHS,width = 4.25, height=7,units="in")


#png for each panel
ggsave(filename="graphs/conceptual/Fig_ConceptTrack.LHS.svg", plot=Fig.Concept.LHS,width = 61, height=141,units="mm",bg="transparent")
ggsave(filename="graphs/conceptual/Fig_ConceptTrack.Fitness.main.svg", plot=Fig.Fitness.main,width = 100, height=70,units="mm",bg="transparent")
ggsave(filename="graphs/conceptual/Fig_ConceptTrack.Fitness.inset.svg", plot=Fig.Fitness.inset,width = 42, height=32,units="mm",bg="transparent")
ggsave(filename="graphs/conceptual/Fig_ConceptTrack.Fitness.svg", plot=Fig.Fitness,width = 100, height=70,units="mm",bg="transparent")
ggsave(filename="graphs/conceptual/Fig_ConceptTrack.EventCue.svg", plot=Fig.EventCue,width = 100, height=70,units="mm",bg="transparent")
ggsave(filename="graphs/conceptual/Fig_ConceptTrack.Measurement.svg", plot=Fig.Measurement,width = 100, height=70,units="mm",bg="transparent")
ggsave(filename="graphs/conceptual/Fig_ConceptTrack.airT.svg", plot=Fig.DailyAirT,width = 61, height=55,units="mm",bg="transparent")
ggsave(filename="graphs/conceptual/Fig_ConceptTrack.precip.svg", plot=Fig.Precip,width = 61, height=55,units="mm",bg="transparent")
ggsave(filename="graphs/conceptual/Fig_ConceptTrack.PAR.svg", plot=Fig.PAR7d,width = 61, height=55,units="mm",bg="transparent")




