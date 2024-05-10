#Priority Effects - main code
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

if(length(grep("danielbuonaiuto", getwd())>0)) {
  setwd("~/Documents/git/temporalvar/")
}
library(here)
here()


nruns<-1000
g_notxi <- 1
c_warm<- 1
makeplots <- FALSE
outputy<-data.frame()
#outputy2<-data.frame()
#outputy3<-data.frame()
for (j in c(1:nruns)) {
  
  nyrs <- 300
  #define the environment for this run
  source(here("R","sourcefiles","PriEff_Envt.R"))
  
  #define the species in this run
  source(here("R","sourcefiles","PriEff_Species.R"))
  
  #run the model for nyrs
  source(here("R","sourcefiles","PriEff_Comp.R"))
  source(here("R","sourcefiles","PriEff_Model.R"))
  
  # #write out the results for this run
  if (g_notxi==1) {
    source(here("R","sourcefiles","shiftoutput.R"))
    
  }else {
    source(here("R","sourcefiles","PriEff_OutputwFract.R"))
  }
}

check<-read.csv("R/output/coexshifts_params.csv")
check$coexist<-NA
check$coexist[which(check$sp1_ex==0 & check$sp2_ex==0)]<-"both extinct"
check$coexist[which(check$sp1_ex!=0 & check$sp2_ex==0)]<-"sp1 win"
check$coexist[which(check$sp1_ex==0 & check$sp2_ex!=0)]<-"sp2 win"
check$coexist[which(check$sp1_ex!=0 & check$sp2_ex!=0)]<-"coexist"

coex<-filter(check, coexist=="coexist")



c_warm<- 0
outputy<-data.frame()


for (j in c(1:nruns)) {
cox<-check[j,]

  source(here("R","sourcefiles","PriEff_Envt.R"))
  
  #define the species in this run
  source(here("R","sourcefiles","PriEf_Species_shifter.R"))
  ##make new version of this and comment out xi.tau and Rstar
  
  #run the model for nyrs
  source(here("R","sourcefiles","PriEff_Comp.R"))
  source(here("R","sourcefiles","PriEff_Model.R"))
  
  # #write out the results for this run
  if (g_notxi==1) {
    source(here("R","sourcefiles","shiftoutput.R"))
    
  }else {
    source(here("R","sourcefiles","PriEff_OutputwFract.R"))
  }
}

d1<-read.csv("R/output/coexshifts_params.csv")
d2<-read.csv("R/output/coexshifts_params2.csv")
d1<-d1[match(d2$sp1_Rstar, d1$sp1_Rstar),]
d2$run<-1:1000
head(d1)
head(d2)
d1$coexist<-NA
d1$coexist[which(d1$sp1_ex==0 & d1$sp2_ex==0)]<-"both extinct"
d1$coexist[which(d1$sp1_ex!=0 & d1$sp2_ex==0)]<-"sp1 win"
d1$coexist[which(d1$sp1_ex==0 & d1$sp2_ex!=0)]<-"sp2 win"
d1$coexist[which(d1$sp1_ex!=0 & d1$sp2_ex!=0)]<-"coexist"

d2$coexist<-NA
d2$coexist[which(d2$sp1_ex==0 & d2$sp2_ex==0)]<-"both extinct"
d2$coexist[which(d2$sp1_ex!=0 & d2$sp2_ex==0)]<-"sp1 win"
d2$coexist[which(d2$sp1_ex==0 & d2$sp2_ex!=0)]<-"sp2 win"
d2$coexist[which(d2$sp1_ex!=0 & d2$sp2_ex!=0)]<-"coexist"

head(d1)
head(d2)

coex1<-filter(d1, coexist=="coexist")
die1<-filter(d2,run %in% coex1$run)

coex2<-filter(d2, coexist=="coexist")
die2<-filter(d1,run %in% coex2$run)



dd<-rbind(coex2,die2,coex1,die1)


dd$`sp1_Rstar-sp2_Rstar`<-dd$sp1_Rstar-dd$sp2_Rstar
dd$`sp1_xi_tau-sp2_xi_tau`<-dd$sp1_xi_tau-dd$sp2_xi_tau
dd$`sp1_mean_tau_g50-sp2_mean_tau_g50`<-dd$sp1_mean_tau_g50-dd$sp2_mean_tau_g50
dd$ave_chill<-ifelse(dd$xi.mu>2,"12 weeks","6 weeks")
dd$cox<-ifelse(dd$coexist=="coexist","coexistence","no coexistence")

dd2<-filter(dd,run %in% c(842,272))
dd.e<-filter(dd,cox=="coexistence")

ggplot(data=dd
       ,aes(x=`sp1_mean_tau_g50-sp2_mean_tau_g50`,y=`sp1_Rstar-sp2_Rstar`))+
  geom_point(aes(shape=ave_chill,color=ave_chill,size=cox))+
  scale_size_manual(values=c(3,1))+geom_line(aes(group = run),size=0.1)+ylim(-.5,.5)+coord_cartesian(ylim=c(-.1,0.1))+
  geom_smooth(data=dd.e,aes(x=`sp1_mean_tau_g50-sp2_mean_tau_g50`,y=`sp1_Rstar-sp2_Rstar`,color=ave_chill),method="lm",fullrange=TRUE)


ggplot(data=dd2
       ,aes(x=`sp1_mean_tau_g50-sp2_mean_tau_g50`,y=`sp1_Rstar-sp2_Rstar`))+
  geom_point(aes(shape=ave_chill,color=ave_chill,size=cox))+
  scale_size_manual(values=c(3,1))+geom_line(aes(group = run),size=0.1)

write.csv(dd2,"R/output/casepairs.csv") 





