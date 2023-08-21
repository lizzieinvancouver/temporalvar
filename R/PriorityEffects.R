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
#2:38 pm7:29 am 
#load("firstrun.Rda")

dev.off()

#define the run - Consider creating a dataframe with combinations of parms to test
nruns<-1
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
   source(here("R","sourcefiles","PriEff_Output.R"))
 
}


sapply(g_cumulative, tail, 1) ## Dan think this should be the realized germinaiton fraction
#9:44 am 
#save.image("firstrun.Rda")
### Thought is the tradeoff about phenology (frost risk) acrually in this model in any way?
#####now plot and work with data
check<-read.csv("R/output/prieff_params.csv")
colnames(check)
check$trial<-NA
check$trial[which(check$sp1_meangmax==0.8 & check$sp2_meangmax==0.8)]<-"time only"

check$trial[which(check$sp1_meangmax!=0.8 & check$sp2_meangmax!=0.8 & check$sp1_SDgmax==0 & check$sp2_SDgmax==0)]<- "both fixed"

check$trial[which(check$sp1_meangmax!=0.8 & check$sp2_meangmax!=0.8 & check$sp1_SDgmax!=0 & check$sp2_SDgmax==0)]<- "sp2 fixed"

check$trial[which(check$sp1_meangmax!=0.8 & check$sp2_meangmax!=0.8 & check$sp1_SDgmax==0 &check$sp2_SDgmax!=0)]<- "sp1 fixed"

check$trial[which(check$sp1_meangmax!=0.8 & check$sp2_meangmax!=0.8
                  & check$sp1_SDgmax!=0 & check$sp2_SDgmax!=0)]<- "varying"





library(ggplot2)


check$dif<-check$sp1_xi_tau-check$sp2_xi_tau
ggpubr::ggarrange(ggplot(check,aes(sp1_xi_tau))+geom_histogram(bins=30),
ggplot(check,aes(sp2_xi_tau))+geom_histogram(bins=30),
ggplot(check,aes(dif))+geom_histogram(bins=30))





check$coexist<-NA
check$coexist[which(check$sp1_ex==0 & check$sp2_ex==0)]<-"both extinct"
check$coexist[which(check$sp1_ex!=0 & check$sp2_ex==0)]<-"sp1 win"
check$coexist[which(check$sp1_ex==0 & check$sp2_ex!=0)]<-"sp2 win"
check$coexist[which(check$sp1_ex!=0 & check$sp2_ex!=0)]<-"coexist"



check$`R1/R2`<-check$sp1_Rstar/check$sp2_Rstar
check$`sen1/sen2`<-check$sp1_xi_tau/check$sp2_xi_tau


check$logR1R2<-log(check$`R1/R2`)
check$logsens1sens2<-log(check$`sen1/sen2`)



jpeg("plots/coEx_prieff.jpeg")
ggpubr::ggarrange(ggplot(check,aes(logsens1sens2,logR1R2))+
  geom_point(aes(color=coexist),size=1)+geom_vline(xintercept = 0)+geom_hline(yintercept=0),
ggplot(check,aes(logsens1sens2,logR1R2))+
  geom_point(aes(color=coexist),size=1)+geom_vline(xintercept = 0)+geom_hline(yintercept=0)+
  facet_wrap(~coexist),common.legend=TRUE)

dev.off()

head(check)
ggpubr::ggarrange(ggplot(check,aes(sp1_meangmax))+geom_histogram(),
ggplot(check,aes(sp2_meangmax))+geom_histogram())

ggpubr::ggarrange(ggplot(check,aes(sp1_mean_tau_g50))+geom_histogram(),
                  ggplot(check,aes(sp2_mean_tau_g50))+geom_histogram())

ggpubr::ggarrange(ggplot(check,aes(sp1_Rstar))+geom_histogram(),
                  ggplot(check,aes(sp2_Rstar))+geom_histogram())

ggpubr::ggarrange(ggplot(check,aes(sp1_xi100))+geom_histogram(),
                  ggplot(check,aes(sp2_xi100))+geom_histogram())
check2<-filter(check,coexist=="coexist")
check2<-filter(check2,trial %in% c("time only", "varying"))

check2$ave_chill<-ifelse(check2$xi.mu>2,"8.5","4.5")
jpeg("plots/coexistance_runner.jpeg")

ggplot(check2,aes(logsens1sens2,logR1R2))+
  geom_point(aes(),size=1)+geom_smooth(method="lm")+geom_vline(xintercept = 0)+geom_hline(yintercept=0)+
  facet_grid(trial~ave_chill)+ggthemes::theme_base()
dev.off()

check2 %>%group_by(trial,ave_chill) %>% count()

unique(check2$xi.mu)
                         
ggplot(check,aes(logsens1sens2,logR1R2))+
  geom_point(aes(color=coexist),size=1)+geom_smooth(method="lm")+geom_vline(xintercept = 0)+geom_hline(yintercept=0)+
  facet_grid(trial~xi.mu)



table(check2$trial,check2$xi.mu)
