#Priority Effects - main code
rm(list=ls()) 
#libraries and working directories, oh my!
library(deSolve)
library(scales)
require(MultiRNG)
if(length(grep("lizzie", getwd())>0)) {
    setwd("~/Documents/git/projects/temporalvar/R")
}

if(length(grep("danielbuonaiuto", getwd())>0)) {
  setwd("~/Documents/git/temporalvar/")
}
library(here)
here()
#2:38 pm7:29 am 
load("firstrun.Rda")

#define the run - Consider creating a dataframe with combinations of parms to test
nruns<-1:2000
outputy<-data.frame()
for (n in c(1,nruns)) {

  nyrs <- 500
  
  #define the environment for this run
  source(here("R","sourcefiles","PriEff_Envt.R"))
  
  #define the species in this run
  source(here("R","sourcefiles","PriEff_Species.R"))
  
  #run the model for nyrs
  source(here("R","sourcefiles","PriEff_Comp.R"))
  source(here("R","sourcefiles","PriEff_Model.R"))

  # #write out the results for this run
  # source(here("R","sourcefiles","PriEff_Output.R"))
  mydat<-data.frame(sp1_Rstar=Rstar[1],sp2_Rstar=Rstar[2],
                    sp1_gmax=mean(gmax[,1]),sp1_gmax_sd=sd(gmax[,1]),
                    sp2_gmax=mean(gmax[,2]),sp2_gmax_sd=sd(gmax[,2]),
                    sp1_sense=xi_tau[1],sp2_sense=xi_tau[2],
                    sp1_t50=mean(tau_g50[,1]),sp1_t50_sd=sd(tau_g50[,1]),
                    sp2_t50=mean(tau_g50[,2]),sp2_t50_sd=sd(tau_g50[,2]), 
                    sp1.Bfin=Bfin[500,1],sp2.Bfin=Bfin[500,2])
  outputy<-rbind(mydat,outputy)
}
#2:38 pm7:29 am 
save.image("firstrun.Rda")
outputy$trial<-NA
outputy$trial[which(outputy$sp1_gmax==0.8 & outputy$sp2_gmax==0.8)]<-"time only"
outputy$trial[which(outputy$sp1_gmax_sd==0.0 & outputy$sp2_gmax_sd!=0.0)]<-"sp1 fixed"
outputy$trial[which(outputy$sp1_gmax_sd!=0.0 & outputy$sp2_gmax_sd==0.0)]<-"sp2 fixed"
outputy$trial[which(outputy$sp1_gmax_sd!=0.0 & outputy$sp2_gmax_sd!=0.0)]<-"both sensitive"
outputy$trial[which(outputy$sp1_gmax_sd==0.0& outputy$sp1_gmax!=0.8 & outputy$sp2_gmax_sd==0.0 & outputy$sp2_gmax!=0.8)]<-"both fixed-different fract"


outputy$coexistence<-NA
outputy$coexistence[which(outputy$sp1.Bfin==0 & outputy$sp2.Bfin==0)]<-"both extinct"
outputy$coexistence[which(outputy$sp1.Bfin!=0 & outputy$sp2.Bfin==0)]<-"sp1 win"
outputy$coexistence[which(outputy$sp1.Bfin==0 & outputy$sp2.Bfin!=0)]<-"sp2 win"
outputy$coexistence[which(outputy$sp1.Bfin!=0 & outputy$sp2.Bfin!=0)]<-"coexistence"

outputy$sp1_fin<-NA
outputy$sp2_fin<-NA
outputy$sp1_fin<-ifelse(outputy$sp1.Bfin==0,0,1)
outputy$sp2_fin<-ifelse(outputy$sp2.Bfin==0,0,1)
outputy$happy_together<-outputy$sp1_fin+outputy$sp2_fin

outputy$`R1/R2`<-outputy$sp1_Rstar/outputy$sp2_Rstar
outputy$`sen1/sen2`<-outputy$sp1_sense/outputy$sp2_sense
outputy$`T1/T2`<-outputy$sp1_t50/outputy$sp2_t50
outputy$`gmax1/gmax2`<-outputy$sp1_gmax/outputy$sp2_gmax

outputy$logR1R2<-log(outputy$`R1/R2`)
outputy$logsens1sens2<-log(outputy$`sen1/sen2`)
outputy$loggmax1gmax2<-log(outputy$`gmax1/gmax2`)

library(ggplot2)
outputy2<-dplyr::filter(outputy,trial!="both fixed-different fract")


jpeg("plots/firstrun_plots.jpeg",height=5,width=7, unit="in", res=300)
ggplot(outputy2,aes(logsens1sens2,logR1R2))+
  geom_point(aes(color=happy_together),size=1)+#ylim(0,10)+xlim(0,10)+
  facet_wrap(~trial)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+ggthemes::theme_few()+
  scale_color_viridis_c(option = "C")
dev.off()

ggplot(outputy2,aes(logsens1sens2,logR1R2))+
  geom_point(aes(color=as.factor(happy_together)),size=1)+#ylim(0,10)+xlim(0,10)+
  facet_wrap(~trial)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+ggthemes::theme_few()+
  scale_color_viridis_d(option = "C")

#ggplot(outputy2,aes(`sen1/sen2`,`R1/R2`))+
 # geom_point(aes(color=coexistence),size=1)+ylim(0,10)+xlim(0,10)+
  #facet_wrap(~trial)+geom_hline(yintercept=1)+geom_vline(xintercept=1)+ggthemes::theme_few()+
  #scale_color_viridis_d(option = "C")
#dev.off()

jpeg("plots/firstrun_plotsgmax.jpeg",height=5,width=7, unit="in", res=300)
ggplot(outputy2,aes(loggmax1gmax2,logR1R2))+
  geom_point(aes(color=happy_together),size=1)+#ylim(-0,10)+xlim(0,2)+
  facet_wrap(~trial)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+ggthemes::theme_few()+
  scale_color_viridis_c(option = "C")
dev.off()


jpeg("plots/firstrun_plots_altview.jpeg",height=5,width=7, unit="in", res=300)
ggplot(outputy,aes(logsens1sens2,logR1R2))+
  geom_point(aes(color=as.factor(happy_together)),size=1)+#ylim(0,50)+xlim(0,50)+
  facet_grid(happy_together~trial)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+ggthemes::theme_few()+
  scale_color_viridis_d(option = "C")
dev.off()

ggplot(outputy2,aes(`T1/T2`,`R1/R2`))+
  geom_point(aes(color=coexistence),size=1)+ylim(0,10)+xlim(0,10)+
  facet_wrap(~trial)+geom_hline(yintercept=1)+geom_vline(xintercept=1)+ggthemes::theme_few()+
  scale_color_viridis_d(option = "C")

ggplot(outputy,aes(`sen1/sen2`,`R2/R1`))+
  geom_point(aes(color=coexistence,shape=trial),size=1)+ylim(0,10)+xlim(0,10)+
  geom_hline(yintercept=1)+geom_vline(xintercept=1)+ggthemes::theme_few()+scale_color_viridis_d(option = "C")

ggplot(outputy,aes(`T1/T2`,`R2/R1`))+
  geom_point(aes(color=coexistence,shape=trial),size=1)+ylim(0,10)+xlim(0,10)+
  geom_hline(yintercept=1)+geom_vline(xintercept=1)+ggthemes::theme_few()+scale_color_viridis_d(option = "C")


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


