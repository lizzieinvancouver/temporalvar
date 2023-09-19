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
nruns<-3000
g_notxi <- 0
c_warm<- 0.5
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
   source(here("R","sourcefiles","PriEff_Output.R"))
 
}else {
  source(here("R","sourcefiles","PriEff_OutputwFract.R"))
}
}
  
  sapply(g_cumulative, tail, 1) ## Dan think this should be the realized germinaiton fraction
#9:44 am 
#save.image("firstrun.Rda")
### Thought is the tradeoff about phenology (frost risk) acrually in this model in any way?
#####now plot and work with data
check<-read.csv("R/output/prieff_params.csv")

colnames(check)

check$g_time_dif<-abs(check$sp1_mean_tau_g50-check$sp2_mean_tau_g50)
check$g_time_dif

quantile(check$xi.mu)

check$ave_chill<-ifelse(check$xi.mu>2,"12 weeks","6 weeks")
p1<-ggplot(check,aes(g_time_dif))+geom_density(aes(fill=ave_chill),position="identity",alpha=0.4)+ggthemes::theme_base()+
  scale_fill_viridis_d()

check$`R1/R2`<-check$sp1_Rstar/check$sp2_Rstar
check$`sen1/sen2`<-check$sp1_xi_tau/check$sp2_xi_tau
check$`t501/t502`<-check$sp1_mean_tau_g50/check$sp2_mean_tau_g50
check$logR1R2<-log(check$`R1/R2`)
check$logsens1sens2<-log(check$`sen1/sen2`)
check$logt501t502<-log(check$`t501/t502`)


p2<-ggplot(check,aes(logsens1sens2,logt501t502))+geom_point(size=0.2)+facet_wrap(~ave_chill)+
  ggthemes::theme_base()
jpeg("plots/coexistance_chilldiffs.jpeg")
ggpubr::ggarrange(p1,p2,ncol=1)
dev.off()

if(FALSE){
check$trial<-NA
check$trial[which(check$sp1_meangmax==0.8 & check$sp2_meangmax==0.8)]<-"time only"

check$trial[which(check$sp1_meangmax!=0.8 & check$sp2_meangmax!=0.8 & check$sp1_SDgmax==0 & check$sp2_SDgmax==0)]<- "both fixed"

check$trial[which(check$sp1_meangmax!=0.8 & check$sp2_meangmax!=0.8 & check$sp1_SDgmax!=0 & check$sp2_SDgmax==0)]<- "sp2 fixed"

check$trial[which(check$sp1_meangmax!=0.8 & check$sp2_meangmax!=0.8 & check$sp1_SDgmax==0 &check$sp2_SDgmax!=0)]<- "sp1 fixed"

check$trial[which(check$sp1_meangmax!=0.8 & check$sp2_meangmax!=0.8
                  & check$sp1_SDgmax!=0 & check$sp2_SDgmax!=0)]<- "varying"
}







check$coexist<-NA
check$coexist[which(check$sp1_ex==0 & check$sp2_ex==0)]<-"both extinct"
check$coexist[which(check$sp1_ex!=0 & check$sp2_ex==0)]<-"sp1 win"
check$coexist[which(check$sp1_ex==0 & check$sp2_ex!=0)]<-"sp2 win"
check$coexist[which(check$sp1_ex!=0 & check$sp2_ex!=0)]<-"coexist"




if(FALSE){
ggpubr::ggarrange(ggplot(check,aes(sp1_meangmax))+geom_histogram(),
ggplot(check,aes(sp2_meangmax))+geom_histogram())

ggpubr::ggarrange(ggplot(check,aes(sp1_mean_tau_g50))+geom_histogram(),
                  ggplot(check,aes(sp2_mean_tau_g50))+geom_histogram())

ggpubr::ggarrange(ggplot(check,aes(sp1_Rstar))+geom_histogram(),
                  ggplot(check,aes(sp2_Rstar))+geom_histogram())

ggpubr::ggarrange(ggplot(check,aes(sp1_xi100))+geom_histogram(),
                  ggplot(check,aes(sp2_xi100))+geom_histogram())


}

#jpeg("plots/coEx_prieff.jpeg")
ggplot(check,aes(logsens1sens2,logR1R2))+geom_point(aes(color=coexist),size=1)+geom_vline(xintercept = 0)+
  geom_hline(yintercept=0)+facet_wrap(~ave_chill)
        
table(check$ave_chill)

ggplot(check,aes(logsens1sens2,logR1R2))+
  geom_point(aes(color=coexist),size=1)+geom_vline(xintercept = 0)+geom_hline(yintercept=0)+
  facet_grid(ave_chill~coexist)+geom_smooth(method="lm")
dev.off()

check$timecat<-NA
check$timecat[which(check$g_time_dif<=5)]<-".0-5"
check$timecat[which(check$g_time_dif>5 & check$g_time_dif<=10)]<-".5-10"
check$timecat[which(check$g_time_dif>10 & check$g_time_dif<=15)]<-"10-15"
check$timecat[which(check$g_time_dif>15 & check$g_time_dif<=20)]<-"15-20"
check$timecat[which(check$g_time_dif>20 )]<-"20+"
check2<-filter(check,coexist=="coexist")


table(check2$ave_chill)
table(check$ave_chill)
lm(check2$logR1R2~check2$logsens1sens2*check2$ave_chill)

round(65/1513,2)
round(49/1487,2)
jpeg("plots/coexistance_runner.jpeg")
ggplot(check2,aes(logsens1sens2,logR1R2))+
  geom_point(aes(color=timecat),size=3)+geom_smooth(method="lm",fullrange=TRUE)+geom_vline(xintercept = 0)+geom_hline(yintercept=0)+
  facet_wrap(~ave_chill)+ggthemes::theme_base()+scale_color_viridis_d()
dev.off()

ggplot(check2,aes(logt501t502,logR1R2))+
  geom_point(aes(color=timecat),size=2)+geom_smooth(method="lm")+geom_vline(xintercept = 0)+geom_hline(yintercept=0)+
  facet_wrap(~ave_chill)+ggthemes::theme_base()



dev.off()

check3<-filter(check,coexist=="sp1 win" & logR1R2>0|coexist=="sp2 win" & logR1R2<0)
ggplot(check3,aes(logt501t502,logR1R2))+
  geom_point(aes(color=timecat,shape=coexist),size=2)+geom_vline(xintercept = 0)+geom_hline(yintercept=0)+
  facet_grid(~ave_chill)+ggthemes::theme_base()

ggplot(check3,aes(logsens1sens2,logR1R2))+
  geom_point(aes(color=timecat,shape=coexist),size=2)+geom_vline(xintercept = 0)+geom_hline(yintercept=0)+
  facet_grid(~ave_chill)+ggthemes::theme_base()+scale_color_viridis_d()

jpeg("plots/dominance.jpeg")
ggplot()+
  geom_point(data=check,aes(x=logsens1sens2,y=logR1R2),size=1,color="grey")+
  geom_point(data=check3,aes(x=logsens1sens2,y=logR1R2,color=timecat),size=3)+
  facet_wrap(~ave_chill)+ggthemes::theme_base()+
  geom_vline(xintercept = 0)+geom_hline(yintercept=0)+scale_color_viridis_d()
dev.off()
check2 %>%group_by(ave_chill) %>% count()
check3 %>%group_by(ave_chill) %>% count()

unique(check2$xi.mu)
                         
ggplot(check,aes(logsens1sens2,logR1R2))+
  geom_point(aes(color=coexist),size=1)+geom_smooth(method="lm")+geom_vline(xintercept = 0)+geom_hline(yintercept=0)+
  facet_grid(trial~xi.mu)




table(check$ave_chill,check$coexist)
table(check3$ave_chill,check3$coexist)
57/720
53/723

87/724
93/706


#####with germ fract
fracking<-read.csv("R/output/prieff_params_withfractvar.csv")

fracking$coexist<-NA
fracking$coexist[which(fracking$sp1_ex==0 & fracking$sp2_ex==0)]<-"both extinct"
fracking$coexist[which(fracking$sp1_ex!=0 & fracking$sp2_ex==0)]<-"sp1 win"
fracking$coexist[which(fracking$sp1_ex==0 & fracking$sp2_ex!=0)]<-"sp2 win"
fracking$coexist[which(fracking$sp1_ex!=0 & fracking$sp2_ex!=0)]<-"coexist"

fracking$`R1/R2`<-fracking$sp1_Rstar/fracking$sp2_Rstar
fracking$`sen1/sen2`<-fracking$sp1_xi_tau/fracking$sp2_xi_tau
fracking$`t501/t502`<-fracking$sp1_mean_tau_g50/fracking$sp2_mean_tau_g50
fracking$`xi0_1/xi0_2`<-fracking$sp1_xi0/fracking$sp2_xi0
fracking$`xi100_1/xi100_2`<-fracking$sp1_xi100/fracking$sp2_xi100
fracking$`xi_rng_1/xi_rng_2`<-fracking$sp1_xi_rng/fracking$sp2_xi_rng

fracking$logR1R2<-log(fracking$`R1/R2`)
fracking$logsens1sens2<-log(fracking$`sen1/sen2`)
fracking$logt501t502<-log(fracking$`t501/t502`)
fracking$`logxi0_1/xi0_2`<-log(fracking$`xi0_1/xi0_2`)
fracking$`logxi100_1/xi100_2`<-log(fracking$`xi100_1/xi100_2`)
fracking$`logxi_rng_1/xi_rng_2`<-log(fracking$`xi_rng_1/xi_rng_2`)

fracking$ave_chill<-ifelse(fracking$xi.mu>2,"12 weeks","6 weeks")

ggplot(fracking,aes(logsens1sens2,logR1R2))+
  geom_point(aes(color=coexist),size=3)+geom_vline(xintercept = 0)+geom_hline(yintercept=0)+
  facet_wrap(~ave_chill)+ggthemes::theme_base()+scale_color_viridis_d()

ggplot(fracking,aes(logsens1sens2,`logxi0_1/xi0_2`))+
  geom_point(aes(color=coexist),size=3)+geom_vline(xintercept = 0)+geom_hline(yintercept=0)+
  facet_wrap(~ave_chill)+ggthemes::theme_base()+scale_color_viridis_d()

