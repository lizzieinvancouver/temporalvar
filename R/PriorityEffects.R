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
#fracking<-read.csv("R/output/prieff_params_withfractvar.csv")


check$g_time_dif<-abs(check$sp1_mean_tau_g50-check$sp2_mean_tau_g50)

#fracking$g_time_dif<-abs(fracking$sp1_mean_tau_g50-fracking$sp2_mean_tau_g50)
#fracking$fracdiff<-abs(fracking$sp1_meangmax-fracking$sp2_meangmax)


check$ave_chill<-ifelse(check$xi.mu>2,"12 weeks","6 weeks")
#fracking$ave_chill<-ifelse(fracking$xi.mu>2,"12 weeks","6 weeks")

check$`R1/R2`<-check$sp1_Rstar/check$sp2_Rstar
check$`sen1/sen2`<-check$sp1_xi_tau/check$sp2_xi_tau
check$`t501/t502`<-check$sp1_mean_tau_g50/check$sp2_mean_tau_g50
check$logR1R2<-log(check$`R1/R2`)
check$logsens1sens2<-log(check$`sen1/sen2`)
check$logt501t502<-log(check$`t501/t502`)


check$coexist<-NA
check$coexist[which(check$sp1_ex==0 & check$sp2_ex==0)]<-"both extinct"
check$coexist[which(check$sp1_ex!=0 & check$sp2_ex==0)]<-"sp1 win"
check$coexist[which(check$sp1_ex==0 & check$sp2_ex!=0)]<-"sp2 win"
check$coexist[which(check$sp1_ex!=0 & check$sp2_ex!=0)]<-"coexist"

#fracking$coexist<-NA
#fracking$coexist[which(fracking$sp1_ex==0 & fracking$sp2_ex==0)]<-"both extinct"
#fracking$coexist[which(fracking$sp1_ex!=0 & fracking$sp2_ex==0)]<-"sp1 win"
#fracking$coexist[which(fracking$sp1_ex==0 & fracking$sp2_ex!=0)]<-"sp2 win"
#fracking$coexist[which(fracking$sp1_ex!=0 & fracking$sp2_ex!=0)]<-"coexist"

#fracking$`R1/R2`<-fracking$sp1_Rstar/fracking$sp2_Rstar
#fracking$`sen1/sen2`<-fracking$sp1_xi_tau/fracking$sp2_xi_tau
#fracking$`t501/t502`<-fracking$sp1_mean_tau_g50/fracking$sp2_mean_tau_g50
#fracking$`xi0_1/xi0_2`<-fracking$sp1_xi0/fracking$sp2_xi0
#fracking$`xi100_1/xi100_2`<-fracking$sp1_xi100/fracking$sp2_xi100
#fracking$`xi_rng_1/xi_rng_2`<-fracking$sp1_xi_rng/fracking$sp2_xi_rng
#fracking$`gmax1/gmax2`<-fracking$sp1_meangmax/fracking$sp2_meangmax


#fracking$logR1R2<-log(fracking$`R1/R2`)
#fracking$logsens1sens2<-log(fracking$`sen1/sen2`)
#fracking$logt501t502<-log(fracking$`t501/t502`)
#fracking$`logxi0_1/xi0_2`<-log(fracking$`xi0_1/xi0_2`)
#fracking$`logxi100_1/xi100_2`<-log(fracking$`xi100_1/xi100_2`)
#fracking$`logxi_rng_1/xi_rng_2`<-log(fracking$`xi_rng_1/xi_rng_2`)
#fracking$`loggmax1/gmax2`<-log(fracking$`gmax1/gmax2`)

#fracking$timecat<-abs(fracking$logt501t502)
check$timecat<-abs(check$logt501t502)

#fracking2<-filter(fracking,coexist=="coexist")
check2<-filter(check,coexist=="coexist")

p1<-ggplot(check,aes(g_time_dif))+geom_density(aes(fill=ave_chill),position="identity",alpha=0.8)+ggthemes::theme_few()+
  scale_fill_viridis_d(name="average chilling period",begin=.15,end=.92)+xlab("difference in T50")+ylab("frequency")

p2<-ggplot(check,aes(logsens1sens2,logt501t502))+geom_point(size=.2,aes(color=ave_chill))+
  scale_color_viridis_d(name="average chilling period",begin=.15,end=.92)+#facet_wrap(~ave_chill)+
  ggthemes::theme_few()+ylab(bquote(log(T50[sp1]/T50[sp2])))+xlab(bquote(log(sensitivity[sp1]/sensitivity[sp2])))

#p3<-ggplot(fracking,aes(fracdiff))+geom_density(aes(fill=ave_chill),position="identity",alpha=0.4)+ggthemes::theme_base()+
 # scale_fill_viridis_d()

#p4<-ggplot(fracking,aes(`logxi_rng_1/xi_rng_2`,`loggmax1/gmax2`))+geom_point(size=0.2)+facet_wrap(~ave_chill)+
 # ggthemes::theme_base()
jpeg("plots/coexistance_chilldiffs.jpeg",width = 9,height=5,unit='in',res=200)
ggpubr::ggarrange(p1,p2,ncol=2,nrow=1,common.legend=TRUE,labels = c("a)","b)"))
dev.off()


round(table(check2$ave_chill)/table(check$ave_chill),2)
#round(table(fracking2$ave_chill)/table(fracking$ave_chill),2)

summary(lm(check2$logR1R2~check2$logsens1sens2*check2$ave_chill))
summary(lm(check3$logR1R2~check3$logsens1sens2*check3$ave_chill))

jpeg("plots/coexistance_runner.jpeg",width = 8,height=10,unit='in',res=200)

  
  
coex<-ggplot()+
    geom_point(data=check,aes(x=logsens1sens2,y=logR1R2),size=.01,color="grey")+
    geom_point(data=check2,aes(logsens1sens2,logR1R2,color=ave_chill),size=2)+geom_smooth(data=check2,method="lm",fullrange=TRUE,aes(logsens1sens2,logR1R2,color=ave_chill,fill=ave_chill))+geom_vline(xintercept = 0)+geom_hline(yintercept=0)+
  #facet_wrap(~ave_chill)
    ggthemes::theme_few()+scale_color_viridis_d(name="average chilling period",begin=.15,end=.92)+scale_fill_viridis_d(name="average chilling period",begin=.15,end=.92)+
    xlab(bquote(log(sensitivity[sp1]/sensitivity[sp2])))+
  ylab(bquote(log(Rstar[sp1]/Rstar[sp2])))+theme(legend.position = "top")




check3<-filter(check,coexist=="sp1 win" & logR1R2>0|coexist=="sp2 win" & logR1R2<0)

round(table(check3$ave_chill)/table(check$ave_chill),2)
round(table(fracking3$ave_chill)/table(fracking$ave_chill),2)


jpeg("plots/dominance.jpeg")
dom<-ggplot()+
  geom_point(data=check,aes(x=logsens1sens2,y=logR1R2),size=.01,color="grey")+
  geom_point(data=check3,aes(x=logsens1sens2,y=logR1R2,color=ave_chill,fill=ave_chill),size=2)+
  ggthemes::theme_few()+
  geom_vline(xintercept = 0)+geom_hline(yintercept=0)+
  geom_smooth(data=check3,method="lm",fullrange=TRUE,aes(logsens1sens2,logR1R2,color=ave_chill,fill=ave_chill))+
  scale_color_viridis_d(name="average chilling period",begin=.15,end=.92)+scale_fill_viridis_d(name="average chilling period",begin=.15,end=.92)+
  xlab(bquote(log(sensitivity[sp1]/sensitivity[sp2])))+
  ylab(bquote(log(Rstar[sp1]/Rstar[sp2])))+theme(legend.position = "top")

jpeg("plots/coexistance_runner.jpeg",width = 9,height=6,unit='in',res=200)
ggpubr::ggarrange(coex,dom,common.legend = TRUE,labels=c("a)","b)"))
dev.off()

ggplot()+
  geom_point(data=fracking,aes(x=logsens1sens2,y=logR1R2),size=1,color="grey")+
  geom_point(data=fracking3,aes(x=logsens1sens2,y=logR1R2,color=g_time_dif),size=3)+
  facet_wrap(~ave_chill)+ggthemes::theme_base()+
  geom_vline(xintercept = 0)+geom_hline(yintercept=0)+scale_color_viridis_c(),common.legend=TRUE,label=c("SPE","SPE+FRACT"))
dev.off()

ggplot(fracking2,aes(sp2_meangmax))+geom_histogram()+facet_wrap(~ave_chill)

