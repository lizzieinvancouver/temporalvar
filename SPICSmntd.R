#phylogenetic signal in SPCIS

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

#library("hkdatasets")
library(sf)
library(dplyr)
#library("mapview")
#library("tmap")
#library(maps)
#library(usmap)
library(phytools)


setwd("~/Documents/git/bioticHogs/")

d<-read.csv("Data/SPCIS_plant_taxa.csv") ### this is the public data for reproducibility

p<-read.csv("Data/SPCIS_plots.csv")

p<-p %>% #select only the most recent survey for resampled plots
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 

regionz<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")
regionz<-dplyr::select(regionz,Plot,NA_L1NAME,NA_L2NAME,NA_L3NAME,US_L4NAME)
regionz<-distinct(regionz)


d<-left_join(d,regionz)
d<-d %>% #select only the most recent survey for resampled plots
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 


phy<-read.tree("Analyses/Master/Input/FINAL_PHYLOGENY_SpCodes_Apr2023.tr")
d<-filter(d,SpCode %in% c(phy$tip.label)) 

##get only plots with 5% cover invasive
invaders<-d %>% dplyr::group_by(Plot,NativeStatus)%>% dplyr::summarize(RelCov=sum(PctCov_100,na.rm=TRUE))
invaders<-filter(invaders,NativeStatus=="I")
invaders<-filter(invaders,RelCov>5)

d.i<-filter(d,Plot %in% invaders$Plot)

p.i<-filter(p,Plot %in% invaders$Plot)
p.samp<-sample_n(p.i,100)

d.samp<-d

d.samp<-dplyr::select(d.samp,SpCode,NativeStatus,NA_L1NAME)
d.samp$nonnative<-ifelse(d.samp$NativeStatus=="I",1,0)
d.samp<-distinct(d.samp)
d.samp<-filter(d.samp,!is.na(NA_L1NAME))

phy.dist<-cophenetic(phy)


matricize<-function(x){
  temp<-dplyr::select(x,NA_L1NAME,SpCode,nonnative) ### select useful columns
  temp<-dplyr::filter(temp,!is.na(temp$SpCode))
  temp<- tidyr::spread(temp,SpCode,nonnative)
  temp[is.na(temp)] <- 0 
  temp<-as.data.frame(temp)
  rownames(temp)<-temp$NA_L1NAME # convert plot names to row names
  temp<-dplyr::select(temp,-NA_L1NAME)
  temp<-as.matrix.data.frame(temp)#
}


L1<-matricize(d.samp)


library(caper)
library(picante)
#goo.mpd<-mpd(goo, phy.dist)
#USA.ses.mpd<-ses.mpd(L1, phy.dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 9999)
USA.ses.mntd<-ses.mntd(L1, phy.dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 9999)

uno<-data.frame(NA_L1NAME=c(rownames(USA.ses.mpd)),Estimate=c(USA.ses.mpd$mpd.rand.mean))

uno$lower<-c(uno$Estimate-USA.ses.mpd$mpd.rand.sd)
uno$upper<-c(uno$Estimate+USA.ses.mpd$mpd.rand.sd)
uno$lower2<-c(uno$Estimate-2*USA.ses.mpd$mpd.rand.sd)
uno$upper2<-c(uno$Estimate+2*USA.ses.mpd$mpd.rand.sd)

uno$metric<-c("MPD")
uno$cat<-"randomly generated"

uno2<-data.frame(NA_L1NAME=c(rownames(USA.ses.mpd)),Estimate=USA.ses.mpd$mpd.obs)
uno2$lower<-uno2$Estimate
uno2$upper<-uno2$Estimate
uno2$lower2<-uno2$Estimate
uno2$upper2<-uno2$Estimate
uno2$metric<-c("MPD")
uno2$cat<-"observed"
uno<-rbind(uno,uno2)
uno<-filter(uno,NA_L1NAME!="TROPICAL WET FORESTS")


mpd.plot<-ggplot(uno,aes(NA_L1NAME,Estimate))+
  geom_errorbar(aes(ymin=lower,ymax=upper,color=cat),width=0.0,size=.5)+
  geom_errorbar(aes(ymin=lower2,ymax=upper2,color=cat),width=0,size=.5,linetype="dotdash")+
  geom_point(aes(shape=cat,color=cat),size=5)+
  scale_x_discrete(name="")+ggthemes::theme_few(base_size = 9)+
  scale_shape_manual(name="",values=c(16,1))+scale_color_manual(name="",values=c("black","firebrick4"))+
  theme(legend.position = "top")+ylab("Mean phylogenetic distance")+
  theme(axis.text.x = element_text(angle = 355, vjust = .5, hjust=.5))
  #theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())


uno<-data.frame(NA_L1NAME=c(rownames(USA.ses.mntd)),Estimate=c(USA.ses.mntd$mntd.rand.mean))

uno$lower<-c(uno$Estimate-USA.ses.mntd$mntd.rand.sd)
uno$upper<-c(uno$Estimate+USA.ses.mntd$mntd.rand.sd)
uno$lower2<-c(uno$Estimate-2*USA.ses.mntd$mntd.rand.sd)
uno$upper2<-c(uno$Estimate+2*USA.ses.mntd$mntd.rand.sd)

uno$metric<-c("MNTD")
uno$cat<-"randomly generated"

uno2<-data.frame(NA_L1NAME=c(rownames(USA.ses.mntd)),Estimate=USA.ses.mntd$mntd.obs)
uno2$lower<-uno2$Estimate
uno2$upper<-uno2$Estimate
uno2$lower2<-uno2$Estimate
uno2$upper2<-uno2$Estimate
uno2$metric<-c("MNTD")
uno2$cat<-"observed"
uno<-rbind(uno,uno2)
uno<-filter(uno,NA_L1NAME!="TROPICAL WET FORESTS")


mntd.plot<-ggplot(uno,aes(reorder(NA_L1NAME,Estimate),Estimate))+
  geom_errorbar(aes(ymin=lower,ymax=upper,color=cat),width=0.0,size=.5)+
  geom_errorbar(aes(ymin=lower2,ymax=upper2,color=cat),width=0,size=.5,linetype="dotdash")+
  geom_point(aes(shape=cat,color=cat),size=3)+
  scale_x_discrete(name="")+ggthemes::theme_few(base_size = 9)+
  scale_shape_manual(name="",values=c(16,1))+scale_color_manual(name="",values=c("black","firebrick4"))+
  theme(legend.position = "top")+ylab("Mean nearest taxon distance")+
  theme(axis.text.x = element_text(angle = 355, vjust = .5, hjust=.5))

mntd.small<-filter(uno,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

mntd.plot.small<-ggplot(mntd.small,aes(reorder(NA_L1NAME,Estimate),Estimate))+
  geom_errorbar(aes(ymin=lower,ymax=upper,color=cat),width=0.0,size=.5)+
  geom_errorbar(aes(ymin=lower2,ymax=upper2,color=cat),width=0,size=.5,linetype="dotdash")+
  geom_point(aes(shape=cat,color=cat),size=4)+
  scale_x_discrete(name="")+ggthemes::theme_few(base_size = 9)+
  scale_shape_manual(name="",values=c(16,1))+scale_color_manual(name="",values=c("black","firebrick4"))+
  theme(legend.position = "top")+ylab("Mean nearest taxon distance")+coord_cartesian(ylim=c(10,85))

jpeg("Analyses/Master/Plots/spcis_phylosig.jpeg",width = 11,height=3,unit='in',res=300)
mntd.plot.small
dev.off()

jpeg("Analyses/Master/Plots/spcis_phylosigall.jpeg",width = 11,height=4,unit='in',res=300)
mntd.plot
dev.off()


ETF<-L1[1,]
NAD<-L1[5,]
L4<-rbind(ETF,NAD)

 par(mfrow = c(1, 2))
 for (i in row.names(L4)) {
  plot(phy, show.tip.label = FALSE, main = i)
  tiplabels(tip = which(phy$tip.label %in% names(which(L4 [i, ] > 0))),
               pch=15, cex=.5,col="red")
   }
