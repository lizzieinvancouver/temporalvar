library(ggplot2)
runnum <- 1:50 # each run gets a number
Bout <- mapply(cbind, Bout, "RunID"=runnum, SIMPLIFY=F) #  assign that number
Bout.df<-do.call(rbind.data.frame, Bout) ### make the list a data frame

a<-ggplot(Bout.df,aes(time,R))+geom_smooth()+facet_wrap(~as.factor(RunID)) ##plot

bout2<-tidyr::gather(Bout.df,"species","biomass",3:4) #clean

b<-ggplot(bout2,aes(time,biomass,color=species))+geom_smooth()+facet_wrap(~as.factor(RunID)) #plot2
jpeg("plots/withinseas_upchillmeansd.jpeg")
ggpubr::ggarrange(a,b)
dev.off()
