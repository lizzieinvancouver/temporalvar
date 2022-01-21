#quick plots of within season dynamics

par(mfrow=c(1,2))
plot(Bout[[y]]$R~Bout[[y]]$time, type="l",
     xlab="days",ylab=NA, main="Resource")
plot(Bout[[y]]$B1~Bout[[y]]$time, type="l", ylim=c(0,max(Bout[[y]]$B1,Bout[[y]]$B2)),
     xlab="days",ylab=NA,main="Sp1 & Sp2 Density")
lines(Bout[[y]]$B2~Bout[[y]]$time, type="l",col="blue")

###################################################
### Dan is bad at lists, make a data frame######
###############################################
runnum <- 1:10 # each run gets a number
Bout <- mapply(cbind, Bout, "RunID"=runnum, SIMPLIFY=F) #  assign that number
Bout.df<-do.call(rbind.data.frame, Bout) ### make the list a data frame

a<-ggplot(Bout.df,aes(time,R))+geom_smooth()+facet_wrap(~as.factor(RunID)) ##plot

bout2<-tidyr::gather(Bout.df,"species","biomass",3:4) #clean

b<-ggplot(bout2,aes(time,biomass,color=species))+geom_smooth()+facet_wrap(~as.factor(RunID)) #plot2
jpeg("plots/withinseas_firsttime.jpeg")
ggpubr::ggarrange(a,b)
dev.off()
