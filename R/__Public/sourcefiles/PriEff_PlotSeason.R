#quick plots of within season dynamics

#Use plot function for deSolve; note that we 'harvest' when R=0

Bout.trim <- Bout[[y]]

for (i in seq(max(ind.Rstar),days+1,dt)) Bout.trim[[i]] <- NA

if (y%%10==3)plot(Bout.trim)

 par(mfrow=c(1,2))
 Bout.df<-as.data.frame(Bout[[y]])
 
 plot(Bout.df$R~Bout.df$time, type="l",
      xlab="days",ylab=NA, main="Resource")
 plot(Bout.df$B1~Bout.df$time, type="l", ylim=c(0,max(Bout.df$B1,Bout.df$B2)),
      xlab="days",ylab=NA,main="Sp1 & Sp2 Density")
 lines(Bout.df$B2~Bout.df$time, type="l",col="blue")

###################################################
### Dan is bad at lists, make a data frame######
###############################################

 #a<-ggplot(Bout.df,aes(time,R))+geom_smooth()+facet_wrap(~as.factor(RunID)) ##plot

#bout2<-tidyr::gather(Bout.df,"species","biomass",3:4) #clean

#b<-ggplot(bout2,aes(time,biomass,color=species))+geom_smooth()+facet_wrap(~as.factor(RunID)) #plot2
#jpeg("plots/withinseas_firsttime.jpeg")
#ggpubr::ggarrange(a,b)
#dev.off()
