#####  getPlots
# When doing a few runs locally, plots basic info about each run to a one page pdf

##if this is being called from getOutput, rather than PlotCoexist, need to define IDx's
if (exists("runIDx")) {
  jobID<-c(jobIDx,taskIDx)
  j <- runIDx
}

pdfplot <- paste0(PlotOut_loc, "/Coexist_",jobID[1],"-",jobID[2],"_",j,".pdf")

ttl <-paste0(ifelse(megaDflag,paste0("Megadrought, "),paste0("")),"Sp = ", c(1,2),
                             ", tauI = ",  sprintf("%.2f",tauI), 
                             ", alpha = ", sprintf("%.2f",alpha),
                             ", s = ",     sprintf("%.2f",s),
                             ", Rstar = ", sprintf("%3.2e",Rstar), 
                             ", Job = ", jobID[1],"-",jobID[2],"_",j)
pdf(file=pdfplot, width=7.5,height=10,paper="letter",onefile=TRUE)

layout(matrix(c(5,5,1,3,2,4), 3, 2, byrow = TRUE))

col_1=rgb(0,0,1,.5)
col_2=rgb(1,0,0,.5)

##Print out overlapping histograms of tauIP for sp1 and sp 2 for stationary and nonstationary

h1t1=hist(tauP[seq(1,nonsta[1])]-tauIhat[seq(1,nonsta[1]),1],plot=FALSE)
h2t1=hist(tauP[seq(1,nonsta[1])]-tauIhat[seq(1,nonsta[1]),2],plot=FALSE)
h1t2=hist(tauP[seq(nonsta[1]+1,sum(nonsta))]-tauIhat[seq(nonsta[1]+1,sum(nonsta)),1],plot=FALSE)
h2t2=hist(tauP[seq(nonsta[1]+1,sum(nonsta))]-tauIhat[seq(nonsta[1]+1,sum(nonsta)),2],plot=FALSE)
ymax=max(h1t1$counts,h2t1$counts,h1t2$counts,h2t2$counts)
xmin=min(h1t1$breaks,h2t1$breaks,h1t2$breaks,h2t2$breaks)
xmax=max(h1t1$breaks,h2t1$breaks,h1t2$breaks,h2t2$breaks)

plot(h1t1$mids,h1t1$counts,type="s",col=col_1, ylim=c(0,ymax), xlim=c(xmin,xmax),
     xlab="TauIP",ylab=("Frequency"),lwd=3,main="Stationary")
lines(h2t1$mids,h2t1$counts,type="s",col=col_2,lwd=3)
legend("topleft",legend=c("Sp1","Sp2"), col=c(col_1,col_2),lty=1,lwd=2,bty="n")

plot(h1t2$mids,h1t2$counts,type="s",col=col_1,ylim=c(0,ymax), xlim=c(xmin,xmax),
     xlab="TauIP",ylab=("Frequency"),lwd=3,main="NonStationary")
lines(h2t2$mids,h2t2$counts,type="s",col=col_2,lwd=3)
legend("topleft",legend=c("Sp1","Sp2"), col=c(col_1,col_2),lty=1,lwd=2,bty="n")

##Print out overlapping histograms of germination for sp1 and sp 2 for stationary and nonstationary

g1t1=hist(g[seq(1,nonsta[1]),1],plot=FALSE)
g2t1=hist(g[seq(1,nonsta[1]),2],plot=FALSE)
g1t2=hist(g[seq(nonsta[1]+1,sum(nonsta)),1],plot=FALSE)
g2t2=hist(g[seq(nonsta[1]+1,sum(nonsta)),2],plot=FALSE)

ymax=max(g1t1$counts,g2t1$counts,g1t2$counts,g2t2$counts)
xmin=min(g1t1$breaks,g2t1$breaks,g1t2$breaks,g2t2$breaks)
xmax=max(g1t1$breaks,g2t1$breaks,g1t2$breaks,g2t2$breaks)

gmean100sta <- colMeans(g[seq(nonsta[1]-100,nonsta[1]),])
gmean100ns <-  colMeans(g[seq((length(g[,1]) - 100),length(g[,1])),])

plot(g1t1$mids,g1t1$counts,type="s",col=col_1, ylim=c(0,ymax), xlim=c(xmin,xmax),lwd=3,
     xlab="Germination",ylab=("Frequency"),
     main=paste0("Stationary: gmean[100] = ",paste0(sprintf("%3g",gmean100sta),collapse=", ")))
lines(g2t1$mids,g2t1$counts,type="s",col=col_2,lwd=3)
legend("top",legend=c("Sp1","Sp2"), col=c(col_1,col_2),lty=1,lwd=2,bty="n")

plot(g1t2$mids,g1t2$counts,type="s",col=col_1,ylim=c(0,ymax), xlim=c(xmin,xmax),
     xlab="Germination",ylab=("Frequency"),lwd=3,
     main = paste0("Nonstationary: gmean[100] = ",paste0(sprintf("%3g",gmean100ns),collapse=", ")))
lines(g2t2$mids,g2t2$counts,type="s",col=col_2,lwd=3)
legend("top",legend=c("Sp1","Sp2"), col=c(col_1,col_2),lty=1,lwd=2,bty="n")

#print out interannual variation for both species through time

plot(log(N[,1])~seq(1,sum(nonsta)),col=c(col_1),type="l",ylab="log(Seed Bank)", xlab="Time",lwd=3,
     main=ttl,ylim=c(0,max(log(N),na.rm=TRUE)))
lines(log(N[,2])~seq(1,sum(nonsta)),col=c(col_2),type="l",lwd=3)
legend("top",legend=c("Sp1","Sp2"), col=c(col_1,col_2),lty=1,lwd=2,bty="n")

dev.off()
