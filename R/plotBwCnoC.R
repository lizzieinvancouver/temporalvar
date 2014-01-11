#compare biomass within season with and without competition

dev.new(width=14, height=10)
par(mfrow=c(3,2))
for (i in seq(1,length(plotyrs),3)){
  q=plotyrs[i]
  #plot dynamics WITH comptition
  plot(Bout[[q]][,3]~Bout[[q]]$time, ylim=c(0, max(Bfin[plotyrs,])),
       xlab="step, step, step", ylab="Biomass", type="n",main=plotyrs[i])
  for (j in 1:nsp){
    lines(Bout[[q]][,2+j]~Bout[[q]]$time, col=colerz[j], lty=lspbyrs, lwd=lwd)
  }
  #overlay R on a second y axis
  par(new=TRUE)
  plot(Bout[[q]]$R~Bout[[q]]$time,type="l",col="black",lty=2,ylim=c(0,max(R0)),xaxt="n",yaxt="n",xlab="",ylab="")
  axis(4)
  #mtext("Resource",side=4,line=3,cex=0.7)
  
  #plot dynamics WIHTOUT competition
  plot(BnoCout[[q]][,3]~BnoCout[[q]]$time, ylim=c(0, max(BnoC[plotyrs,])),
       xlab="step, step, step", ylab="Biomass", type="n",main=plotyrs[i])
  for (j in 1:nsp){
    lines(BnoCout[[q]][,2+j]~BnoCout[[q]]$time, col=colerz[j], lty=lspbyrs, lwd=lwd)
  }
  #overlay R on a second y axis
  par(new=TRUE)
  plot(BnoCout[[q]]$R~BnoCout[[q]]$time,type="l",col="black",lty=2,ylim=c(0,max(R0)),xaxt="n",yaxt="n",xlab="",ylab="")
  axis(4)
  #mtext("Resource",side=4,line=3,cex=0.7)
  
}
