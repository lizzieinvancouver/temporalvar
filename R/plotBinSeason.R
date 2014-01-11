#plot B within Season

dev.new(width=14, height=10)
par(mfrow=c(3,3))
for (i in 1:length(plotyrs)){
  q=plotyrs[i]
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
}
