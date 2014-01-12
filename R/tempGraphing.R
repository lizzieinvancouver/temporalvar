plotyrs<- floor(seq(10,nyrs-1,length.out=9))

par(mfrow=c(3,3))
for (i in 1:length(plotyrs)){
  q=plotyrs[i]
  plot(B[q,1,]~c(1:tsteps), ylim=c(0,max(Bfin[plotyrs,])),
       xlab="step, step, step", ylab="Biomass", type="n",main=plotyrs[i])
  for (j in 1:nsp){
    lines(B[q,j,]~c(1:tsteps), col=colerz[j], lty=lspbyrs, lwd=lwd)
  }
  #overlay R on a second y axis
  par(new=TRUE)
  plot(R[q,]~c(1:tsteps),type="l",col="black",lty=2,ylim=c(0,max(R0[plotyrs])),xaxt="n",yaxt="n",xlab="",ylab="")
  axis(4)
  #mtext("Resource",side=4,line=3,cex=0.7)
}


par(mfrow=c(3,3))
for (i in 1:length(plotyrs)){
  q=plotyrs[i]
  plot(BnoC[q,1,]~c(1:tsteps), ylim=c(0, max(BnoCfin[plotyrs,])),
       xlab="step, step, step", ylab="Biomass", type="n",main=plotyrs[i])
  for (j in 1:nsp){
    lines(BnoC[q,j,]~c(1:tsteps), col=colerz[j], lty=lspbyrs, lwd=lwd)
  }
  #overlay R on a second y axis
  par(new=TRUE)
  plot(RnoC[q,]~c(1:tsteps),type="l",col="black",lty=2,ylim=c(0,max(R0[plotyrs])),xaxt="n",yaxt="n",xlab="",ylab="")
  axis(4)
  #mtext("Resource",side=4,line=3,cex=0.7)
}
