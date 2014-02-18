# between years plot
dev.new(width=14, height=10)
plot(BnoCfin[,1]~c(1:nyrs), ylim=c(min(BnoCfin), max(BnoCfin)),
     xlab="year", ylab="Abundance", type="n")
for (i in 1:nsp) {
  plot(BnoCfin[,i]~c(1:nyrs), col=colerz[i], lty=lspbyrs, lwd=lwd,type="l")
  #lines(BnoCfin[,i]~c(1:nyrs), col=colerz[i], lty=lspbyrs, lwd=lwd)
}

dev.new(width=14, height=10)
par(mfrow=c(3,3))
par(mar=c(5, 4, 4, 8) + 0.1)
selectyrs <- seq(1, nyrs, by=floor(nyrs/8))
yrlength <- ndays/dt
for (yr in seq_along(selectyrs)){
  allsp <- BnoC[selectyrs[yr], ,]
  plot(BnoC[selectyrs[yr], 1,]~c(1:yrlength), ylim=c(min(allsp), max(BnoC)), 
       xlab="step, step, step",  ylab="Abundance", type="n",
       main=paste("year: ", selectyrs[yr], sep=""))
  for (sp in c(1:nsp)){
    lines(BnoC[selectyrs[yr], sp ,]~c(1:yrlength), ylim=c(min(allsp), max(BnoC)),
          col=colerz[sp])
  }
  
}
#   par(new=TRUE)
#   plot(RnoC[selectyrs[yr],1,]~c(1:yrlength), axes=FALSE, xlab="", ylab="",
#        ylim=c(min(R[selectyrs[yr],]), max(R[selectyrs[yr],])), type="n",
#        col=rcol, lty=lresbyrs, lwd=lwd)
#   for (sp in c(1:nsp)){
#     lines(RnoC[selectyrs[yr], sp ,]~c(1:yrlength), ylim=c(0, max(RnoC)),
#           col=colerz[sp])
#   }
#   #raxis <- seq(0, max(R[selectyrs[yr],]), by=max(R[selectyrs[yr],])/10)
#   #axis(4, at=raxis, labels=round(raxis, digits=2))
#   #mtext("Resource", side=4, line=3, cex=0.75)
# }



# plotyrs<- floor(seq(10,nyrs-1,length.out=9))
# 
# par(mfrow=c(3,3))
# for (i in 1:length(plotyrs)){
#   q=plotyrs[i]
#   plot(B[q,1,]~c(1:tsteps), ylim=c(0,max(Bfin[plotyrs,])),
#        xlab="step, step, step", ylab="Biomass", type="n",main=plotyrs[i])
#   for (j in 1:nsp){
#     lines(B[q,j,]~c(1:tsteps), col=colerz[j], lty=lspbyrs, lwd=lwd)
#   }
#   #overlay R on a second y axis
#   par(new=TRUE)
#   plot(R[q,]~c(1:tsteps),type="l",col="black",lty=2,ylim=c(0,max(R0[plotyrs])),xaxt="n",yaxt="n",xlab="",ylab="")
#   axis(4)
#   #mtext("Resource",side=4,line=3,cex=0.7)
# }
# 
# 
# par(mfrow=c(3,3))
# for (i in 1:length(plotyrs)){
#   q=plotyrs[i]
#   plot(BnoC[q,1,]~c(1:tsteps), ylim=c(0, max(BnoCfin[plotyrs,])),
#        xlab="step, step, step", ylab="Biomass", type="n",main=plotyrs[i])
#   for (j in 1:nsp){
#     lines(BnoC[q,j,]~c(1:tsteps), col=colerz[j], lty=lspbyrs, lwd=lwd)
#   }
#   #overlay R on a second y axis
#   par(new=TRUE)
#   plot(RnoC[q,]~c(1:tsteps),type="l",col="black",lty=2,ylim=c(0,max(R0[plotyrs])),xaxt="n",yaxt="n",xlab="",ylab="")
#   axis(4)
#   #mtext("Resource",side=4,line=3,cex=0.7)
# }
