dev.new(width=14, height=10)
par(mfrow=c(3,3))
par(mar=c(5, 4, 4, 8) + 0.1)
yrlength <- ndays/dt
for (yr in seq_along(plotyrs)){
    plot(Bout[[plotyrs[yr]]][,3]~Bout[[plotyrs[yr]]]$time, 
         ylim=c(0, max(Bout[[plotyrs[yr]]][3:(nsp+2)])), 
         xlab="step, step, step",  ylab="Abundance", type="n",
         main=paste("year: ", plotyrs[yr], sep=""))
    for (sp in c(1:nsp)){
        lines(Bout[[plotyrs[yr]]][,2+sp]~Bout[[plotyrs[yr]]]$time, 
              ylim=c(0, max(Bout[[plotyrs[yr]]][,nsp+2])),
              col=colerz[sp])
    }
    par(new=TRUE)
    plot(Bout[[plotyrs[yr]]]$R~Bout[[plotyrs[yr]]]$time, axes=FALSE, xlab="", ylab="",
        ylim=c(min(Bout[[plotyrs[yr]]]$R), max(Bout[[plotyrs[yr]]]$R)), type="l",
        col=rcol, lty=lresbyrs, lwd=lwd)
    raxis <- seq(0, max(Bout[[plotyrs[yr]]]$R), by=max(Bout[[plotyrs[yr]]]$R)/10)
    axis(4, at=raxis, labels=round(raxis, digits=2))
    mtext("Resource", side=4, line=3, cex=0.75)
}

## Did program ever jump out of loop?
for (i in c(1:100)){
   print(length(Bout[[i]]$time))
}
