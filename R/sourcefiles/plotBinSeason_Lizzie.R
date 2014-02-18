dev.new(width=14, height=10)
par(mfrow=c(3,3))
par(mar=c(5, 4, 4, 8) + 0.1)
selectyrs <- seq(1, nyrs, by=floor(nyrs/8))
yrlength <- ndays/dt
for (yr in seq_along(selectyrs)){
    plot(Bout[[selectyrs[yr]]][,3]~Bout[[selectyrs[yr]]]$time, 
         ylim=c(0, max(Bout[[selectyrs[yr]]][3:(nsp+2)])), 
         xlab="step, step, step",  ylab="Abundance", type="n",
         main=paste("year: ", selectyrs[yr], sep=""))
    for (sp in c(1:nsp)){
        lines(Bout[[selectyrs[yr]]][,2+sp]~Bout[[selectyrs[yr]]]$time, 
              ylim=c(0, max(Bout[[selectyrs[yr]]][,nsp+2])),
              col=colerz[sp])
    }
    par(new=TRUE)
    plot(Bout[[selectyrs[yr]]]$R~Bout[[selectyrs[yr]]]$time, axes=FALSE, xlab="", ylab="",
        ylim=c(min(Bout[[selectyrs[yr]]]$R), max(Bout[[selectyrs[yr]]]$R)), type="l",
        col=rcol, lty=lresbyrs, lwd=lwd)
    raxis <- seq(0, max(Bout[[selectyrs[yr]]]$R), by=max(Bout[[selectyrs[yr]]]$R)/10)
    axis(4, at=raxis, labels=round(raxis, digits=2))
    mtext("Resource", side=4, line=3, cex=0.75)
}
