#Between Years Plots
iyr = 20  #initial year in plot, to remove transient
# between years plot
par(mfrow=c(1,1))
dev.new(width=14, height=10)
par(mar=c(5,4,4,5)+.1)
plot(Bfin[iyr:nyrs,1]~c(iyr:nyrs), ylim=c(0, max(Bfin[iyr:nyrs])),
     xlab="year", ylab="Abundance", type="n")
for (i in 1:nsp) {
  lines(Bfin[iyr:nyrs,i]~c(iyr:nyrs), col=colerz[i], lty=lspbyrs, lwd=lwd)
}
#overlay germination fraction
par(new=TRUE)
plot(g[iyr:nyrs,i]~c(iyr:nyrs),type="n",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("germination fraction (+)",side=4,line=3)
for (i in 1:nsp) {
  points(g[iyr:nyrs,i]~c(iyr:nyrs),col=colerz[i],pch="+")
}
#OVerlays Resource level, R0, each year
# par(new=TRUE)
# plot(R0~c(1:nyrs),type="l",col="black",lty=2,xaxt="n",yaxt="n",xlab="",ylab="")
# axis(4)
# mtext("R0",side=4,line=3)