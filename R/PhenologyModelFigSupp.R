library(RColorBrewer)
nonsta <- c(500,500,0)
nyrs <- sum(nonsta)
megaDflag <- 0
R0ns_flag <- 0.5

setwd("~/Documents/git/projects/temporalvar/R")
source("sourcefiles/getEnvt.R")

## Two panel figure (shifting R0 timing)
pdf("graphs/modelruns/manuscript/modelsupp.pdf", width=8, height=3.5)
par(mfrow=c(1,2))
#Plot shifting tauP PDF
x <- seq(0, 1, length = 5000)
plot(x,dbeta(x, pns[1], qns[1]),type="l", ylim=c(0,4.5), 
     col= brewer.pal(9,"YlGn")[9],
     ylab=expression(paste("PDF(",tau[p],")")),
     xlab=expression(paste("Timing of Resource Pulse, ",tau[p])))
for (i in seq(1,8,1)) {
  points(x,dbeta(x,pns[60*i],qns[60*i]),type="l",col=brewer.pal(9,"YlGn")[9-i])
}
#plot realization of shifting tauP
plot(seq(1,nyrs,1),tauP,type="n",xlab="year", ylab=expression(paste(tau[p])))
rect(500,0.09,1035,0.83,col='light gray', border=NA)
points(seq(1,nyrs,1),tauP,type="l",col="dark gray")
points(seq(1,nyrs,1),tauP,type="p",bty="o",pch=16,
       col=c(rep(brewer.pal(9,"YlGn")[9],500),rep(rev(brewer.pal(9,"YlGn")),length.out=500, each=56)))
dev.off()



## Four panel (shifting R0 timing and size)
pdf("graphs/modelruns/manuscript/modelsupp4panel.pdf", width=7.5, height=6)
par(mfrow=c(2,2))
#Plot shifting tauP PDF
x <- seq(0, 1, length = 5000)
plot(x,dbeta(x, pns[1], qns[1]),type="l", ylim=c(0,4.5), 
     col= brewer.pal(9,"YlGn")[9],
     ylab=expression(paste("PDF(",tau[p],")")),
     xlab=expression(paste("Timing of Resource Pulse, ",tau[p])))
for (i in seq(1,8,1)) {
  points(x,dbeta(x,pns[60*i],qns[60*i]),type="l",col=brewer.pal(9,"YlGn")[9-i])
}
#plot realization of shifting tauP
plot(seq(1,nyrs,1),tauP,type="n",xlab="year", ylab=expression(paste(tau[p])))
rect(500,0.09,1035,0.83,col='light gray', border=NA)
points(seq(1,nyrs,1),tauP,type="l",col="dark gray")
points(seq(1,nyrs,1),tauP,type="p",bty="o",pch=16,
       col=c(rep(brewer.pal(9,"YlGn")[9],500),rep(rev(brewer.pal(9,"YlGn")),length.out=500, each=56)))

#plot shifting R0 PDF
xR <- seq(.5,4,length=10000)
plot(xR,dlnorm(xR, mu_R0ns[1], sigma),type="l", ylim=c(0,1.5),
     col= brewer.pal(9,"Purples")[9],
     ylab=expression(paste("PDF(R0)")),
     xlab=expression(paste("Amount of Resource Pulse, ",R(0))))
for (i in seq(1,8,1)){
  points(xR,dlnorm(xR,mu_R0ns[60*i],sigma),type="l",col=brewer.pal(9,"Purples")[9-i])
}

#plot realizaton of shifting R0
plot(seq(1,nyrs,1),R0,type="n",xlab="year")
rect(500,0.75,1035,3.75,col='gainsboro', border=NA)
points(seq(1,nyrs,1),R0,type="l",col="dark gray")
points(seq(1,nyrs,1),R0,type="p",bty="o",pch=16,
       col=c(rep(brewer.pal(9,"Purples")[9],500),rep(rev(brewer.pal(9,"Purples")),length.out=500, each=56)))
dev.off()

