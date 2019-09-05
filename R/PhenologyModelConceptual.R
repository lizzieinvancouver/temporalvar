## Started 12 August 2019 ##
## By Lizzie ##

## Working on idealized plotting ## 

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## SpA is non-tracking; spB is tracking 
taui.A <- 0.55
taui.B <- 0.75
alpha.A <- 0
alpha.B <- 0.5
tauP <- c(0.45,0.18)
tauihat.A <- alpha.A*tauP+(1-alpha.A)*taui.A
tauihat.B <- alpha.B*tauP+(1-alpha.B)*taui.B
spcol <- c("purple","green")

## tauP in STA and NS; and two draws
x <- seq(0, 1, length = 10000)
plot.start.tauP <- dbeta(x, 10, 10)
plot.end.tauP <- dbeta(x, 5, 15)
tauPcol <- c("dark blue","light blue")

##germination for sp A and B
gmax <- 0.5
h <- 100
g.A <- gmax*exp(-h*(matrix(rep(x,2),nrow=length(x),ncol=2,byrow=FALSE)-matrix(rep(tauihat.A,length(x)),nrow=length(x), ncol=2,byrow=TRUE))^2)
g.B <- gmax*exp(-h*(matrix(rep(x,2),nrow=length(x),ncol=2,byrow=FALSE)-matrix(rep(tauihat.B,length(x)),nrow=length(x), ncol=2,byrow=TRUE))^2)

#Set up for four panels
def.par <- par(no.readonly = TRUE) # save default, for resetting...
nf<-layout(matrix(c(1,1,2,3,4,4),6,1,byrow=TRUE),widths=2,heights=1)
layout.show(nf)

## Plot the ideal tauP distributions ...
par(mar=c(2,0,1,0))
plot(x,plot.end.tauP, type="l", axes=FALSE, bty="o",
     ylab = expression(paste("P(",tau[p],")")),xlab="", 
     ylim=c(-2, 5),col=tauPcol[2],xaxp=c(0,1,1))
lines(x,plot.start.tauP, ylab="", xlab="", yaxt="n", col=tauPcol[1])
#legend("topright",legend=c("Stationary Period","End of NonStationary Period"), col=c("black","cyan"), cex=.7, lty=1, bty="n")
arrows(.5,4.75,.225,4.75,length=0.1,col="black")
arrows(taui.A,-1,taui.A,-0.05,length=0.1, col = spcol[1])
text(taui.A,-1,pos=1,expression(paste(tau[i[A]])),cex=2,font=3, col=spcol[1])
arrows(taui.B,-1,taui.B,-0.05,length=0.1, col = spcol[2])
text(taui.B,-1,pos=1,expression(paste(tau[i[B]])),cex=2,font=3, col = spcol[2])
axis(1,at=c(0,1),pos=0)

## Plotting two instances of tau P
ycol=c("dark","light")
for (y in c(1,2)) {
  par(mar=c(0,0,0,0))
  plot(x,rep(0,length(x)),type="n",ylim = c(-3.5,3),axes=FALSE, ylab="",xlab="")
  axis(1,at=c(0,1),pos=0)
  arrows(tauP[y],1.5,tauP[y],0,length=0.1, col="black")
  text(tauP[y],1.5,pos=3,expression(paste(tau[p[t]])),cex=2,font=3, col = tauPcol[y])
  arrows(tauihat.A[y],-1.5,tauihat.A[y],-0.05,length=0.1, col = spcol[1])
  text(tauihat.A[y],-1.5,pos=1,expression(paste(hat(tau[i[A]]))),cex=2,font=3, col = spcol[1])
  arrows(tauihat.B[y],-1.5,tauihat.B[y],-0.05,length=0.1, col = spcol[2])
  text(tauihat.B[y],-1.5,pos=1,expression(paste(hat(tau[i[B]]))),cex=2,font=3, col = paste(ycol[y],spcol[2]))
}

#plot germination
par(mar=c(0,0,1,0))
plot(x,g.A[,1], type="l", axes=FALSE,bty="n",
     ylim=c(-2, .5),col=spcol[1],xaxp=c(0,1,1))
for (y in c(1,2)) {
  lines(x,g.B[,y], ylab="", xlab="", yaxt="n", col=paste(ycol[y],spcol[2]))
}
axis(1,at=c(0,1),pos=0)
text(.5,-1,pos=3,expression(paste("Timing of Resource Pulse",tau[p])),cex=1.5)
text(0,0.5,pos=3,"germination",cex=1.5,srt=90)
     
     