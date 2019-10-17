## Started 12 August 2019 ##
## By Lizzie ##

## Working on idealized plotting ## 

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/projects/temporalvar/R")

## SpA is non-tracking; spB is tracking
taui.A <- 0.58
taui.B <- 0.7
alpha.A <- 0
alpha.B <- 0.5
tauP <- c(0.6,0.3)
tauihat.A <- alpha.A*tauP+(1-alpha.A)*taui.A
tauihat.B <- alpha.B*tauP+(1-alpha.B)*taui.B
spcol <- list(spA =c("darkorange","darkorange"),spB=c("dodgerblue4","dodgerblue"))

## tauP in STA and NS; and two draws
x <- seq(0, 1, length = 10000)
tauP.STA <- dbeta(x, 10, 10)
tauP.NST <- dbeta(x, 5, 15)
tauPcol <- c("darkgreen","chartreuse3")

##germination for sp A and B
gmax <- 0.5
h <- 100
#germination given tauP
g.A <- gmax*exp(-h*(x -(alpha.A*x+(1-alpha.A)*taui.A))^2)
g.B.tr <- gmax*exp(-h*(x -(alpha.B*x+(1-alpha.B)*taui.B))^2)
g.B.ntr <- gmax*exp(-h*(x -(0*x+(1-0)*taui.B))^2)

#germination times P(tauP)

g.A.STA <- g.A*tauP.STA
g.A.NST <- g.A*tauP.NST
g.B.tr.STA <- g.B.tr*tauP.STA
g.B.tr.NST <- g.B.tr*tauP.NST
g.B.ntr.STA <- g.B.ntr*tauP.STA
g.B.ntr.NST <- g.B.ntr*tauP.NST

#Set up for three panels (tauP dist, germ sp A, germ sp B)
pdf("R/graphs/conceptual/TauP_Germination.pdf", width=6, height=9)
def.par <- par(no.readonly = TRUE) # save default, for resetting...
#nf<-layout(matrix(c(1,2,3),3,1,byrow=TRUE),widths=2,heights=1)
#layout.show(nf)
par(mfrow=c(3,1), oma = c(2, 0, 0, 0))

## PLOT ONE:  Plot the ideal tauP distributions ...

lwdhere <- 2.5
par(mar=c(0,4,0,0),mgp=c(1.5,1,0))
plot(x,tauP.NST, type="l", axes=FALSE, bty="o",cex.lab=1.5,
     ylab = expression(paste("PDF(",tau[p],")")),xlab="", 
     ylim=c(-2, 5),col=tauPcol[2],xaxp=c(0,1,1), lwd=lwdhere)
lines(x,tauP.STA, ylab="", xlab="", yaxt="n", col=tauPcol[1], lwd=lwdhere)
#tauI arrows
arrows(taui.A,-.75,taui.A,-0.05,length=0.075, col = spcol[[1]][1])
text(taui.A,-.75,pos=1,expression(paste(tau[i[A]])),cex=1.75,font=3, col=spcol[[1]][1])
arrows(taui.B,-.75,taui.B,-0.05,length=0.075, col = spcol[[2]][2])
text(taui.B,-.75,pos=1,expression(paste(tau[i[B]])),cex=1.75,font=3, col=spcol[[2]][2])
axis(1,at=c(0,1),pos=0)
#tauP arrows
arrows(tauP[1],1.25,tauP[1],0.05,length=0.1,col=tauPcol[1])
text(tauP[1],1.25,pos=3,expression(paste(tau[p[t]])),cex=1.75,font=3, col=tauPcol[1])
arrows(tauP[2],1.25,tauP[2],0.05,length=0.1,col=tauPcol[2])
text(tauP[2],1.25,pos=3,expression(paste(tau[p[t]])),cex=1.75,font=3, col=tauPcol[2])

#PLOT TWO: plot effective germination sp A
par(mar=c(0,4,2,0),mgp=c(1.5,1,0))
plot(x,g.A.STA, type="l", col=spcol[[1]][1], 
     ylim=c(-0.6,1.75), axes=FALSE, bty="n", cex.lab=1.5,
     ylab = "Germination Rate",
     xlab="", lwd=lwdhere)
lines(x,g.B.tr.STA,  col=spcol[[2]][1],lty=1, lwd=lwdhere)
lines(x,g.B.ntr.STA,  col=spcol[[2]][2],lty=1, lwd=lwdhere)
axis(1,at=c(0,1),pos=-0.01)
axis(2,at=c(0,0.75,1.5),pos=0)

#plot tauIhat arrows
arrows(tauihat.A[1],-.25,tauihat.A[1],-0.01,length=0.075, col = spcol[[1]][1])
text(tauihat.A[1],-.25,pos=1,expression(paste(widehat(tau)[i[A]])),cex=1.75,font=3, col = spcol[[1]][1])

arrows(tauihat.B[1],-.25,tauihat.B[1],-0.01,length=0.075, col = spcol[[2]][1])
text(tauihat.B[1],-.25,pos=1,expression(paste(widehat(tau)[i[B]])),cex=1.75,font=3, col = spcol[[2]][1])

arrows(taui.B,-.2,tauihat.B[1]+.005,-0.2,length=0.075, col = spcol[[2]][2])
text(taui.B[1],-.3,pos=4,expression(paste(tau[i[B]])),cex=1.75,font=3, col = spcol[[2]][2])

#plot tauP arrow
arrows(tauP[1],0.5,tauP[1],0,length=0.1, col=tauPcol[1])
text(tauP[1],0.5,pos=3,expression(paste(tau[p[t]])),cex=1.75,font=3, col = tauPcol[1])


#PLOT THREE: Germination at end of nonstationary
par(mar=c(0,4,0,0),mgp=c(1.5,1,0))
plot(x,g.A.NST, type="l", col=spcol[[1]][2],  
     ylim=c(-0.05,.2), axes=FALSE, bty="n", cex.lab=1.5,
     ylab = "Germination Rate", lwd=lwdhere)
lines(x,g.B.ntr.NST,  col=spcol[[2]][2],lty=1, lwd=lwdhere)
lines(x,g.B.tr.NST,  col=spcol[[2]][1],lty=1, lwd=lwdhere)
axis(1,at=c(0,1),pos=-0.001)
axis(2,at=c(0,0.075,.15),pos=0)
#tauP arrows
arrows(tauP[2],.05,tauP[2],0,length=0.1, col=tauPcol[2])
text(tauP[2],.05,pos=3,expression(paste(tau[p[
  t]])),cex=1.75,font=3, col = tauPcol[2])
#tauIhat arrow
arrows(tauihat.A[2],-.025,tauihat.A[2],-0.001,length=0.075, col = spcol[[1]][2])
text(tauihat.A[2],-.025,pos=1,expression(paste(widehat(tau)[i[A]])),cex=1.75,font=3, col = spcol[[1]][2])

arrows(tauihat.B[2],-.025,tauihat.B[2],-0.001,length=0.075, col = spcol[[2]][1])
text(tauihat.B[2],-.025,pos=1,expression(paste(widehat(tau)[i[B]])),cex=1.75,font=3, col = spcol[[2]][1])

arrows(taui.B,-.015,tauihat.B[2]+0.005,-0.015,length=0.075, col = spcol[[2]][2])
text(taui.B,-.025,pos=4,expression(paste(tau[i[B]])),cex=1.75,font=3, col = spcol[[2]][2])

#overall x-axis title at bottom
mtext(expression(paste("Timing of Resource Pulse, ",tau[p])),outer=TRUE,side=1,line=1,cex=1)

dev.off()

# ## Plotting two instances of tau P
# for (y in c(1,2)) {
#   par(mar=c(0,0,0,0))
#   plot(x,rep(0,length(x)),type="n",ylim = c(-3.5,3),axes=FALSE, ylab="",xlab="")
#   axis(1,at=c(0,1),pos=0)
#   arrows(tauP[y],1.5,tauP[y],0,length=0.1, col=tauPcol[y])
#   text(tauP[y],1.5,pos=3,expression(paste(tau[p[t]])),cex=2,font=3, col = tauPcol[y])
#   arrows(tauihat.A[y],-1.5,tauihat.A[y],-0.05,length=0.1, col = spcol[[1]][y])
#   text(tauihat.A[y],-1.5,pos=1,expression(paste(widehat(tau)[i[A]])),cex=2,font=3, col = spcol[[1]][y])
#   arrows(tauihat.B[y],-1.5,tauihat.B[y],-0.05,length=0.1, col = spcol[[2]][y])
#   text(tauihat.B[y],-1.5,pos=1,expression(paste(widehat(tau)[i[B]])),cex=2,font=3, col = spcol[[2]][y])
# }
     
     
