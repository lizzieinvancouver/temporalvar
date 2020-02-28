## Started 21 Feb 2020 ##
## By Megan ##

## Working on plot to compare the difference between priority effects and interannual variation in germination##
##Six Panel Plot (3 rows, 2 cols):
# Left column:  Priority Effects
#Right column:  Interannual Variation
# Top row:  %new germinants v Day of Year
#Second row:  % of population tat has germinated vs Year
#Third row:  Priroity effects within year (left panel) and interannual veriation in germiantion (right panel) 

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
#setwd("~/Documents/git/projects/temporalvar/R")
setwd(paste0(getwd(),"/R"))


#Set up calculations for Panels A = Germination v DOY  & B Germination v Year for germination model
tauP.alpha = 5
tauP.beta = 15
yrs <- seq(1,100,1)

## SpA is non-tracking; spB is tracking
taui.A <- 0.58
taui.B <- 0.7
tauP.t <- rbeta(length(yrs), tauP.alpha, tauP.beta)
spcol <- list(spA =c("darkorange"),spB=c("dodgerblue"))

## tauP in STA and NS; and two draws
doy.max <- 120
doy <- seq(0, doy.max, 1)
x <- doy/doy.max
tauP <- dbeta(x, tauP.alpha,tauP.beta)
tauP.col <- c("chartreuse3")

##germination for sp A and B
gmax <- 0.5
h <- 100
#germination conditional on tauP=x
g.A.cond <- gmax*exp(-h*(x - taui.A)^2)
g.B.cond <- gmax*exp(-h*(x - taui.B)^2)

#marginal germination = (germination|tauP=X) times P(tauP)
g.A <- g.A.cond*tauP
g.B <- g.B.cond*tauP


#time series of germination within a year (t=1)
g.A.doy <- g.A*(x==tauP.t[1])                          
g.B.doy <- g.B*(x==tauP.t[1])
g.A.cum <- cumsum(g.A.doy)
g.B.cum <- cumsum(g.B.doy)


#time series of germination between years (t=1:100)
g.A.yr <- g.A[x==tauP.t]

###PICK UP HERE:  NEED TO FIND AWAY TO LOOKUP VALUES OF g.A such tath the series g.A.yr[t=1] is 
#   g.A[tauP==tauP.t]
#for tauP.t[1], find index of tauP s.t. tauP[index-1]<tauP.t[1] AND tauP[index]>= tauP.t[1]
#g.A[index] is the germibation in t==1

## PLOT ONE:  Plot the ideal tauP distributions ...
function f1(key,arr,index){
  #return index of array arr that is closest in value to key
  for (i in seq(1,length(a)))
}


#Set up for six panels 
pdf("graphs/conceptual/PriorityEff_BetHedge.pdf", width=6, height=9)
def.par <- par(no.readonly = TRUE) # save default, for resetting...
#nf<-layout(matrix(c(1,2,3),3,1,byrow=TRUE),widths=2,heights=1)
#layout.show(nf)
par(mfrow=c(3,1), oma = c(2, 0, 0, 0))



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
     
     
