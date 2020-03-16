## Started 21 Feb 2020 ##
## By Megan ##

## Working on plot to compare the difference between priority effects and interannual variation in germination##
##Four Panel Plot (2 rows, 2 cols):
# Left column:  Interannual Variation
# Right column:  Priority Effects
# Top row:   germination v Day of Year
#Second row:  germination by year

##QUESTION:  do we want a synthesis plot, showing interannual variation in germination +
###           priority effects within year


## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("C:/Users/Megan/Documents/GitHub/temporalvar/R")

#Set up calculations for Panels A = Germination v DOY  & B Germination v Year for germination model
tauP.alpha = 5
tauP.beta = 15
yrs <- seq(1,60,1)
tauP.t <- rbeta(length(yrs), tauP.alpha, tauP.beta)

## tauP in STA and NS; and two draws
doy.max <- 100
step=.25
doy <- seq(0, doy.max, step)
x <- doy/doy.max
tauP <- dbeta(x, tauP.alpha,tauP.beta)
tauPcol <- "chartreuse3"

##germination for sp A and B
taui.A <- 0.22
taui.B <- 0.35
spcol <- list(spA =c("darkorange"),spB=c("dodgerblue"))
gmax <- 0.5
h <- 100

#germination distribution --> germination} tauP=x 
g.A.cond <- gmax*exp(-h*(x - taui.A)^2)
g.B.cond <- gmax*exp(-h*(x - taui.B)^2)

#priority effects -> normalize annual germination to 1
g.A.priority <-g.A.cond
g.B.priority <- g.B.cond
g.A.priority.cum <- cumsum(g.A.priority)/sum(g.A.priority)
g.B.priority.cum <- cumsum(g.B.priority)/sum(g.B.priority)

#marginal germination = (germination|tauP=X) * P(tauP)
#Expected germination distribution
g.A <- g.A.cond*tauP
g.B <- g.B.cond*tauP

#time series of germination within a year (t=1)
#germination v time -->  all germination happens on one day, but the amount depends on tauP - tauI
#need a sequence of zeros with germination on, say, day 20

t= 10 #germation year ot plot
d=20/step #day germination happens
g.A.doy <- c(rep(0,d-1), gmax*exp(-h*(tauP.t[t] - taui.A)^2),rep(0,length(doy)-d))                         
g.B.doy <- c(rep(0,d-1), gmax*exp(-h*(tauP.t[t] - taui.B)^2),rep(0,length(doy)-d)) 
g.A.cum <- cumsum(g.A.doy)
g.B.cum <- cumsum(g.B.doy)

#time series of germination between years (t=1:100)
lookup.germ <- function(tauPx,x,g){
  #for an array of values of x in the tauP distribution, look up the germination rate and return array
  #x is index for g (i.e., germination at x)
  #tauPx and x are defined ont he same range (in this case, 0,1) but may be different length arrays
  ret = NULL
  for (i in seq(1,length(tauPx))){
    ret <- c(ret,g[length(x[x<tauPx[i]])])
  }
  return(ret)
}
g.A.yr <- lookup.germ(tauP.t,x,g.A)
g.B.yr <- lookup.germ(tauP.t,x,g.B)

#Set up for four panels 
pdf("graphs/conceptual/PriorityEff_BetHedge.pdf", width=8, height=6)
def.par <- par(no.readonly = TRUE) # save default, for resetting...
nf<-layout(matrix(c(1,2,3,4),2,2,byrow=TRUE))
#layout.show(nf)
par(mfrow=c(2,2), oma = c(2, 0, 0, 0))

#Top Right Panel: Germ Frac & Cumulative v DOY 
par(mar=c(2,2,2,2))
yup <- max(c(g.A.doy,g.B.doy))*1.2
plot(doy,g.A.doy, type="l", lwd=3, lty=3,col=spcol$spA,
     axes=FALSE,ylim=c(-0.05,yup+.05))
 axis(1,at=c(0,doy.max),pos=-0.001,labels=FALSE,lty=1,lwd.ticks=0)
 axis(2,at=c(0,yup),pos=0,labels=FALSE,lwd.ticks=0)
 axis(3,at=c(0,doy.max),pos=yup,labels=FALSE,lty=1,lwd.ticks=0)
 axis(4,at=c(0.001,yup),pos=doy.max,labels=FALSE,lty=1,lwd.ticks=0)
 
points(doy,g.A.cum, type="l",lwd=2, lty=1, col=spcol$spA,
       ylab = "germination", xlab="day of year")
points(doy[d],g.A.doy[d],type="p",pch=20,col=spcol$spA)
points(doy,g.B.doy, type="l", lty=3,lwd=3,col=spcol$spB,)
points(doy,g.B.cum, type="l",lty=1, lwd=2,col=spcol$spB)
points(doy[d],g.B.doy[d],type="p",pch=20,col=spcol$spB)
legend(x=65,y=yup*.65,
       legend=c("sp A daily","sp B daily","sp A cum.","sp cum."),
       col=c(spcol$spA,spcol$spB,spcol$spA,spcol$spB), lty=c(3,3,1,1), lwd=c(2,2,2,2),
       pch=c(20,20,NA,NA),bty="n",cex=0.8)
##Add tauP v doy to top of panel
points(doy,tauP/max(tauP)*.05+yup*1.005,type="l",col=tauPcol, lwd=1,bty="n")
arrows(x0=d*step-1, x1=d*step-1,y0=yup*1.005,y1=tauP[d-1]/max(tauP)*.05+yup*1.005,length=.075,col=tauPcol,cex=0.5)
##Add tauI v doy to bottom of panel
points(doy,-g.A/max(g.A)*.05,type="l",col=spcol$spA, lwd=1)
points(doy,-g.B/max(g.B)*.05,type="l",col=spcol$spB, lwd=1)
points(d*step-1,-g.A[d-1]/max(g.A)*.05,pch=8,col=spcol$spA,cex=0.5)
points(d*step-1,-g.B[d-1]/max(g.B)*.05,pch=8,col=spcol$spB,cex=0.5)
mtext("germination", side=2, line=0)
mtext("day of year",side=1,line=0)

#Middle Right Panel:  Germ Frac v year
par(mar=c(2,2,2,2))
yup <- max(c(g.A.yr,g.B.yr))*1.05
plot(yrs,g.A.yr,type="b",pch=20, col=spcol$spA,
     axes=FALSE, ylim=c(0,yup+.3))
points(yrs,g.B.yr,type="b",pch=20, col=spcol$spB)
#legend("topright",legend=c("sp A", "sp B"),col=c(spcol$spA,spcol$spB),lty=1,pch=20,bty="n")
axis(1,at=c(0,max(yrs)),pos=-0.01,labels=FALSE,lty=1,lwd.ticks=0)
axis(2,at=c(0,yup),pos=0,labels=FALSE,lwd.ticks=0)
axis(3,at=c(0,max(yrs)),pos=yup,labels=FALSE,lty=1,lwd.ticks=0)
axis(4,at=c(-0.01,yup),pos=max(yrs),labels=FALSE,lty=1,lwd.ticks=0)
#plot tauP.t above axis
##Add tauP v y to top of panel
points(yrs,tauP.t/max(tauP.t)*.25+yup*1.005,pch=20,cex=0.5,col=tauPcol)
points(yrs,tauP.t/max(tauP.t)*.25+yup*1.005,type="l",col=tauPcol,lty=1, lwd=1)
mtext("germination",side=2,line=0)
mtext("year",side=1,line=0)

#PANEL for Priority Effect v doy
par(mar=c(2,2,2,2))
yup <- max(c(g.A.priority.cum))*1.05
plot(doy,g.A.priority,type="l",pch=20, lty=3,lwd=3,col=spcol$spA, ylim=c(0,yup),
     axes=FALSE)
points(doy,g.B.priority,type="l",lty=3, lwd=3,col=spcol$spB)
points(doy,g.A.priority.cum,type="l",lty=1,lwd=2,pch=20, col=spcol$spA)
points(doy,g.B.priority.cum,type="l",lty=1,lwd=2,pch=20, col=spcol$spB)

axis(1,at=c(0,doy.max),pos=0,labels=FALSE,lty=1,lwd.ticks=0)
axis(2,at=c(0,yup),pos=0,labels=FALSE,lwd.ticks=0)
axis(3,at=c(0,doy.max),pos=yup-.001,labels=FALSE,lty=1,lwd.ticks=0)
axis(4,at=c(0,yup),pos=doy.max,labels=FALSE,lty=1,lwd.ticks=0)
mtext("day of year",side=1,line=0)
mtext("germination",side=2,line=0)
mtext("cum. germination",side=4,line=0)

#PANEL for Priority Effect v Year
par(mar=c(2,2,2,2))
yup <- max(c(g.B.yr))*1.05
plot(yrs,g.B.yr*.7,type="b",pch=20, lty=1,lwd=1,col=spcol$spA, ylim=c(0,yup+.25),
     axes=FALSE)
points(yrs,g.B.yr,type="b",pch=20, lty=1,lwd=1,col=spcol$spB)
axis(1,at=c(0,max(yrs)),pos=0,labels=FALSE,lty=1,lwd.ticks=0)
axis(2,at=c(0,yup),pos=0,labels=FALSE,lwd.ticks=0)
axis(3,at=c(0,max(yrs)),pos=yup-.001,labels=FALSE,lty=1,lwd.ticks=0)
axis(4,at=c(0,yup),pos=max(yrs),labels=FALSE,lty=1,lwd.ticks=0)
mtext("year",side=1,line=0)
mtext("germination",side=2,line=0)
points(yrs,tauP.t/max(tauP.t)*.2+yup*1.005,pch=20,cex=0.5,col=tauPcol)
points(yrs,tauP.t/max(tauP.t)*.2+yup*1.005,type="l",col=tauPcol,lty=1, lwd=1)

dev.off()



