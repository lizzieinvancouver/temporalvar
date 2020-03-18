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
library(plotrix)
setwd("C:/Users/Megan/Documents/GitHub/temporalvar/R")

#Set up calculations for Panels A = Germination v DOY  & B Germination v Year for germination model
tauP.alpha = 5
tauP.beta = 15
yrs <- seq(1,40,1)
tauP.t <- rbeta(length(yrs), tauP.alpha, tauP.beta)

## tauP
doy.max <- 100
step=.1
doy <- seq(0, doy.max, step)
x <- doy/doy.max
tauP <- dbeta(x, tauP.alpha,tauP.beta)
tauPcol <- "chartreuse3"
tauP.ex = 0.25
d=tauP.ex*100/step     #tauP day for figure

##germination for sp A and B
taui.A <- 0.27
taui.B <- 0.32
spcol <- list(spA =c("darkorange"),spB=c("dodgerblue"))
gmax <- 0.5
h <- 100

#germination distribution --> germination| tauP=x 
g.A.cond <- gmax*exp(-h*(x - taui.A)^2)
g.B.cond <- gmax*exp(-h*(x - taui.B)^2)

#priority effects
g.A.priority <-g.A.cond*.7
g.B.priority <- g.B.cond
g.A.priority.cum <- cumsum(g.A.priority)/sum(g.B.priority)
g.B.priority.cum <- cumsum(g.B.priority)/sum(g.B.priority)
g.A.priority.yr <- rep(g.A.priority.cum[length(x)],length(yrs))
g.B.priority.yr <- rep(g.B.priority.cum[length(x)],length(yrs))

#marginal germination = (germination|tauP=X) * P(tauP)
#Expected germination distribution
g.A <- g.A.cond*tauP
g.B <- g.B.cond*tauP

#time series of germination within a year (t=1)
#germination v time -->  all germination happens on one day, but the amount depends on tauP - tauI
#need a sequence of zeros with germination on, say, day 20

g.A.doy <- c(rep(0,d-1), g.A[d],rep(0,length(doy)-d)) 
g.B.doy <- c(rep(0,d-1), g.B[d],rep(0,length(doy)-d)) 
# g.A.doy <- c(rep(0,d-1), gmax*exp(-h*(tauP.ex - taui.A)^2),rep(0,length(doy)-d))                         
# g.B.doy <- c(rep(0,d-1), gmax*exp(-h*(tauP.ex - taui.B)^2),rep(0,length(doy)-d)) 
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
pdf("graphs/conceptual/PriorityEff_BetHedge.pdf", 
    title="Testing Testing",
    width=8, height=7)
def.par <- par(no.readonly = TRUE) # save default, for resetting...
nf<-layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE),widths=c(4,4),heights=c(2,2,1))
#layout.show(nf)
par(mfrow=c(3,2), oma = c(2, 0, 0, 0))

#Top Left: Priority Effect v doy
par(mar=c(2,2,2,2))
scale=0.1
sos <- 8
yup <- max(c(g.B.priority.cum))*1.05
plot(doy,g.A.priority,type="l",pch=20, lty=2,lwd=3,col=spcol$spA, 
     ylim=c(-scale,yup+scale), axes=FALSE)
points(doy,g.B.priority,type="l",lty=2, lwd=3,col=spcol$spB)
points(doy,g.A.priority.cum,type="l",lty=1,lwd=2,pch=20, col=spcol$spA)
points(doy,g.B.priority.cum,type="l",lty=1,lwd=2,pch=20, col=spcol$spB)
axis(1,at=c(0,doy.max),pos=0,labels=FALSE,lty=1,lwd.ticks=0)
axis(2,at=c(0,yup),pos=0,labels=FALSE,lwd.ticks=0)
axis(3,at=c(0,doy.max),pos=yup-.001,labels=FALSE,lty=1,lwd.ticks=0)
axis(4,at=c(0,yup),pos=doy.max,labels=FALSE,lty=1,lwd.ticks=0)
arrows(x0=sos, x1=sos,y0=scale+yup,y1=yup,
       length=.075,col=tauPcol,cex=0.5)
mtext("Day of Year",side=1,line=0)
mtext("Germination",side=2,line=0)
mtext("Cum. Germination",side=4,line=0)
legend(x=65,y=yup*.45,
       legend=c("Sp1 Daily","Sp2 Daily","Sp1 Cum.","Sp2 Cum."),
       col=c(spcol$spA,spcol$spB,spcol$spA,spcol$spB), lty=c(3,3,1,1), lwd=c(2,2,2,2),
       pch=c(20,20,NA,NA),bty="n",cex=0.8)
text(x=0,y=yup+scale,"A",adj=c(0,0.5),font=2)


#Top Right Panel for Priority Effect v Year
scale=0.1
par(mar=c(2,2,2,2))
yup <- max(c(g.B.priority.yr))*1.05
plot(yrs,g.A.priority.yr,type="l", lty=1,lwd=2,col=spcol$spA, ylim=c(-scale,yup+scale),
     axes=FALSE)
points(yrs,g.B.priority.yr,type="l",lty=1,lwd=2,col=spcol$spB)
points(yrs,g.A.priority.yr,pch=20,col=spcol$spA)
points(yrs,g.B.priority.yr,pch=20,col=spcol$spB)
axis(1,at=c(0,max(yrs)+1),pos=0,labels=FALSE,lty=1,lwd.ticks=0)
axis(2,at=c(0,yup),pos=0,labels=FALSE,lwd.ticks=0)
axis(3,at=c(0,max(yrs)+1),pos=yup-.001,labels=FALSE,lty=1,lwd.ticks=0)
axis(4,at=c(0,yup),pos=max(yrs)+1,labels=FALSE,lty=1,lwd.ticks=0)
mtext("Year",side=1,line=0)
mtext("Germination",side=2,line=0)
points(yrs,rep(tauP.ex,length(yrs))*scale+yup*1.005,pch=20,cex=0.5,col=tauPcol)
points(yrs,rep(tauP.ex,length(yrs))*scale+yup*1.005,type="l",col=tauPcol,lty=1, lwd=1)
text(x=0,y=yup+scale,"B",adj=c(0,0.5),font=2)

#Bottom Left Panel: Germ Frac & Cumulative v DOY 
par(mar=c(2,2,2,2))
scale <- 0.2
yup <- max(c(g.A.doy,g.B.doy))*1.2
plot(doy,g.A.doy, type="l", lwd=3, lty=2,col=spcol$spA,
     axes=FALSE,ylim=c(-scale,yup+scale))
axis(1,at=c(0,doy.max),pos=-0.001,labels=FALSE,lty=1,lwd.ticks=0)
axis(2,at=c(0,yup),pos=0,labels=FALSE,lwd.ticks=0)
axis(3,at=c(0,doy.max),pos=yup,labels=FALSE,lty=1,lwd.ticks=0)
axis(4,at=c(0.001,yup),pos=doy.max,labels=FALSE,lty=1,lwd.ticks=0)
mtext("Cum. Germination",side=4,line=0)
points(doy,g.A.cum, type="l",lwd=2, lty=1, col=spcol$spA,
       ylab = "Germination", xlab="Day of Year")
points(doy[d],g.A.doy[d],type="p",pch=20,col=spcol$spA)
points(doy,g.B.doy, type="l", lty=3,lwd=2,col=spcol$spB,)
points(doy,g.B.cum, type="l",lty=1, lwd=2,col=spcol$spB)
points(doy[d],g.B.doy[d],type="p",pch=20,col=spcol$spB)
##Add tauP v doy to top of panel
points(doy,tauP/max(tauP)*scale+yup*1.005,type="l",col=tauPcol, lwd=1,bty="n")
arrows(x0=d*step, x1=d*step,y0=tauP[d]/max(tauP)*scale+yup*1.005,y1=yup*1.005,
       length=.075,col=tauPcol,cex=0.5)
##Add tauI v doy to bottom of panel
points(doy,-g.A/max(c(g.A,g.B))*scale,type="l",col=spcol$spA, lwd=1)
points(doy,-g.B/max(c(g.A,g.B))*scale,type="l",col=spcol$spB, lwd=1)
arrows(x0=d*step-step/2,x1=d*step-step/2,y0=-g.A[d]/max(c(g.A,g.B))*scale,y1=0,
       length=.075,col=spcol$spA,cex=0.5)
arrows(x0=d*step+step/2,x1=d*step+step/2,y0=-g.B[d]/max(c(g.A,g.B))*scale,y1=0,
       length=.075,col=spcol$spB,cex=0.5)
mtext("Germination", side=2, line=0)
mtext("Day of Year",side=1,line=0)
text(x=0,y=yup+scale,"C",adj=c(0,0.5),font=2)

#Bottom Right:  Germ Frac v year
par(mar=c(2,2,2,2))
scale=0.2
yup <- max(c(g.A.yr,g.B.yr))*1.05
plot(yrs,g.A.yr,type="b",pch=20, col=spcol$spA,
     axes=FALSE, ylim=c(-scale,yup+scale))
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
mtext("Germination",side=2,line=0)
mtext("Year",side=1,line=0)
text(x=0,y=yup+scale,"D",adj=c(0,0.5),font=2)

##Add caption to bottom of pdf
plot(c(0,1),c(0,1),type="n",axes=FALSE, xlab="",ylab="")
textbox(x=c(0,1),y=1,
     textlist=c("Figure 4. Tracking can be conceptualized as changes in priority effects or changes in storage effects.",
     "In a priority effect model (A,B), the coexistence mechanism is a within-year tradeoff between",
     "an early-germinating species that pre-empts resources (sp 1) and late-germinating species that is a superior resource", 
     "competitor (sp 2) (A, where green arrow indicates the start of season); no between-year variation is required",
     "to maintain coexistence. In a storage effect model (C,D), variation in the timing of the start of season",
     "(indicated by the distribution in green, top of C) results results in differential species-response to the environment", 
     "(illustrated by species-specific germination curves, bottom of C); this interannual variation in", 
     "species-response to the environment (D) - along with a seedbank or other interannual storage mechanism - can maintain coexistence",
     "through reduced interspecific competition."),box=FALSE)

dev.off()



