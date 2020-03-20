## Started 20 March 2020 ##
## By Lizzie ##

## Quick plots for environmental tracking definition ##
if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/temporalvar/R") 
} else setwd("~/Documents/git/boopboop")


# Normal distribution for fundamental tracking
pdf(paste("graphs/conceptual/normal.pdf", sep=""), width = 7.5, height = 7)
x <- seq(-12, 12, length=1000)
hx <- dnorm(x, mean=0, sd=3)
plot(x, hx, type="l", lwd=2,
  ylab="Total fitness", xlab="(Event day - ideal day)")
abline(v=0, lty=2, lwd=2)
dev.off()

pdf(paste("graphs/conceptual/normalscatter.pdf", sep=""), width = 7.5, height = 7)
xnew <- seq(0, 30, length=30)
ynew <- xnew + rnorm(length(xnew), 0, 0.1)
plot(rev(xnew), ynew, pch=16,
   xlab="Ideal day", ylab="Event day", cex=1.7)
abline(lm(rev(xnew)~ynew))
dev.off()

# Some dead simple scatterplots
cues <- data.frame(thermsum=seq(200,500, length.out=300), photo=seq(6,18,
    length.out=300), doy=seq(60,90, length.out=300))
xhere <- seq(8, 12, length.out=41)
doyhere <- 120 - 2*xhere + rnorm(length(xhere), 0, 4)
meas <- data.frame(year=seq(1980, 2020), doy=doyhere)

pdf(paste("graphs/conceptual/simplescatterplots.pdf", sep=""), width = 14, height = 5)
par(mfrow=c(1,3))
plot(cues$photo, cues$thermsum, type="n",
   xlab="Daylength", ylab="Required thermal sum")
plot(cues$photo, cues$doy, type="n",
  xlab="Cue system in environment", ylab="Event day")
plot(meas$year, meas$doy, pch=16,
  xlab="Metric of environment (e.g., MST)", ylab="Event day", cex=1.7)
abline(lm(doy~year, data=meas))
dev.off()
