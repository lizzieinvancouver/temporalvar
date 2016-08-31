### Started 22 August 2016 ###
### By Lizzie ###

## Analyzing the runs data for the storage effect model! ##

## Questions: ##
# None just now! #


## To do ##
# Change all the file paths so they use the jobID? #

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## libraries
library(ggplot2)
require(plyr); require(dplyr); require(tidyr) # data formatting

## set working directory
setwd("~/Documents/git/projects/temporalvar/R")

load('..//notposting/workon/Track_varR_2spp_60120935.Rdata')

# Here's the structure:
# list(jobID=jobID, arrayNum=a, runNum=r, sppvars=sppvars,
#    tauI=tauI, tauP=tauP, Bfin=Bfin,Bout=Bout (within Bout is time, R, B1, B2) )

## some plots

## first, check out a couple runs
runstouse <- c(1, 777, 851)

pdf("graphs/Track_varR_2spp_60120935_3runs.pdf", width=5, height=7)
par(mfrow=c(3,2))

for (whichrun in seq_along(runstouse)){
yhere <- c(0, 65)
colorz <- c("black", "red", "blue")

xhere <- c(1:length(modelruns[[1]]$tauP))
xname <- "time"
yname <- "tauP or Bfin"
leg.txt <- c("tauP", "Bfin sp1", "Bfin sp2")
lty.here <- c(rep(1, 3))

plot(modelruns[[whichrun]]$tauP~xhere, type="l", ylim=yhere, col=colorz[1], lty=lty.here[1],
    xlab=xname, ylab=yname, main=paste("run ", whichrun))
lines(modelruns[[whichrun]]$Bfin[,1]~xhere, ylim=yhere, col=colorz[2], lty=lty.here[2],
    xlab=xname, ylab=yname)
lines(modelruns[[whichrun]]$Bfin[,2]~xhere, ylim=yhere, col=colorz[3], lty=lty.here[3],
    xlab=xname, ylab=yname)
legend("topright", leg.txt, lwd=2, col=colorz)

yhere <- c(0, 1)
yname <- "tauP or tauI"
leg.txt <- c("tauP", "tauI 1", "tauI 2")
lty.here <- c(rep(1, 1), rep(2, 2))

plot(modelruns[[whichrun]]$tauP~xhere, type="l", ylim=yhere, col=colorz[1], lty=lty.here[1],
    xlab=xname, ylab=yname)
lines(modelruns[[whichrun]]$tauI[,1]~xhere, ylim=yhere, col=colorz[2], lty=lty.here[2],
    xlab=xname, ylab=yname)
lines(modelruns[[whichrun]]$tauI[,2]~xhere, ylim=yhere, col=colorz[3], lty=lty.here[3],
    xlab=xname, ylab=yname)
legend("topright", leg.txt, lwd=2, lty=lty.here, col=colorz)
}

dev.off()

## now the next fig from 20160420_2sppCoexistFigs.pdf
par(mfrow=c(1,2))

whichrun.now <- 1
plot(modelruns[[whichrun.now]]$tauI[,1]~modelruns[[whichrun.now]]$tauI[,2])
whichrun.now <- 111
plot(modelruns[[whichrun.now]]$tauI[,1]~modelruns[[whichrun.now]]$tauI[,2])


## Lizzie's random notes and such on how to use lists
# A never-ending battle ...

modelruns[[1]][["sppvars"]]$alpha

# this does something
lapply(modelruns, "[", c("sppvars"))

## from Dan:
# missing extra bracket! sppvars is the list entry name, so need double brackets around it

lapply(modelruns, function(x) x[["sppvars"]]$alpha)

# unstructured data
unlist(lapply(modelruns, function(x) x[["sppvars"]]$alpha))

# structured in data frame
xx <- data.frame(lapply(modelruns, function(x) x[["sppvars"]]$alpha))

(names(xx) <- paste("run", 1:ncol(xx), sep=""))

# structured with runs as rows
xx <- t(xx)
colnames(xx) <- paste("a", 1:2, sep = "")
xx

##
# me working
goo <- unlist(lapply(modelruns, function(x) x[["sppvars"]]$tauIPini))
