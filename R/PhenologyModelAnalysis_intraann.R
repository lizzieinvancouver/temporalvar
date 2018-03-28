## Started in March 2018 ##
## By Lizzie ##

## Plots of interannual dynamics ##

## Code taken from PhenologyModelAnalysisOLD.R ##
## NEXT UP! Pull the Bout runs from the subsamples ... #
# ... and adjust code so that it runs over the folders ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## libraries
library(ggplot2)
require(plyr); require(dplyr); require(tidyr) # data formatting

## set working directory
setwd("~/Documents/git/projects/temporalvar/R")

source("sourcefiles/analyses/multiplot.R")
source("sourcefiles/analyses/runanalysisfxs.R")

## flags for what to do
getenvrt <- FALSE 

# cheap loop over the files for now
sruns <- c("36426477", "36511349","36691943", "36691954", "36691955")
nsruns <- c("36511352", "36511384", "36691956")

folderID <- 36511349
samplerun <-  read.table(paste("output/", folderID, "/SummaryOut_", folderID,
    "-1.txt", sep=""), header=TRUE)
filenamestart <- c(paste("SummaryOut_", folderID, "-", sep=""))
colnameshere <- colnames(samplerun)
numhere <- c(1:20)

runs <- getfiles(folderID, filenamestart, numhere, colnameshere)
runs$taskrunID <- paste(runs$taskID, runs$runID, sep="-")

# and get the subsample file ... (for which runs are pulled)
# note that by default this run includes 10 coexisting runs (first 10) and 10 not coexisting
subsample <-  read.csv(paste("output/", folderID, "/subsample.csv", sep=""), header=TRUE)
subsample$taskrunID <- paste(subsample$taskID, subsample$runID, sep="-")


# Get envt params (slow)
if(getenvrt){
samplerunep<-  read.table(paste("output/", folderID, "/EnvtParms_", folderID,
    "-1.txt", sep=""), header=TRUE)
filenamestartep <- c(paste("EnvtParms_", folderID, "-", sep=""))
colnameshereep <- colnames(samplerunep)

runsep <- getfiles(folderID, filenamestartep, numhere, colnameshereep)
runsep$taskrunID <- paste(runsep$arrayID, runsep$runID, sep="-")
}

# Get Bout (slower) 
samplerunbout <-  read.table(paste("output/", folderID, "/Bout_", folderID,
    "-1-1.txt", sep=""), header=TRUE)
filenamestartbout <- c(paste("Bout_", folderID, "-", sep=""))
colnamesherebout <- colnames(samplerunbout)
numherebout <- subsample$taskrunID

runsbout <- getBoutfiles(folderID, filenamestartbout, numherebout, colnamesherebout)
runsbout$taskrunID <- paste(runsbout$taskID, runsbout$runID, sep="-")

##
## Formatting task-run-yr
##
if(getenvrt){
runsep$taskrunIDyr <- paste(runsep$taskrunID, runsep$yr, sep="-")
}
runsbout$taskrunIDyr <- paste(runsbout$taskrunID, runsbout$yr, sep="-")


###################
# Plot a few runs #
###################
runbouthere <- runsbout
boutcoexist <- runsbout[which(runsbout$taskrunID %in% subsample$taskrunID[1:10]),]
if(getenvrt){
epcoexist <- runsep[which(runsep$taskrunID %in% subsample$taskrunID[1:10]),]
}

##
## Check out a couple coexisting runs
## 
runstouse <- subsample$taskrunID[1:5]

pdf(paste("graphs/modelruns/intraannual/runs_", folderID, "_coexist_sampleruns.pdf", sep=""), width=10, height=7)
# par(mfrow=c(3,2))

for (whichrunnum in c(1:length(runstouse))){
whichrun <- runstouse[whichrunnum]
yhere <- c(-10, 3)
colorz <- c("black", "red", "blue")

bouthere <- subset(runbouthere, taskrunID==whichrun)

xhere <- c(1:nrow(bouthere))
xname <- "time"
yname <- "R or Bout"
leg.txt <- c("R", "Bout sp1", "Bout sp2")
lty.here <- c(rep(1, 3))

# plot the whole time series
plot(bouthere$R~xhere, type="l", ylim=yhere, col=colorz[1], lty=lty.here[1],
    xlab=xname, ylab=yname, main=paste("run ", whichrun))
lines(log10(as.numeric(bouthere$B1))~xhere, ylim=yhere, col=colorz[2], lty=lty.here[2],
    xlab=xname, ylab=yname, lwd=2)
lines(log10(as.numeric(bouthere$B2))~xhere, ylim=yhere, col=colorz[3], lty=lty.here[3],
    xlab=xname, ylab=yname)
legend("topright", leg.txt, lwd=2, col=colorz)

# plot a shorter part time series
howshort <- as.integer(nrow(bouthere)/20)
xhereshort <- c(1:howshort)
plot(bouthere$R[1:howshort]~xhereshort, type="l", ylim=yhere, col=colorz[1], lty=lty.here[1],
    xlab=xname, ylab=yname, main=paste("run ", whichrun))
lines(log10(as.numeric(bouthere$B1[1:howshort]))~xhereshort, ylim=yhere, col=colorz[2], lty=lty.here[2],
    xlab=xname, ylab=yname, lwd=2)
lines(log10(as.numeric(bouthere$B2[1:howshort]))~xhereshort, ylim=yhere, col=colorz[3], lty=lty.here[3],
    xlab=xname, ylab=yname)
}

dev.off()


##
## And now check out a couple coexisting runs
## 
runstouse.nocoexist <- subsample$taskrunID[11:14]

pdf(paste("graphs/modelruns/intraannual/runs_", folderID, "_nocoexist_sampleruns.pdf", sep=""), width=10, height=7)
# par(mfrow=c(3,2))

for (whichrunnum in c(1:length(runstouse.nocoexist))){
whichrun <- runstouse.nocoexist[whichrunnum]
yhere <- c(-15, 3)
colorz <- c("black", "red", "blue")

bouthere <- subset(runbouthere, taskrunID==whichrun)

xhere <- c(1:nrow(bouthere))
xname <- "time"
yname <- "R or Bout"
leg.txt <- c("R", "Bout sp1", "Bout sp2")
lty.here <- c(rep(1, 3))

# plot the whole time series
plot(bouthere$R~xhere, type="l", ylim=yhere, col=colorz[1], lty=lty.here[1],
    xlab=xname, ylab=yname, main=paste("run ", whichrun))
lines(log10(as.numeric(bouthere$B1))~xhere, ylim=yhere, col=colorz[2], lty=lty.here[2],
    xlab=xname, ylab=yname, lwd=2)
lines(log10(as.numeric(bouthere$B2))~xhere, ylim=yhere, col=colorz[3], lty=lty.here[3],
    xlab=xname, ylab=yname)
legend("topright", leg.txt, lwd=2, col=colorz)
}

dev.off()
