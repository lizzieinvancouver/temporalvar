## Started in March 2018 ##
## By Lizzie ##

## Plots of interannual dynamics ##

## Code taken from PhenologyModelAnalysisOLD.R ##
## NEXT UP! Find out where are the COEXISTING runs! ##

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
filenamestartbout <- c(paste("Bout_", folderID, "-1-", sep=""))
colnamesherebout <- colnames(samplerunbout)
numherebout <- c(1:13)

runsbout <- getBoutfiles(folderID, filenamestartbout, numherebout, colnamesherebout)
runsbout$taskrunID <- paste(runsbout$taskID, runsbout$runID, sep="-")

##
## Get list of coexisting
##

df <- makediffs(runs)
df$taskrunID <- paste(df$taskID, df$runID, sep="-")
df.coexist <- subset(df, ncoexist==2)

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

## Check out a couple runs
# Note: Not sure that I have any that co-exist....
runbouthere <- runsbout
boutcoexist <- runsbout[which(runsbout$taskrunID %in% df.coexist$taskrunID),]
if(getenvrt){
epcoexist <- runsep[which(runsep$taskrunID %in% df.coexist$taskrunID),]
}

runstouse <- df$taskrunID[1:10]

pdf(paste("graphs/modelruns/intraannual/runs_", folderID, "_sampleruns.pdf", sep=""), width=10, height=7)
# par(mfrow=c(3,2))

for (whichrunnum in c(1:length(runstouse))){
whichrun <- runstouse[whichrunnum]
yhere <- c(-20, 3)
colorz <- c("black", "red", "blue")

bouthere <- subset(runbouthere, taskrunID==whichrun)

xhere <- c(1:nrow(bouthere))
xname <- "time"
yname <- "R or Bout"
leg.txt <- c("R", "Bout sp1", "Bout sp2")
lty.here <- c(rep(1, 3))

plot(bouthere$R~xhere, type="l", ylim=yhere, col=colorz[1], lty=lty.here[1],
    xlab=xname, ylab=yname, main=paste("run ", whichrun))
lines(log10(as.numeric(bouthere$B1))~xhere, ylim=yhere, col=colorz[2], lty=lty.here[2],
    xlab=xname, ylab=yname, lwd=2)
lines(log10(as.numeric(bouthere$B2))~xhere, ylim=yhere, col=colorz[3], lty=lty.here[3],
    xlab=xname, ylab=yname)
legend("topright", leg.txt, lwd=2, col=colorz)
}

dev.off()
