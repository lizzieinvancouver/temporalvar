## Started in March 2018 ##
## By Lizzie ##

## Plots of interannual dynamics ##

## Code taken from PhenologyModelAnalysisOLD.R ##
## See START HERE below for what to do next ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## libraries
library(ggplot2)
require(plyr); require(dplyr); require(tidyr) # data formatting

## set working directory
setwd("~/Documents/git/projects/temporalvar/R")

source("sourcefiles/multiplot.R")
source("sourcefiles/runanalysisfxs.R")

## flags for what to do
# runbfin <- TRUE 

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
samplerunep<-  read.table(paste("output/", folderID, "/EnvtParms_", folderID,
    "-1.txt", sep=""), header=TRUE)
filenamestartep <- c(paste("EnvtParms_", folderID, "-", sep=""))
colnameshereep <- colnames(samplerunep)

runsep <- getfiles(folderID, filenamestartep, numhere, colnameshereep)
runsep$taskrunID <- paste(runsep$arrayID, runsep$runID, sep="-")


# Get Bout (slower) ## START HERE! Can I make into DF??
getBoutfiles <- function(folderID, boutfilenamestart){
    boutfilenamestart <- 
    filename <- paste("output/", folderID, "/", "Bout/", "Bout*", ".txt", sep="")
    dat <- lapply(Sys.glob(filename), function(i) read.table(i, header=TRUE))
    datahere <- do.call("rbind", dat)
    return(dat)
}

runsbout <- getBoutfiles(folderID, "Bout_36511349")


##
## Get list of coexisting
##

df <- makediffs(runs)
df$taskrunID <- paste(df$taskID, df$runID, sep="-")
df.coexist <- subset(df, ncoexist==2)

##
## Formatting task-run-yr
##
runsbfin$taskrunIDyr <- paste(runsbfin$taskrunID, runsbfin$yr, sep="-")
runsep$taskrunIDyr <- paste(runsep$taskrunID, runsep$yr, sep="-")


##############################
# Plot a few coexisting runs #
##############################

## Step 1: check out a couple coexisting runs
runbfinhere <- runsbfin
bfincoexist <- runsbfin[which(runbfinhere$taskrunID %in% df.coexist$taskrunID),]
epcoexist <- runsep[which(runsep$taskrunID %in% df.coexist$taskrunID),]

runstouse <- df.coexist$taskrunID[1:6]

pdf(paste("graphs/modelruns/interannual/runs_", folderID, "_coexist_sampleruns.pdf", sep=""), width=5, height=7)
par(mfrow=c(3,2))

for (whichrunnum in c(1:length(runstouse))){
whichrun <- runstouse[whichrunnum]
yhere <- c(-6, 3)
colorz <- c("black", "red", "blue")

ephere <- subset(epcoexist, taskrunID==whichrun)
bfinhere <- subset(bfincoexist, taskrunID==whichrun)

xhere <- c(1:length(ephere$tauP))
xheretauI <- c(1:nrow(bfinhere))
xname <- "time"
yname <- "tauP or Bfin"
leg.txt <- c("tauP", "Bfin sp1", "Bfin sp2")
lty.here <- c(rep(1, 3))

plot(ephere$tauP~xhere, type="l", ylim=yhere, col=colorz[1], lty=lty.here[1],
    xlab=xname, ylab=yname, main=paste("run ", whichrun))
lines(log10(bfinhere$Bfin1)~xheretauI, ylim=yhere, col=colorz[2], lty=lty.here[2],
    xlab=xname, ylab=yname, lwd=2)
lines(log10(bfinhere$Bfin2)~xheretauI, ylim=yhere, col=colorz[3], lty=lty.here[3],
    xlab=xname, ylab=yname)
legend("topright", leg.txt, lwd=2, col=colorz)

yhere <- c(0, 1)
yname <- "tauP or tauIhat"
leg.txt <- c("tauP", "tauIhat 1", "tauIhat 2")
lty.here <- c(rep(1, 1), rep(2, 2))

plot(ephere$tauP~xhere, type="l", ylim=yhere, col=colorz[1], lty=lty.here[1],
    xlab=xname, ylab=yname)
lines(ephere$tauIhat1~xhere, ylim=yhere, col=colorz[2], lty=lty.here[2],
    xlab=xname, ylab=yname, lwd=2)
lines(ephere$tauIhat2~xhere, ylim=yhere, col=colorz[3], lty=lty.here[3],
    xlab=xname, ylab=yname)
legend("topright", leg.txt, lwd=2, lty=lty.here, col=colorz)
}

dev.off()

## START HERE ##
## Step X: g_i ~ eff|tauI-tauP|
# for a couple coexisting runs and a couple not coexisting

gi.runstouse <- as.numeric(row.names(subset(as.data.frame(whocoexisted), whocoexisted>1))[1:3])
gi.runstouse.nocoexist <- as.numeric(row.names(subset(as.data.frame(whocoexisted),
    whocoexisted==1))[1:3])

xlim <- c(0,1)

pdf(paste("graphs/modelruns/gi_dynamics/Track_varR_2spp_", jobIDname, "_gi_3coexistruns.pdf", sep=""), width=5, height=7)
par(mfrow=c(3,2))
for (whichrun in seq_along(gi.runstouse)){
    runhere <- gi.runstouse[whichrun]
    colorz <- c("black", "blue")
    ylab <- "gi"
    xlab <- "tauIhat"
    plot(modelruns[[runhere]]$g[,1]~modelruns[[runhere]]$tauIhat[,1], ylab=ylab,
        xlab=xlab, xlim=xlim)
    points(modelruns[[runhere]]$g[,2]~modelruns[[runhere]]$tauIhat[,2], col=colorz[2],
        xlim=xlim)
    hist(modelruns[[runhere]]$envtvars[["tauP"]], main=paste("run", runhere), xlim=xlim)
}
dev.off()

pdf(paste("graphs/modelruns/gi_dynamics/Track_varR_2spp_", jobIDname, "_gi_3nocoexistruns.pdf", sep=""), width=5, height=7)
par(mfrow=c(3,2))
for (whichrun in seq_along(gi.runstouse.nocoexist)){
    runhere <- gi.runstouse.nocoexist[whichrun]
    colorz <- c("black", "blue")
    ylab <- "gi"
    xlab <- "tauIhat"
    plot(modelruns[[runhere]]$g[,1]~modelruns[[runhere]]$tauIhat[,1], ylab=ylab,
        xlab=xlab, xlim=xlim)
    points(modelruns[[runhere]]$g[,2]~modelruns[[runhere]]$tauIhat[,2], col=colorz[2],
        xlim=xlim)
    hist(modelruns[[runhere]]$envtvars[["tauP"]], main=paste("run", runhere), xlim=xlim)
}
dev.off()
