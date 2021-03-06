## Started in March 2018 ##
## By Lizzie ##

## Plots of interannual dynamics ##
## Code taken from PhenologyModelAnalysisOLD.R ##

## This code makes Bfin by tauP plots and g_i dynamics plots ##

## Right now running models with no coexisting runs (otheruns) is a bit cheap ...#
# You have to change sruns -> otheruns in two places and ...
# change runningotheruns to TRUE ## 

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
runbfin <- TRUE # Note that reading in these files is SLOW!
runnonstat <- FALSE # if you are using the nsruns list, this will check for coexistence at end of stationary
runningotheruns <- FALSE # since there are no coexisting runs to ref...
# this just pulls the Bout files I have (see below)

# cheap loop over the files for now
# sruns <- c("36691943")
sruns <- c("36426477", "36511349","36691943","36511352", "36511384")
otheruns <- c("36691954") # have not pulled 36691955 or 36691956

if(runnonstat){
nsruns <- c("36511352", "36511384", "36691954")
nsruns.yrs <- c(500, 1000, 1000) # for nsruns, how many initial stationary years 
}


runlist <- sruns # sruns <- c("36691954") 

for (folderIDhere in c(1:length(runlist))){
    
folderID <- runlist[folderIDhere] # 36511349
samplerun <-  read.table(paste("output/", folderID, "/SummaryOut_", folderID,
    "-1.txt", sep=""), header=TRUE)
filenamestart <- c(paste("SummaryOut_", folderID, "-", sep=""))
colnameshere <- colnames(samplerun)
numhere <- c(1:20)

runs <- getfiles(folderID, filenamestart, numhere, colnameshere)
runs$taskrunID <- paste(runs$taskID, runs$runID, sep="-")

# Get bfin (slow!) 
if(runbfin){
samplerunbfin <-  read.table(paste("output/", folderID, "/BfinN_", folderID,
    "-1.txt", sep=""), header=TRUE)
filenamestartBfin <- c(paste("BfinN_", folderID, "-", sep=""))
colnamesherebfin <- colnames(samplerunbfin)

runsbfin <- getfiles(folderID, filenamestartBfin, numhere, colnamesherebfin)
runsbfin$taskrunID <- paste(runsbfin$taskID, runsbfin$runID, sep="-")
}

# Get envt params (also slow)
samplerunep<-  read.table(paste("output/", folderID, "/EnvtParms_", folderID,
    "-1.txt", sep=""), header=TRUE)
filenamestartep <- c(paste("EnvtParms_", folderID, "-", sep=""))
colnameshereep <- colnames(samplerunep)

runsep <- getfiles(folderID, filenamestartep, numhere, colnameshereep)
runsep$taskrunID <- paste(runsep$arrayID, runsep$runID, sep="-")

# Finish making taskrun and taskrunyr
runsbfin$taskrunIDyr <- paste(runsbfin$taskrunID, runsbfin$yr, sep="-")
runsep$taskrunIDyr <- paste(runsep$taskrunID, runsep$yr, sep="-")
runs$taskrunID <- paste(runs$taskID, runs$runID, sep="-")

##
## Get list of coexisting
##

if(runnonstat){
runs$stat.ncoexist <- NA
runs$stat.coexist1 <- NA
runs$stat.coexist2 <- NA
idshere <- unique(runsbfin$taskrunID)

for (onerun in c(1:length(idshere))){
    subby <- subset(runsbfin, taskrunID==idshere[onerun])
    ifelse(nrow(subby)<nsruns.yrs[folderIDhere],
           runs$stat.coexist1[runs$taskrunID==idshere[onerun]] <- subby$Bfin1[nrow(subby)],
       ifelse(subby$Bfin1[nsruns.yrs[folderIDhere]]>0,
           runs$stat.coexist1[runs$taskrunID==idshere[onerun]] <- 1,
           runs$stat.coexist1[runs$taskrunID==idshere[onerun]] <- 0))
    ifelse(nrow(subby)<nsruns.yrs[folderIDhere],
           runs$stat.coexist2[runs$taskrunID==idshere[onerun]] <- subby$Bfin2[nrow(subby)],
       ifelse(subby$Bfin2[nsruns.yrs[folderIDhere]]>0,
           runs$stat.coexist2[runs$taskrunID==idshere[onerun]] <- 1,
           runs$stat.coexist2[runs$taskrunID==idshere[onerun]] <- 0))
     }                         
runs$stat.ncoexist <- runs$stat.coexist1+runs$stat.coexist2
table(runs$stat.ncoexist)
}

df <- makediffs(runs)
df.coexist <- subset(df, ncoexist==2)
df.nocoexist <- subset(df, ncoexist<2)

##
## Checking R0 and run lengths
##
hiLtstp <- subset(runsbfin, Ltstp==1001)
lowLtstp <- subset(runsbfin, Ltstp<900)

envt.hiLtstp <- runsep[which(runsep$taskrunIDyr %in% hiLtstp$taskrunIDyr),]
envt.lowLtstp <- runsep[which(runsep$taskrunIDyr %in% lowLtstp$taskrunIDyr),]

par(mfrow=c(1,2))
hist(envt.hiLtstp$R0, main="Ltstp is max")
hist(envt.lowLtstp$R0, main="Ltstp < 900")

##############################
# Plot a few coexisting runs #
##############################

## Step 1: check out a couple coexisting runs
if(runningotheruns){
runstouse <- c("1-4", "1-46", "1-54", "1-107", "1-154") # for 36691954 only
}

if(!runningotheruns){
runstouse <- df.coexist$taskrunID[1:6]
}

runbfinhere <- runsbfin
bfincoexist <- runsbfin[which(runbfinhere$taskrunID %in% runstouse),]
epcoexist <- runsep[which(runsep$taskrunID %in% runstouse),]

    
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

## Step 2: g_i ~ eff|tauI-tauP|
# for a couple coexisting runs and a couple not coexisting
gi.runstouse <- runsep[which(runsep$taskrunID %in% df.coexist$taskrunID[1:3]),]
gi.runstouse.nocoexist <- runsep[which(runsep$taskrunID %in% df.nocoexist$taskrunID[1:3]),]


xlim <- c(0,1)

pdf(paste("graphs/modelruns/gi_dynamics/runs_", folderID, "_gi_3coexistruns.pdf", sep=""), width=5, height=7)
par(mfrow=c(3,2))
for (whichrun in seq_along(unique(gi.runstouse$taskrunID))){
    runhere <- unique(gi.runstouse$taskrunID)[whichrun]
    dathere <- subset(gi.runstouse, taskrunID==runhere)
    colorz <- c("black", "blue")
    ylab <- "gi"
    xlab <- "tauIhat"
    plot(dathere$g1~dathere$tauIhat1, ylab=ylab,
        xlab=xlab, xlim=xlim)
    points(dathere$g2~dathere$tauIhat2, col=colorz[2],
        xlim=xlim)
    hist(dathere$tauP, main=paste("run", runhere), xlim=xlim)
}
dev.off()

pdf(paste("graphs/modelruns/gi_dynamics/runs_", folderID, "_gi_3nocoexistruns.pdf", sep=""), width=5, height=7)
par(mfrow=c(3,2))
for (whichrun in seq_along((unique(gi.runstouse.nocoexist$taskrunID)))){
    runhere <- unique(gi.runstouse.nocoexist$taskrunID)[whichrun]
    dathere <- subset(gi.runstouse.nocoexist, taskrunID==runhere)
    colorz <- c("black", "blue")
    ylab <- "gi"
    xlab <- "tauIhat"
    plot(dathere$g1~dathere$tauIhat1, ylab=ylab,
        xlab=xlab, xlim=xlim)
    points(dathere$g2~dathere$tauIhat2, col=colorz[2],
        xlim=xlim)
    hist(dathere$tauP, main=paste("run", runhere), xlim=xlim)
}
dev.off()

}
