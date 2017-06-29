### Started 22 August 2016 ###
### By Lizzie ###

## Analyzing the runs data for the storage effect model! ##
## Last updated 9 January 2017 (flight from Oahu) ##

## Notes ##
## This is the initial file for analyzing the model output ##
## It analyzes runs from January-May 2017 ##
## After that we changed the length of runs and some other factors ##
## Keeping this in case ever needed ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## libraries
library(ggplot2)
require(plyr); require(dplyr); require(tidyr) # data formatting

## set working directory
setwd("~/Documents/git/projects/temporalvar/R")

# source some stuff
source("sourcefiles/multiplot.R")

jobID <- "87783087" # "78476247"
load(paste("..//ModelRuns/Track_varR_2spp_", jobID, ".Rdata", sep=""))

# Checking the 0 track runs.... 
# load('..//ModelRuns/Track0_varR_2spp_87839120.Rdata')
# load('..//ModelRuns/Track0_varR_2spp_88003925.Rdata')

# Here's the structure:
# list(jobID=jobID, arrayNum=a, runNum=r, sppvars=sppvars (the unchanging ones and means of the changing ones),
#    envtvars=envtvars (within which is R0, tauP, eps),
#    tauIhat, Bfin=Bfin, g,
#    Bout=Bout (within Bout is time, R, B1, B2),
# To look at the structure try:
# head(modelruns[[1]])
# ... though, err, not always working so also check the top of what you get with:
# str(modelruns)
# tauIPini: tauP-tauIhat using initial tauP (from stationary period)
# tauIPns:  tauP-tauIhat using tauP from non-stationary period
# tauIPfin: tauP-tauIhat using last stage of tauP (from second stationary period, when there is one)

# checking out Rstar ratios of the two species model
rstarz <- c()
for (i in c(1:length(modelruns))){
   rstarz[i] <- modelruns[[i]]$sppvars[["Rstar"]][1]/modelruns[[i]]$sppvars[["Rstar"]][2]
}
hist(rstarz, breaks=100)

###############################################
### Setting up some basic stuff we may want ###
###############################################

# Checking out which runs coexisted
# note: last timestep of Bfin is 199
whocoexisted <- c()
for (j in c(1:length(modelruns))){
    if ((modelruns[[j]]$Bfin[199,1]>0)==TRUE &
        (modelruns[[j]]$Bfin[199,2]>0)==TRUE)
        whocoexisted[j] <- 2
    if(
        (modelruns[[j]]$Bfin[199,1]>0)==FALSE &
        (modelruns[[j]]$Bfin[199,2]>0)==TRUE)
        whocoexisted[j] <- 1
    if(
        (modelruns[[j]]$Bfin[199,1]>0)==TRUE &
        (modelruns[[j]]$Bfin[199,2]>0)==FALSE)
        whocoexisted[j] <- 1
    if(
        (modelruns[[j]]$Bfin[199,1]>0)==FALSE &
        (modelruns[[j]]$Bfin[199,2]>0)==FALSE)
        whocoexisted[j] <- 0
    }
        
length(whocoexisted[which(whocoexisted>1)])
# preserve row numbers: subset(as.data.frame(whocoexisted), whocoexisted>1)

# make a dataframe with necessary info
# this isn't really necessary but c'mon, it makes quick plotting easier!

df <- data.frame(modelrun=rep(0, length(modelruns)),
    tauI.sp1=rep(0, length(modelruns)), tauI.sp2=rep(0, length(modelruns)),
    alpha.sp1=rep(0, length(modelruns)), alpha.sp2=rep(0, length(modelruns)), 
    tauIhat.sp1=rep(0, length(modelruns)), tauIhat.sp2=rep(0, length(modelruns)),
    tauIPini.sp1=rep(0, length(modelruns)), tauIPini.sp2=rep(0, length(modelruns)),
    gi.sp1= rep(0, length(modelruns)), gi.sp2= rep(0, length(modelruns)),
    c.sp1=rep(0, length(modelruns)), c.sp2=rep(0, length(modelruns)),
    rstar.sp1=rep(0, length(modelruns)), rstar.sp2=rep(0, length(modelruns)),
    coexist=whocoexisted)

## impt note! This only uses tauIP initial just now (ignores nonstationary runs)
# tauIPini is defined in getSpecies.R and is the average difference between tauP and ...tauIhat for stationary perion
for (k in c(1:length(modelruns))){
    df$modelrun[k] <- k
    df$tauI.sp1[k] <- modelruns[[k]]$sppvars[["tauI"]][1]
    df$tauI.sp2[k] <- modelruns[[k]]$sppvars[["tauI"]][2]
    df$alpha.sp1[k] <- modelruns[[k]]$sppvars[["alpha"]][1]
    df$alpha.sp2[k] <- modelruns[[k]]$sppvars[["alpha"]][2]
    df$tauIhat.sp1[k] <- mean(modelruns[[k]]$tauIhat[,1])
    df$tauIhat.sp2[k] <- mean(modelruns[[k]]$tauIhat[,2])
    df$tauIPini.sp1[k] <- modelruns[[k]]$sppvars[["tauIPini"]][1]
    df$tauIPini.sp2[k] <- modelruns[[k]]$sppvars[["tauIPini"]][2]
    df$gi.sp1[k] <- mean(modelruns[[k]]$g[,1])
    df$gi.sp2[k] <- mean(modelruns[[k]]$g[,2])
    df$c.sp1[k] <- modelruns[[k]]$sppvars[["c"]][1]
    df$c.sp2[k] <- modelruns[[k]]$sppvars[["c"]][2]
    df$rstar.sp1[k] <- modelruns[[k]]$sppvars[["Rstar"]][1]
    df$rstar.sp2[k] <- modelruns[[k]]$sppvars[["Rstar"]][2]
    }

df$coexist <- as.factor(df$coexist)
df$diff.tauI <- df$tauI.sp1-df$tauI.sp2
df$diff.alpha <- df$alpha.sp1-df$alpha.sp2
df$diff.tauIhat <- df$tauIhat.sp1-df$tauIhat.sp2
df$diff.tauIPini <- df$tauIPini.sp1-df$tauIPini.sp2
df$diff.gi <- df$gi.sp1-df$gi.sp2
df$diff.c <- df$c.sp1-df$c.sp2
df$diff.rstar <- df$rstar.sp1-df$rstar.sp2
df$rstar.ratio <- df$rstar.sp1/df$rstar.sp2
df$tauI.ratio <- df$tauI.sp1/df$tauI.sp2
df$alpha.ratio <- df$alpha.sp1/df$alpha.sp2
df$tauIhat.ratio <- df$tauIhat.sp1/df$tauIhat.sp2
df$tauIPini.ratio <- df$tauIPini.sp1/df$tauIPini.sp2
df$diff.gi <- df$gi.sp1-df$gi.sp2
df$diff.c <- df$c.sp1-df$c.sp2


###############################################
## A few plots ##
## (working off 20160420_2sppCoexistFigs.pdf)
###############################################

## Step 1: check out a couple runs
runstouse <- c(1, 8, 110)

pdf(paste("graphs/modelruns/Track_varR_2spp_", jobID, "_3sampleruns.pdf", sep=""), width=5, height=7)
par(mfrow=c(3,2))

for (whichrun in seq_along(runstouse)){
yhere <- c(0, 65)
colorz <- c("black", "red", "blue")

xhere <- c(1:length(modelruns[[1]]$envtvars[["tauP"]]))
xheretauI <- c(1:length(modelruns[[1]]$Bfin[,1]))
xname <- "time"
yname <- "tauP or Bfin"
leg.txt <- c("tauP", "Bfin sp1", "Bfin sp2")
lty.here <- c(rep(1, 3))

plot(modelruns[[whichrun]]$envtvars[["tauP"]]~xhere, type="l", ylim=yhere, col=colorz[1], lty=lty.here[1],
    xlab=xname, ylab=yname, main=paste("run ", runstouse[whichrun]))
lines(modelruns[[whichrun]]$Bfin[,1]~xheretauI, ylim=yhere, col=colorz[2], lty=lty.here[2],
    xlab=xname, ylab=yname)
lines(modelruns[[whichrun]]$Bfin[,2]~xheretauI, ylim=yhere, col=colorz[3], lty=lty.here[3],
    xlab=xname, ylab=yname)
legend("topright", leg.txt, lwd=2, col=colorz)

yhere <- c(0, 1)
yname <- "tauP or tauIhat"
leg.txt <- c("tauP", "tauIhat 1", "tauIhat 2")
lty.here <- c(rep(1, 1), rep(2, 2))

plot(modelruns[[whichrun]]$envtvars[["tauP"]]~xhere, type="l", ylim=yhere, col=colorz[1], lty=lty.here[1],
    xlab=xname, ylab=yname)
lines(modelruns[[whichrun]]$tauIhat[,1]~xhere, ylim=yhere, col=colorz[2], lty=lty.here[2],
    xlab=xname, ylab=yname)
lines(modelruns[[whichrun]]$tauIhat[,2]~xhere, ylim=yhere, col=colorz[3], lty=lty.here[3],
    xlab=xname, ylab=yname)
legend("topright", leg.txt, lwd=2, lty=lty.here, col=colorz)
}

dev.off()

## Step 2: Look at the differences in parameters between the two species
# and color dots by coexisting yay or nay

df.coexist <- subset(df, coexist==2)

plot.params <- function(df, df.coexist, colname, colname.sp1, colname.sp2){
    pdf(paste("graphs/modelruns/params/Track_varR_2spp_", jobID, colname, "_compare.pdf", sep=""),
        width=9, height=4)
    plot.allruns <- ggplot(data=df, aes(x=df[colname.sp1], y=df[colname.sp2],
        colour=coexist)) +
        geom_point() +
        labs(x = colname.sp1, y=colname.sp2)
    plot.coexistruns <- ggplot(data=df.coexist, aes(x=df.coexist[colname.sp1], 
        y=df.coexist[colname.sp2], colour=coexist)) +
        geom_point() +
        labs(x = colname.sp1, y=colname.sp2)
    multiplot(plot.allruns, plot.coexistruns, cols=2)
    dev.off()
}

plot.params(df, df.coexist, "alpha", "alpha.sp1", "alpha.sp2")
plot.params(df, df.coexist, "tauI", "tauI.sp1", "tauI.sp2")
plot.params(df, df.coexist, "tauIPini", "tauIPini.sp1", "tauIPini.sp2")
plot.params(df, df.coexist, "tauIhat", "tauIhat.sp1", "tauIhat.sp2")
plot.params(df, df.coexist, "c", "c.sp1", "c.sp2")
plot.params(df, df.coexist, "rstar", "rstar.sp1", "rstar.sp2")


## Step 3: g_i ~ eff|tauI-tauP|
# for a couple coexisting runs and a couple not coexisting

gi.runstouse <- as.numeric(row.names(subset(as.data.frame(whocoexisted), whocoexisted>1))[1:3])
gi.runstouse.nocoexist <- as.numeric(row.names(subset(as.data.frame(whocoexisted),
    whocoexisted==1))[1:3])

xlim <- c(0,1)

pdf(paste("graphs/modelruns/gi_dynamics/Track_varR_2spp_", jobID, "_gi_3coexistruns.pdf", sep=""), width=5, height=7)
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

pdf(paste("graphs/modelruns/gi_dynamics/Track_varR_2spp_", jobID, "_gi_3nocoexistruns.pdf", sep=""), width=5, height=7)
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

## Step 4: Look at parameter differences between species

plot.paramdiffs <- function(df, df.coexist, figname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/Track_varR_2spp_", jobID, figname, ".pdf", sep=""),
        width=9, height=4)
    plot.allruns <- ggplot(data=df, aes(x=df[colname.x], y=df[colname.y],
        colour=coexist)) +
        geom_point() +
        labs(x = colname.x, y=colname.y)
    plot.coexistruns <- ggplot(data=df.coexist, aes(x=df.coexist[colname.x], 
        y=df.coexist[colname.y], colour=coexist)) +
        geom_point() +
        labs(x = colname.x, y=colname.y)
    multiplot(plot.allruns, plot.coexistruns, cols=2)
    dev.off()
}

plot.paramdiffs(df, df.coexist, "rstar_vs_tauI", "diff.tauI", "diff.rstar")
plot.paramdiffs(df, df.coexist, "rstar_vs_alpha", "diff.alpha", "diff.rstar")
plot.paramdiffs(df, df.coexist, "rstar_vs_tauIhat", "diff.tauIhat", "diff.rstar")
plot.paramdiffs(df, df.coexist, "rstar_vs_tauIPini", "diff.tauIPini", "diff.rstar")
plot.paramdiffs(df, df.coexist, "gi_vs_tauI", "diff.tauI", "diff.gi")
plot.paramdiffs(df, df.coexist, "gi_vs_tauIPini", "diff.tauIPini", "diff.gi")
plot.paramdiffs(df, df.coexist, "gi_vs_alpha", "diff.alpha", "diff.gi")
plot.paramdiffs(df, df.coexist, "gi_vs_c", "diff.c", "diff.gi")

plot.paramdiffs(df, df.coexist, "rstar_vs_tauI_ratios", "tauI.ratio", "rstar.ratio")
plot.paramdiffs(df, df.coexist, "rstar_vs_alpha_ratios", "alpha.ratio", "rstar.ratio")
plot.paramdiffs(df, df.coexist, "rstar_vs_tauIhat_ratios", "tauIhat.ratio", "rstar.ratio")
plot.paramdiffs(df, df.coexist, "rstar_vs_tauIPini_ratios", "tauIPini.ratio", "rstar.ratio")


####
####
####

stop(print("stop here! Below is just list help notes ..."))

## Wait, first a plot of the alpha-rstar result ...
plot(diff.rstar~diff.alpha, data=df.coexist, pch=16, col="deeppink4",
     ylab="R* difference between 2 species",
     xlab="Tracking difference between 2 species")

plot(diff.rstar~diff.tauIPini, data=df.coexist, pch=16, col="deepskyblue4",
     ylab="R* difference between 2 species",
     xlab="Effective tracking difference between 2 species")

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
