### Started 22 August 2016 ###
### By Lizzie ###

## Analyzing the runs data for the storage effect model! ##
## Last updated 8 January 2017 (from Oahu) ##

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

# source some stuff
source("sourcefiles/multiplot.R")

jobID <- "78476247"
load(paste("..//ModelRuns/Track_varR_2spp_", jobID, ".Rdata", sep=""))

# Here's the structure:
# list(jobID=jobID, arrayNum=a, runNum=r, sppvars=sppvars (the unchanging ones and means of the changing ones),
#    envtvars=envtvars (within which is R0, tauP, eps),
#    tauIhat, Bfin=Bfin, g,
#    Bout=Bout (within Bout is time, R, B1, B2),
# head(modelruns[[1]]) # err, not always working so check the top of what you get with:
# str(modelruns)

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


# make a dataframe with necessary info for step 2-4

tauI.coexist.df <- data.frame(modelrun=rep(0, length(modelruns)),
    tauI.sp1=rep(0, length(modelruns)), tauI.sp2=rep(0, length(modelruns)),
    alpha.sp1=rep(0, length(modelruns)), alpha.sp2=rep(0, length(modelruns)), 
    tauIhat.sp1=rep(0, length(modelruns)), tauIhat.sp2=rep(0, length(modelruns)),
    gi.sp1= rep(0, length(modelruns)), gi.sp2= rep(0, length(modelruns)),
    c.sp1=rep(0, length(modelruns)), c.sp2=rep(0, length(modelruns)),
    rstar.sp1=rep(0, length(modelruns)), rstar.sp2=rep(0, length(modelruns)),
    coexist=whocoexisted)

for (k in c(1:length(modelruns))){
    tauI.coexist.df$modelrun[k] <- k
    tauI.coexist.df$tauI.sp1[k] <- modelruns[[k]]$sppvars[["tauIPini"]][1]
    tauI.coexist.df$tauI.sp2[k] <- modelruns[[k]]$sppvars[["tauIPini"]][2]
    tauI.coexist.df$alpha.sp1[k] <- modelruns[[k]]$sppvars[["alpha"]][1]
    tauI.coexist.df$alpha.sp2[k] <- modelruns[[k]]$sppvars[["alpha"]][2]
    tauI.coexist.df$tauIhat.sp1[k] <- mean(modelruns[[k]]$tauIhat[,1])
    tauI.coexist.df$tauIhat.sp2[k] <- mean(modelruns[[k]]$tauIhat[,2])
    tauI.coexist.df$gi.sp1[k] <- mean(modelruns[[k]]$g[,1])
    tauI.coexist.df$gi.sp2[k] <- mean(modelruns[[k]]$g[,2])
    tauI.coexist.df$c.sp1[k] <- modelruns[[k]]$sppvars[["c"]][1]
    tauI.coexist.df$c.sp2[k] <- modelruns[[k]]$sppvars[["c"]][2]
    tauI.coexist.df$rstar.sp1[k] <- modelruns[[k]]$sppvars[["Rstar"]][1]
    tauI.coexist.df$rstar.sp2[k] <- modelruns[[k]]$sppvars[["Rstar"]][2]
    }

tauI.coexist.df$coexist <- as.factor(tauI.coexist.df$coexist)
tauI.coexist.df$diff.rstar <- tauI.coexist.df$rstar.sp1-tauI.coexist.df$rstar.sp2
tauI.coexist.df$diff.tauIhat <- tauI.coexist.df$tauIhat.sp1-tauI.coexist.df$tauIhat.sp2
tauI.coexist.df$diff.alpha <- tauI.coexist.df$alpha.sp1-tauI.coexist.df$alpha.sp2
tauI.coexist.df$diff.gi <- tauI.coexist.df$gi.sp1-tauI.coexist.df$gi.sp2
tauI.coexist.df$diff.c <- tauI.coexist.df$gi.sp1-tauI.coexist.df$gi.sp2




###############################################
## some plots ##
## (working off 20160420_2sppCoexistFigs.pdf)
###############################################

## Step 1: check out a couple runs
runstouse <- c(1, 110, 851)

pdf(paste("graphs/modelruns/Track_varR_2spp_", jobID, "_3runs.pdf", sep=""), width=5, height=7)
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
    xlab=xname, ylab=yname, main=paste("run ", whichrun))
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

## Step 2: tauI_sp1 vs tauI_sp2 and color dots by coexisting yay or nay
# and while I was at it, I added some rstar comparisons

# here comes your plot!
tauI3 <- ggplot(data=tauI.coexist.df, aes(x=tauI.sp1, y=tauI.sp2, colour=coexist)) +
    geom_point()

tauI1 <- ggplot(data=subset(tauI.coexist.df, coexist==2), aes(x=tauI.sp1, y=tauI.sp2)) +
    geom_point()

tauIhat3 <- ggplot(data=tauI.coexist.df, aes(x=tauIhat.sp1, y=tauIhat.sp2, colour=coexist)) +
    geom_point()

tauIhat1 <- ggplot(data=subset(tauI.coexist.df, coexist==2), aes(x=tauIhat.sp1, y=tauIhat.sp2)) +
    geom_point()

rstar3 <- ggplot(data=tauI.coexist.df, aes(x=rstar.sp1, y=rstar.sp2, colour=coexist)) +
    geom_point()

rstar1 <- ggplot(data=subset(tauI.coexist.df, coexist==2), aes(x=rstar.sp1, y=rstar.sp2)) +
    geom_point()

alpha1 <- ggplot(data=tauI.coexist.df, aes(x=alpha.sp1, y=alpha.sp2, colour=coexist)) +
    geom_point()

alpha2 <- ggplot(data=subset(tauI.coexist.df, coexist==2), aes(x=alpha.sp1, y=alpha.sp2)) +
    geom_point()


# called Track_varR_2spp_78476247_coexist1.pdf now, but need better name
quartz()
multiplot(tauI3, tauIhat3, rstar3, alpha3, tauI1, tauIhat1, rstar1, alpha1, cols=2)

## Step 3: g_i ~ eff|tauI-tauP|
## ADD histogram of tauP as underlay as side by side plot

# try to pick a run that coexisted!
quartz()
plot(modelruns[[1]]$g[,1]~modelruns[[1]]$tauIhat[,1])
points(modelruns[[1]]$g[,2]~modelruns[[1]]$tauIhat[,2], col="blue")

quartz()
plot(modelruns[[3]]$g[,1]~modelruns[[3]]$tauIhat[,1])
points(modelruns[[3]]$g[,2]~modelruns[[3]]$tauIhat[,2], col="blue")

## Step 4: diff in Rstar between spp versus diff in tauI between spp

rstardiff <- ggplot(data=tauI.coexist.df, aes(x=diff.rstar, y=diff.tauIhat, colour=coexist)) +
    geom_point()

tauIhatdiff <- ggplot(data=subset(tauI.coexist.df, coexist==2),
    aes(x=diff.rstar, y=diff.tauIhat, colour=coexist)) +
    geom_point()

# called Track_varR_2spp_78476247_coexist2.pdf now, but need better name
quartz()
multiplot(rstardiff, tauIhatdiff, cols=2)



## Few new ones:

# mean gi versus alpha
 ggplot(data=tauI.coexist.df, aes(x=alpha.sp1, y=gi.sp1, colour=coexist)) +
    geom_point()

# diff gi versus diff rstar
ggplot(data=tauI.coexist.df, aes(x=diff.gi, y=diff.rstar, colour=coexist)) +
    geom_point()
ggplot(data=subset(tauI.coexist.df, coexist==2), aes(x=diff.gi, y=diff.rstar, colour=coexist)) +
    geom_point()

# diff gi versus diff c
ggplot(data=subset(tauI.coexist.df, coexist==2), aes(x=diff.gi, y=diff.c, colour=coexist)) +
    geom_point()

####
####
####

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
