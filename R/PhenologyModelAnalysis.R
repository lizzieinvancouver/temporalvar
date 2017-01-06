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

# source some stuff
source(multiplot)

## set working directory
setwd("~/Documents/git/projects/temporalvar/R")

jobID <- "78345799"
load(paste("..//notposting/workon/Track_varR_2spp_", jobID, ".Rdata", sep=""))

# Here's the structure:
# list(jobID=jobID, arrayNum=a, runNum=r, sppvars=sppvars (the unchanging ones and means of the changing ones),
#    tauI=tauI, envtvars=envtvars (within which is R0, tauP, eps),
#    Bfin=Bfin, Bout=Bout (within Bout is time, R, B1, B2) )

# checking out Rstar ratios of the two species model
rstarz <- c()
for (i in c(1:length(modelruns))){
   rstarz[i] <- modelruns[[i]]$sppvars[["Rstar"]][1]/modelruns[[i]]$sppvars[["Rstar"]][2]
}
hist(rstarz, breaks=100)

# checking out which runs coexisted
# note: something is up with the Bfin so that it's always 0,0 at the last timestep (200)
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

## some plots (working off 20160420_2sppCoexistFigs.pdf)

## Step 1: check out a couple runs
runstouse <- c(1, 110, 851)

pdf(paste("graphs/Track_varR_2spp_", jobID, "_3runs.pdf", sep=""), width=5, height=7)
par(mfrow=c(3,2))

for (whichrun in seq_along(runstouse)){
yhere <- c(0, 65)
colorz <- c("black", "red", "blue")

xhere <- c(1:length(modelruns[[1]]$envtvars[["tauP"]]))
xname <- "time"
yname <- "tauP or Bfin"
leg.txt <- c("tauP", "Bfin sp1", "Bfin sp2")
lty.here <- c(rep(1, 3))

plot(modelruns[[whichrun]]$envtvars[["tauP"]]~xhere, type="l", ylim=yhere, col=colorz[1], lty=lty.here[1],
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

plot(modelruns[[whichrun]]$envtvars[["tauP"]]~xhere, type="l", ylim=yhere, col=colorz[1], lty=lty.here[1],
    xlab=xname, ylab=yname)
lines(modelruns[[whichrun]]$tauI[,1]~xhere, ylim=yhere, col=colorz[2], lty=lty.here[2],
    xlab=xname, ylab=yname)
lines(modelruns[[whichrun]]$tauI[,2]~xhere, ylim=yhere, col=colorz[3], lty=lty.here[3],
    xlab=xname, ylab=yname)
legend("topright", leg.txt, lwd=2, lty=lty.here, col=colorz)
}

dev.off()

## Step 2: tauI_sp1 vs tauI_sp2 and color dots by coexisting yay or nay
# and while I was at it, I added some rstar comparisons

# make a dataframe with necessary info

## START HERE! Need to pull out effective tauI
tauI.coexist.df <- data.frame(modelrun=rep(0, length(modelruns)),
    tauI.sp1=rep(0, length(modelruns)), tauI.sp2=rep(0, length(modelruns)),
    # tauIeff.sp1=rep(0, length(modelruns)), tauIeff.sp2=rep(0, length(modelruns)),
    rstar.sp1=rep(0, length(modelruns)), rstar.sp2=rep(0, length(modelruns)),
    coexist=whocoexisted)

for (k in c(1:length(modelruns))){
    tauI.coexist.df$modelrun[k] <- k
    tauI.coexist.df$tauI.sp1[k] <- modelruns[[k]]$sppvars[["tauIPini"]][1]
    tauI.coexist.df$tauI.sp2[k] <- modelruns[[k]]$sppvars[["tauIPini"]][2]
    ## ADD ME
    tauI.coexist.df$rstar.sp1[k] <- modelruns[[k]]$sppvars[["Rstar"]][1]
    tauI.coexist.df$rstar.sp2[k] <- modelruns[[k]]$sppvars[["Rstar"]][2]
    }
tauI.coexist.df$coexist <- as.factor(tauI.coexist.df$coexist)

# here comes your plot!
tauI3 <- ggplot(data=tauI.coexist.df, aes(x=tauI.sp1, y=tauI.sp2, colour=coexist)) +
    geom_point()

tauI1 <- ggplot(data=subset(tauI.coexist.df, coexist==2), aes(x=tauI.sp1, y=tauI.sp2)) +
    geom_point()

rstar3 <- ggplot(data=tauI.coexist.df, aes(x=rstar.sp1, y=rstar.sp2, colour=coexist)) +
    geom_point()

rstar1 <- ggplot(data=subset(tauI.coexist.df, coexist==2), aes(x=rstar.sp1, y=rstar.sp2)) +
    geom_point()

multiplot(coexist3, rstar3, coexist1, rstar1, cols=2)

# old code .... delete?
whichrun.now <- 1
plot(modelruns[[whichrun.now]]$tauI[,1]~modelruns[[whichrun.now]]$tauI[,2])
whichrun.now <- 111
plot(modelruns[[whichrun.now]]$tauI[,1]~modelruns[[whichrun.now]]$tauI[,2])


## Step 3: g_i ~ eff|tauI-tauP|
# hmmm, we need g_i for this plot, which I don't have


## Step 4: diff in Rstar between spp versus diff in tauI between spp


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
