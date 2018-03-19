## Code taken from PhenologyModelAnalysisOLD.R ##
## To remind me of which figures I still need to make ... ##

#######################
# Have started editing this #
#######################

## Step 1: check out a couple coexisting runs
runbfinhere <- runs1bfin
df.coexist <- df.coexist
bfincoexist <- runs1bfin[which(runbfinhere$taskrunID %in% df.coexist$taskrunID),]

pdf(paste("graphs/modelruns/Track_varR_2spp_", folderID, "_sampleruns.pdf", sep=""), width=5, height=7)
par(mfrow=c(3,2))

for (whichrunnum in c(1:length(df.coexist))){
whichrun <- runstouse[whichrunnum]
yhere <- c(-6, 3)
colorz <- c("black", "red", "blue")

## START HERE  ... 
xhere <- c(1:length(modelruns[[1]]$envtvars[["tauP"]]))
xheretauI <- c(1:length(modelruns[[1]]$Bfin[,1]))
xname <- "time"
yname <- "tauP or Bfin"
leg.txt <- c("tauP", "Bfin sp1", "Bfin sp2")
lty.here <- c(rep(1, 3))

plot(modelruns[[whichrun]]$envtvars[["tauP"]]~xhere, type="l", ylim=yhere, col=colorz[1], lty=lty.here[1],
    xlab=xname, ylab=yname, main=paste("run ", whichrun))
lines(log10(modelruns[[whichrun]]$Bfin[,1])~xheretauI, ylim=yhere, col=colorz[2], lty=lty.here[2],
    xlab=xname, ylab=yname)
lines(log10(modelruns[[whichrun]]$Bfin[,2])~xheretauI, ylim=yhere, col=colorz[3], lty=lty.here[3],
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


## Step 3: g_i ~ eff|tauI-tauP|
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
