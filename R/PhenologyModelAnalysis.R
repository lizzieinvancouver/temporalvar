### Started 13 February 2018 ###
### By Lizzie ###

## Analyzing the runs data for the storage effect model! ##
## Last updated August 2018 ##

######################
### To do items!!! ###
######################
# https://www.theanalysisfactor.com/r-tutorial-part-12/
# (1) Make sure the declining R0 runs with BIGGER declines: 8995922, 8995924 are being handled correctly given the way I use seq(XX, length(runz), 5)
# (2) Double check how I handle the ext.alpha runs
######################

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## libraries
library(ggplot2)
require(plyr); require(dplyr); require(tidyr) # data formatting

## set working directory
setwd("~/Documents/git/projects/temporalvar/R")

source("sourcefiles/analyses/multiplot.R") # used in plot.params
source("sourcefiles/analyses/runanalysisfxs.R")

writeout.someruns.formegan <- FALSE
lookdeeplyatR0runs <- FALSE
inclextalpharuns <- TRUE

runshaveheader <- TRUE

# cheap loop over the files for now
runz <- c("858179", "858221", "858241", "858262", "858282",
          "888288", "888338", "888369", "888380", "888430",
          "888600", "888602", "888605", "888607", "888608",
          "933059", "933107", "933156", "933215", "933272",
          "933566", "933600", "933630", "933682", "933723",
          "1006570", "1006591", "1006608", "1006628", "1006639",
          "8995922", "8995924")
# declining R0 runs with BIGGER declines: 8995922, 8995924 (see alphaRstarR0ext below)
extalpharuns <- c("12519313", "12519314") # extended alpha runs are 12519313, 12519314: first set fills in the alpha - 0-0.3 range; second set is for the full range of alpha (0-1)

# Remember to update below under plotting-related formatting ...
seq(1, length(runz), 5) # 1, are varying everything
seq(2, length(runz), 5) # 2, are NOT varying tauI
seq(3, length(runz), 5) # 3, are NOT varying tracking
seq(4, length(runz), 5) # 4, keeps R* the same across species pairs
seq(5, length(runz), 5) # 5, vary alpha, Rstar, and add in a R0 declines
# (with tauP still getting earlier, as in all other runs)

# runinfo <- read.table("Table_of_RunParms.txt", skip=1, header=TRUE)

# varying (as of 30 August 2018)
# tauI, Rstar <- 440/33148 (1%)
# alpha, Rstar <- 1237/41345 (3%)
# tauI, alpha <- 1166/38808 (3%)
# all varies <-   1214/39014 coexist (3%)

##############################
## Notes on other runs .... ##
##############################


#########################################
## Do some data reading and formatting ##
#########################################
runnow <- c(runz)
if(inclextalpharuns){ # overwrites above if flagged
runnow <- c(runz, extalpharuns)
}

## Setup for pasting runs together into one df (not sure we permanently want this but I want it now)
folderID <- runnow[1]
samplerun <-  read.table(paste("output/SummaryFiles/",folderID,"/SummaryOut_", folderID,
    "-1.txt", sep=""), header=TRUE) # comment.char = "", 
df.all <- data.frame(matrix(ncol=length(colnames(samplerun)), nrow=0))
colnames(df.all) <- colnames(samplerun)

for(folderIDhere in c(1:length(runnow))){
    
folderID <- runnow[folderIDhere] # folderID <- 41801534
if(runshaveheader){
    # Do we need the first line here?
    samplerun <-  read.table(paste("output/SummaryFiles/", folderID, "/SummaryOut_", folderID,
        "-1.txt", sep=""),  header=TRUE)
    file.names <- dir(paste("output/SummaryFiles/", folderID, sep=""), pattern =".txt")
    colnameshere <- colnames(samplerun)
    runs1 <- getfiles(folderID, file.names, colnameshere)
    runs1$taskrunID <- paste(runs1$taskID, runs1$runID, sep="-")
}

if(!runshaveheader){
    file.names <- dir(paste("output/SummaryFiles/", folderID, sep=""), pattern =".txt")
    colnameshere <- colnames(samplerun)
    runs1 <- getfiles(folderID, file.names, colnameshere)
    runs1$taskrunID <- paste(runs1$taskID, runs1$runID, sep="-")
}
 
 
##
## Data formatting to compare species pairs
##
df <- makediffs(runs1)
df <- calcbesttauI(df)
df <- calcsp.besttauI(df)
df <- calcsp.bestrstar(df)
df <- calcsp.bestalpha(df)
df <- calcsp.biggerslopeBfin(df)
df.coexist <- subset(df, ncoexist==2)
print(paste("the current folder ID is", folderID, "the total rows are:", nrow(df), 
    "and the total coexisting rows are:", nrow(df.coexist), sep=" "))

# altogether now! pasting runs together into one df
df.all <- rbind(df.all, df)

##############################
## Set tauP for the graphs! ##
##############################
source("sourcefiles/analyses/tauP.R")


#################################
## Plotting-related formatting ##
#################################

##
## Data formatting to compare species pairs in plots ...
## Build df of just coexisting in period 1 (stat) and keep in wide format 
##
df.coexist1 <- subset(df, ncoexist==2 & period==1)
df.t2 <- subset(df, period==2)
df.t2 <- subset(df.t2, select=c("jobID", "taskID", "runID", "ncoexist", "taskrunID"))
df.plot <- merge(df.coexist1, df.t2, by=c("jobID", "taskID", "runID", "taskrunID"),
    all.x=TRUE, all.y=FALSE, suffixes=c(".t1", ".t2"))

# do again for df.all
# The things that vary by stat/non-stat: R0 (R0_mean, R0_median, R0_autocor), R, tauP, tauIP, Bfin (and related), the gmeans
df.all.coexist1 <- subset(df.all, ncoexist==2 & period==1)
df.all.t2 <- subset(df.all, period==2)
df.all.t2 <- subset(df.all.t2, select=c("jobID", "taskID", "runID", "ncoexist",
    "coexist1", "coexist2", "taskrunID", "R0_mean", "R0_median", "ratio.g", "ratio.tauIP",
    "ratio.tauIPnoalpha", "tauIP1_mean", "tauIP2_mean", "diff.bfinslopes", "slopeBfin1", "slopeBfin2", 
    "minslopeBfin"))
df.all.plot <- merge(df.all.coexist1, df.all.t2, by=c("jobID", "taskID", "runID", "taskrunID"),
    all.x=TRUE, all.y=FALSE, suffixes=c(".t1", ".t2"))

# alternative of above, for when you want the df structure of above but with the ncoexist=(1 or 2) the runs from the stationary (df.all.plot is just ncoexist=2 from the stat period)
df.all.coexist1.justonesp <- subset(df.all, ncoexist==1 & period==1)
df.all.coexist.alt <- rbind(df.all.coexist1, df.all.coexist1.justonesp)
df.all.plot.alt <- merge(df.all.coexist.alt, df.all.t2, by=c("jobID", "taskID", "runID", "taskrunID"),
    all.x=TRUE, all.y=FALSE, suffixes=c(".t1", ".t2"))


##
## Data formatting to compare species 
## Get runs with coexist=2 in period 1 (stat), keep data only in period 2 (ns)
##
df.t2.wp <- subset(df, period==2)
# get just the period 2 runs that had coexistence in period 1
df.long <- df.t2.wp[which(df.t2.wp$taskrunID %in% unique(df.coexist1$taskrunID)),]
df.long.exist <- subset(df.long, coexist1==1 | coexist2==1)
df.long.noexist <- subset(df.long, coexist1==0 | coexist2==0)

# do again for df.all
df.all.t2.wp <- subset(df.all, period==2)
# get just the period 2 runs that had coexistence in period 1
df.all.long <- df.all.t2.wp[which(df.all.t2.wp$taskrunID %in% unique(df.all.coexist1$taskrunID)),]
df.all.long.exist <- subset(df.all.long, coexist1==1 | coexist2==1)
df.all.long.noexist <- subset(df.all.long, coexist1==0 | coexist2==0)
}




##
## Check that when tauI=tauP that species wins ... 
##
# tauRstar.check <-  df.all[which(df.all$jobID %in% runz[c(3,7,11,14)]),]
tauRstar.check <-  df.all[which(df.all$jobID %in% runz[c(3)]),] # tauI and Rstar tradeoff runs
check.sp1 <- subset(tauRstar.check, tauIP1_mean<0.1)
sum(check.sp1$coexist1) # 2776
sum(check.sp1$coexist2) # 657
check.sp2 <- subset(tauRstar.check, tauIP2_mean<0.1)
sum(check.sp2$coexist1) # 632
sum(check.sp2$coexist2) # 2770

try <- subset(check.sp1, coexist2==1)
sum(try$coexist1) # 185, so most of these sp2 wins ... it should have a very low Rstar then compared to sp1.... START HERE and check more thoroughly! Also, check the other species' tauI (i.e., could also be that BOTH species have tauI close to tauP)
hist(try$ratio.rstar)
hist(check.sp1$ratio.rstar)

##
## Group the runs by what type they are so I can plot
##
# See above when I set up runz for where I outline what numbers to pull for each!
tauRstar.runs <- runz[seq(3, length(runz), 5)] # NOT varying tracking: tauI and Rstar tradeoff
alphaRstar.runs <- runz[seq(2, length(runz), 5)] # NOT varying tauI: tracking and Rstar tradeoff
taualpha.runs <- runz[seq(4, length(runz), 5)] # keeps R* the same across species pairs
taualphaRstar.runs <- runz[seq(1, length(runz), 5)] # varying everything (tauI, alpha, Rstar)
alphaRstarR0.runs <- runz[seq(5, length(runz), 5)] # NOT varying tauI: tracking and Rstar tradeoff, add in declining R0!
alphaRstarR0ext.runs <- c("8995922", "8995924")

## Grab just the stationary period
df.all.stat <- subset(df.all, period==1)
# NOT varying tracking
tauRstar.stat.runs.df <- df.all.stat[which(df.all.stat$jobID %in% tauRstar.runs),]
sum(tauRstar.stat.runs.df$diff.alpha) # must equal zero!
# NOT varying tauI
alphaRstar.stat.runs.df <- df.all.stat[which(df.all.stat$jobID %in% alphaRstar.runs),]
sum(alphaRstar.stat.runs.df$diff.tauI) # must equal zero!
# keeps R* the same across species pairs
taualpha.stat.runs.df <- df.all.stat[which(df.all.stat$jobID %in% taualpha.runs),]
sum(taualpha.stat.runs.df$diff.Rstar) # must equal zero!
# varying everything (tauI, alpha, Rstar)
taualphaRstar.stat.runs.df <- df.all.stat[which(df.all.stat$jobID %in% taualphaRstar.runs),]
# add in declining R0
alphaRstarR0.stat.runs.df <- df.all.stat[which(df.all.stat$jobID %in% alphaRstarR0.runs),]
alphaRstarR0ext.stat.runs.df <- df.all.stat[which(df.all.stat$jobID %in% alphaRstarR0ext.runs),]
sum(alphaRstarR0.stat.runs.df$diff.tauI) # must equal zero!
# extended alpha runs
alphaextRstar.stat.runs.df <- df.all.stat[which(df.all.stat$jobID %in% extalpharuns),]



## Grab the stationary coexisting runs and the non-stationary period
# NOT varying tracking
tauRstar.runs.df <- df.all.plot[which(df.all.plot$jobID %in% tauRstar.runs),]
sum(tauRstar.runs.df$diff.alpha) # must equal zero!
# NOT varying tauI
alphaRstar.runs.df <- df.all.plot[which(df.all.plot$jobID %in% alphaRstar.runs),]
sum(alphaRstar.runs.df$diff.tauI) # must equal zero!
# keeps R* the same across species pairs
taualpha.runs.df <- df.all.plot[which(df.all.plot$jobID %in% taualpha.runs),]
sum(taualpha.runs.df$diff.Rstar) # must equal zero!
# varying everything (tauI, alpha, Rstar)
taualphaRstar.runs.df <- df.all.plot[which(df.all.plot$jobID %in% taualphaRstar.runs),]
# add in declining R0
alphaRstarR0.runs.df <- df.all.plot[which(df.all.plot$jobID %in% alphaRstarR0.runs),]
alphaRstarR0ext.runs.df <- df.all.plot[which(df.all.plot$jobID %in% alphaRstarR0ext.runs),]
# extended alpha runs
alphaextRstar.runs.df <- df.all.plot[which(df.all.plot$jobID %in% extalpharuns),]


## Do some counts (these numbers are in manuscript)
# How many runs?
nrow(tauRstar.stat.runs.df)
nrow(alphaRstar.stat.runs.df)
nrow(alphaRstarR0.stat.runs.df)
# How many spp left at of stat period 
sum(tauRstar.stat.runs.df$ncoexist)/(nrow(tauRstar.stat.runs.df)*2)
sum(alphaRstar.stat.runs.df$ncoexist)/(nrow(alphaRstar.stat.runs.df)*2)
sum(alphaRstarR0.stat.runs.df$ncoexist)/(nrow(alphaRstarR0.stat.runs.df)*2)
# How many coexist of non-stat period compared to stat period
nrow(subset(tauRstar.runs.df, ncoexist.t2==2))/nrow(tauRstar.runs.df)
nrow(subset(alphaRstar.runs.df, ncoexist.t2==2))/nrow(alphaRstar.runs.df)
nrow(subset(alphaRstarR0.runs.df, ncoexist.t2==2))/nrow(alphaRstarR0.runs.df)

# How many species left after non-stat? First: (# spp left at end of non-stat)/(# of spp left at end of stat)
sum(tauRstar.runs.df$ncoexist.t2)/sum(tauRstar.stat.runs.df$ncoexist) 
sum(alphaRstar.runs.df$ncoexist.t2)/sum(alphaRstar.stat.runs.df$ncoexist)
sum(taualpha.runs.df$ncoexist.t2)/sum(taualpha.stat.runs.df$ncoexist)
sum(taualphaRstar.runs.df$ncoexist.t2)/sum(taualphaRstar.stat.runs.df$ncoexist)
sum(alphaRstarR0.runs.df$ncoexist.t2)/sum(alphaRstarR0.stat.runs.df$ncoexist)
sum(alphaRstarR0ext.runs.df$ncoexist.t2)/sum(alphaRstarR0ext.stat.runs.df$ncoexist)
# Second: (# spp left at end of non-stat)/(# of co-existing spp left at end of stat)
sum(tauRstar.runs.df$ncoexist.t2)/sum(tauRstar.runs.df$ncoexist.t1) 
sum(alphaRstar.runs.df$ncoexist.t2)/sum(alphaRstar.runs.df$ncoexist.t1)
sum(taualpha.runs.df$ncoexist.t2)/sum(taualpha.runs.df$ncoexist.t1)
sum(taualphaRstar.runs.df$ncoexist.t2)/sum(taualphaRstar.runs.df$ncoexist.t1)
sum(alphaRstarR0.runs.df$ncoexist.t2)/sum(alphaRstarR0.runs.df$ncoexist.t1)
sum(alphaRstarR0ext.runs.df$ncoexist.t2)/sum(alphaRstarR0ext.runs.df$ncoexist.t1)
# From above ... do you see more extinctions after/before non-stat with declining R0?
sum(alphaextRstar.runs.df$ncoexist.t2)/sum(alphaextRstar.runs.df$ncoexist.t1)

# What is average alpha before and after stat?
get.mean.alphavalues(df.all[which(df.all$jobID %in% alphaRstar.runs),])
get.mean.alphavalues.ns(df.all[which(df.all$jobID %in% alphaRstar.runs),])
get.mean.alphavalues(df.all[which(df.all$jobID %in% taualpha.runs),])
get.mean.alphavalues.ns(df.all[which(df.all$jobID %in% taualpha.runs),])
get.mean.alphavalues(df.all[which(df.all$jobID %in% taualphaRstar.runs),])
get.mean.alphavalues.ns(df.all[which(df.all$jobID %in% taualphaRstar.runs),])
get.mean.alphavalues(df.all[which(df.all$jobID %in% alphaRstarR0.runs),])
get.mean.alphavalues.ns(df.all[which(df.all$jobID %in% alphaRstarR0.runs),])
get.mean.alphavalues.ns(df.all[which(df.all$jobID %in% alphaRstarR0ext.runs),])
get.mean.alphavalues.ns(df.all[which(df.all$jobID %in% extalpharuns),])


if(writeout.someruns.formegan){
trstat2spp <- subset(tauRstar.stat.runs.df, ncoexist==2)
trstat2spp.sm <- subset(trstat2spp, select=c("jobID", "taskID", "runID"))
arstat2spp <- subset(alphaRstar.stat.runs.df, ncoexist==2)
arstat2spp.sm <- subset(arstat2spp, select=c("jobID", "taskID", "runID"))
write.csv(trstat2spp.sm, "output/statcoexist.taurstar.csv", row.names=FALSE)
write.csv(arstat2spp.sm, "output/statcoexist.alpharstar.csv", row.names=FALSE)
}

###############
## Plotting! ##
###############
coexist3col <- add.alpha(c("firebrick", "dodgerblue", "seagreen"), alpha=0.4)
coexistmocol<- add.alpha(c("firebrick", "dodgerblue", "seagreen", "yellow", "purple"), alpha=0.4)
tauPcol <- add.alpha(c("yellow", "firebrick"), alpha=0.2)
varhistcol <- add.alpha(c("yellow", "firebrick"), alpha=0.8)
# col2rgb helps here ...
leg.txt <- c("poof", "1 left", "2 survive")
cexhere=0.6
pchhere=16

# For the bfinslope plots
if(FALSE){
library(RColorBrewer)
cols = brewer.pal(4, "RdBu")
# Define colour pallete
pal = colorRampPalette(c("blue", "red"))
# Use the following line with RColorBrewer
colpalettehere = colorRampPalette(cols)
}

library(viridis)
# https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
colpalettehere=viridis


### histograms old code, remove?
if(FALSE){
breaknum <- 20
hist(c(df.long.exist$tauI1, df.long.exist$tau2), xlim=c(0,1), ylim=c(0,30),
     breaks=breaknum, col=coexist3col[3], main="", xlab="number")
par(new=TRUE)
hist(c(df.long.noexist$tauI1, df.long.noexist$tau2), xlim=c(0,1), ylim=c(0,30),
     breaks=breaknum, col=coexist3col[1], main="", xlab="", ylab="") 

plot.histograms(df.long.exist, df.long.noexist, "tauI", "tauI1", "tauI2",
    coexist3col, seq(from=0, to=1, by=0.05), c(0,1), c(0,40))
ggplot(df.plot, aes(x=ratio.tauIP.t2, color=as.factor(ncoexist.t2), fill=as.factor(ncoexist.t2))) +
   geom_histogram(alpha=0.5, position="identity")

ggplot() + geom_density(data=tauP.plot, aes(x=tauP), alpha=0.25) +
    geom_histogram(data=df.plot, aes(x=tauI1, color=as.factor(ncoexist.t2), fill=as.factor(ncoexist.t2)))
}



##
## Histograms and tauI
##
# not varying alpha
plot.histograms.bothspp(tauRstar.runs.df, "tauRstar", "tauI1", "tauI2",
    tauPcol, varhistcol, ylim=c(0,18))
plot.histograms.onespp.skipnonstat(tauRstar.runs.df, "tauRstar", "besttauI",
    tauPcol, varhistcol, ylim=c(0,18)) # doesn't run for nonstat period because there is too little data!
plot.histograms.bars.onespp.skipnonstat(tauRstar.runs.df, "tauRstar", "besttauI",
    tauPcol, varhistcol, 10)
# NOT varying tauI
plot.histograms.bothspp(alphaRstar.runs.df, "alphaRstar", "alpha1", "alpha2",
    tauPcol, varhistcol, ylim=c(0,5))
plot.histograms.onespp(alphaRstar.runs.df, "alphaRstar", "besttauI",
    tauPcol, varhistcol, ylim=c(0,5))
# keeps R* the same across species pairs
plot.histograms.bothspp(taualpha.runs.df, "taualpha", "alpha1", "alpha2",
    tauPcol, varhistcol, ylim=c(0,5))
plot.histograms.onespp(taualpha.runs.df, "taualpha", "besttauI",
    tauPcol, varhistcol, ylim=c(0,5))
# varying everything (tauI, alpha, Rstar)
plot.histograms.bothspp(taualphaRstar.runs.df, "taualphaRstar", "alpha1", "alpha2",
    tauPcol, varhistcol, ylim=c(0,5))
plot.histograms.max(taualphaRstar.runs.df, "taualphaRstar", "alpha1", "alpha2",
    tauPcol, varhistcol, ylim=c(0,5))
plot.histograms.min(taualphaRstar.runs.df, "taualphaRstar", "alpha1", "alpha2",
    tauPcol, varhistcol, ylim=c(0,5))
plot.histograms.onespp(taualphaRstar.runs.df, "taualphaRstar", "besttauI",
    tauPcol, varhistcol, ylim=c(0,5))

# ggplot(tauP.plot, aes(x=tauP, fill=when)) + geom_density(alpha=0.25)


##
## Paramdiff plots
##

### First, some stationary period plots
plot.paramdiffs.stat.onepanel(tauRstar.stat.runs.df, "tauRstar.runs", "_tauIP.rstar", "ratio.tauIP",
    "ratio.rstar", cexhere, pchhere)
plot.paramdiffs.stat.onepanel(alphaRstar.stat.runs.df, "alphaRstar.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere)
plot.paramdiffs.stat.onepanel(taualpha.stat.runs.df, "taualpha.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP", cexhere, pchhere)

plot.rstar.winnersp.stat(tauRstar.stat.runs.df, "tauRstar.runs", "_tauIP.rstar", "ratio.tauIP",
    "ratio.rstar", cexhere, pchhere)
plot.rstar.winnersp.stat(alphaRstar.stat.runs.df, "alphaRstar.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere)
plot.tauI.winnersp.stat(tauRstar.stat.runs.df, "tauRstar.runs", "_tauIP.rstar", "ratio.tauIP",
    "ratio.rstar", cexhere, pchhere)
plot.tauI.winnersp.stat(taualpha.stat.runs.df, "taualpha.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP", cexhere, pchhere)
plot.tauI.winnersp.stat.alt(tauRstar.stat.runs.df, "tauRstar.runs", "_tauIP.rstar", "ratio.tauIP",
    "ratio.rstar", cexhere, pchhere)
plot.tauI.winnersp.stat.alt(taualpha.stat.runs.df, "taualpha.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP", cexhere, pchhere)
plot.alpha.winnersp.stat(alphaRstar.stat.runs.df, "alphaRstar.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere)
plot.alpha.winnersp.stat(taualpha.stat.runs.df, "taualpha.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP", cexhere, pchhere)
plot.alpha.winnersp.stat.alt(alphaRstar.stat.runs.df, "alphaRstar.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere)
plot.alpha.winnersp.stat.alt(taualpha.stat.runs.df, "taualpha.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP", cexhere, pchhere)

plot.alpha.winnersp.stat(alphaextRstar.stat.runs.df, "alphaRextstar.runs", "_alphaext.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere)

# Look at bfin at end of stationary for tauRstar and alphaRstar runs, little extra work
tauRstar.runs.df.alt <- df.all.plot.alt[which(df.all.plot.alt$jobID %in% tauRstar.runs),]
plot.paramdiffs.stat.bfin(tauRstar.runs.df.alt, "tauRstar.runs", "_tauIP.t1.rstar", "ratio.tauIP.t1",
    "ratio.rstar", cexhere, pchhere,  "sp1 wins", "bottomleft", "sp2 wins", "topright", colpalettehere)
alphaRstar.runs.df.alt <- df.all.plot.alt[which(df.all.plot.alt$jobID %in% alphaRstar.runs),]
plot.paramdiffs.stat.bfin(alphaRstar.runs.df.alt, "alphaRstar.runs", "_tauIP.t1.rstar", "ratio.tauIP.t1",
    "ratio.rstar", cexhere, pchhere,  "sp1 wins", "bottomleft", "sp2 wins", "topright", colpalettehere)
plot.paramdiffs.stat.bfin(alphaRstar.runs.df.alt, "alphaRstar.runs", "alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere,  "sp1 wins", "bottomright", "sp2 wins", "topleft", colpalettehere)


### Now, including the non-stationary period (no longer showing ncoexist=0 or 1 from stat period)

# Two things vary ...
plot.paramdiffs.onepanel(tauRstar.runs.df, "tauRstar.runs", "_tauIP.t1.rstar", "ratio.tauIP.t1",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.onepanel(tauRstar.runs.df, "tauRstar.runs", "_tauIP.rstar", "ratio.tauIP.t2",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.twopanel(tauRstar.runs.df, "tauRstar.runs", "_tauIP.rstar", "ratio.tauIP.t2",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.twopanel(tauRstar.runs.df, "tauRstar.runs", "_tauIP.t1.rstar", "ratio.tauIP.t1",
    "ratio.rstar", cexhere, pchhere,  "sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.twopanel(tauRstar.runs.df, "tauRstar.runs", "_tauI.rstar", "ratio.tauI",
    "ratio.rstar", cexhere, pchhere,  "sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.manypanel.bfin(tauRstar.runs.df, "tauRstar.runs", "_tauI.rstar", "ratio.tauI",
    "ratio.rstar", cexhere, pchhere,  "sp1 wins", "bottomleft", "sp2 wins", "topright", colpalettehere)
plot.paramdiffs.manypanel.bfin(tauRstar.runs.df, "tauRstar.runs", "_tauIP.t2.rstar", "ratio.tauIP.t2",
    "ratio.rstar", cexhere, pchhere,  "sp1 wins", "bottomleft", "sp2 wins", "topright", colpalettehere)
plot.paramdiffs.manypanel.bfin(tauRstar.runs.df, "tauRstar.runs", "_tauIP.t1.rstar", "ratio.tauIP.t1",
    "ratio.rstar", cexhere, pchhere,  "sp1 wins", "bottomleft", "sp2 wins", "topright", colpalettehere)
plot.paramdiffs.onesp.bfin(tauRstar.runs.df, "tauRstar.runs", "_tauI.rstar", "ratio.tauI",
    "ratio.rstar", cexhere, pchhere,  "sp1 wins", "bottomleft", "sp2 wins", "topright", colpalettehere)
plot.paramdiffs.onepanel(alphaRstar.runs.df, "alphaRstar.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(alphaRstar.runs.df, "alphaRstar.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(alphaRstar.runs.df, "alphaRstar.runs", "_tauIP.t1.rstar", "ratio.tauIP.t1",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.twopanel.fixedxy(alphaRstar.runs.df, "alphaRstar.runs", "_tauIP.t1.rstar", "ratio.tauIP.t1",
    "ratio.rstar", cexhere, pchhere, c(0.5,2), c(0.5,2), "sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.manypanel.bfin(alphaRstar.runs.df,"alphaRstar.runs", "_tauIP.t1.rstar", "ratio.tauIP.t1",
    "ratio.rstar", cexhere, pchhere,  "sp1 wins", "bottomleft", "sp2 wins", "topright", colpalettehere)
plot.paramdiffs.twopanel(alphaRstar.runs.df, "alphaRstar.runs", "_tauIP.t2.rstar", "ratio.tauIP.t2",
    "ratio.rstar", cexhere, pchhere,"sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.twopanel(alphaRstar.runs.df, "alphaRstar.runs", "_tauIPnoalpha.t1.rstar", "ratio.tauIPnoalpha.t1",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.manypanel.bfin(alphaRstar.runs.df, "alphaRstar.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft", colpalettehere)
plot.paramdiffs.onesp.bfin(alphaRstar.runs.df, "alphaRstar.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft", colpalettehere)


# extended alpha runs
plot.paramdiffs.twopanel(alphaextRstar.runs.df, "alphaextRstar.runs", "_alphaext.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.manypanel.bfin(alphaextRstar.runs.df, "alphaextRstar.runs", "_alphaext.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft", colpalettehere)
plot.paramdiffs.twopanel.fixedxy(alphaextRstar.runs.df, "alphaextRstar.runs", "_tauIP.t1.rstar", "ratio.tauIP.t1",
    "ratio.rstar", cexhere, pchhere, c(0.5,2), c(0.5,2), "sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.manypanel.bfin(alphaextRstar.runs.df, "alphaextRstar.runs", "_tauIP.t1.rstar", "ratio.tauIP.t1",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomleft", "sp2 wins", "topright", colpalettehere)

plot.paramdiffs.onepanel(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP.t1.",
    "ratio.alpha", "ratio.tauIP.t1", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.onepanel(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP.t2", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP.t2", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP.t1.",
    "ratio.alpha", "ratio.tauIP.t1", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel.fixedxy(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP.t2.",
    "ratio.alpha", "ratio.tauIP.t2", cexhere, pchhere, c(0,3), c(0,5),
    "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel.fixedxy(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP.t1", cexhere, pchhere, c(0,3), c(0,5),
    "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.manypanel.bfin(taualpha.runs.df, "taualpha.runs", "_alpha.tauIPnoalpha.t2.",
    "ratio.alpha", "ratio.tauIPnoalpha.t2", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft",
    colpalettehere)
plot.paramdiffs.manypanel.bfin(taualpha.runs.df, "taualpha.runs", "_alpha.tauIPnoalpha.t1.",
    "ratio.alpha", "ratio.tauIPnoalpha.t1", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft",
    colpalettehere)

# some extras for this tricky run ....
taualpha.runs.df$ratio.tauIP.t1t2 <- (taualpha.runs.df$tauIP1_mean.t1-taualpha.runs.df$tauIP1_mean.t2)/
    (taualpha.runs.df$tauIP2_mean.t1-taualpha.runs.df$tauIP2_mean.t2)
plot.paramdiffs.twopanel(taualpha.runs.df, "tauIalphaextras/taualpha.runs", "_alpha.tauIP.t1t2",
    "ratio.alpha", "ratio.tauIP.t1t2", cexhere, pchhere, "sp1 wins", "bottomright",
    "sp2 wins", "topleft")
plot.paramdiffs.twopanel.fixedxy(taualpha.runs.df, "tauIalphaextras/taualpha.runs", "_alpha.tauIP.t1t2",
    "ratio.alpha", "ratio.tauIP.t1t2", cexhere, pchhere, c(0,3), c(-20,20),
    "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualpha.runs.df, "tauIalphaextras/taualpha.runs", "_alpha.tauIPnoalpha.t1.",
    "ratio.alpha", "ratio.tauIPnoalpha.t1", cexhere, pchhere, "sp1 wins", "bottomright",
    "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualpha.runs.df, "tauIalphaextras/taualpha.runs", "_alpha.tauIPnoalpha.t2.",
    "ratio.alpha", "ratio.tauIPnoalpha.t2", cexhere, pchhere, "sp1 wins", "bottomright",
    "sp2 wins", "topleft")
plot.paramdiffs.twopanel.fixedxy(taualpha.runs.df, "tauIalphaextras/taualpha.runs", "_alpha.tauIPnoalpha.t1.",
    "ratio.alpha", "ratio.tauIPnoalpha.t1", cexhere, pchhere, c(-0.5, 4), c(-1, 10),
    "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel.fixedxy(taualpha.runs.df, "tauIalphaextras/taualpha.runs", "_alpha.tauIPnoalpha.t2.",
    "ratio.alpha", "ratio.tauIPnoalpha.t2", cexhere, pchhere, c(-0.75, 4), c(-1, 10),
    "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualpha.runs.df, "tauIalphaextras/taualpha.runs", "_tauIP.t1.g.t1.",
    "ratio.tauIP.t1t2", "ratio.g.t1", cexhere, pchhere, "sp1 wins", "bottomright",
    "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualpha.runs.df, "tauIalphaextras/taualpha.runs", "_tauIP.t1.g.t2.",
    "ratio.tauIP.t1t2", "ratio.g.t2", cexhere, pchhere, "sp1 wins", "bottomright",
    "sp2 wins", "topleft")
plot.paramdiffs.twopanel.fixedxy(taualpha.runs.df, "tauIalphaextras/taualpha.runs", "_tauIP.t1.g.t1.",
    "ratio.tauIP.t1t2", "ratio.g.t1", cexhere, pchhere, c(-20,20), c(0,5),
    "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel.fixedxy(taualpha.runs.df, "tauIalphaextras/taualpha.runs", "_tauIP.t1.g.t2.",
    "ratio.tauIP.t1t2", "ratio.g.t2", cexhere, pchhere, c(-20,20), c(0,5),
    "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualpha.runs.df, "tauIalphaextras/taualpha.runs", "_alpha.g.t1.",
    "ratio.alpha", "ratio.g.t1", cexhere, pchhere, "sp1 wins", "bottomright",
    "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualpha.runs.df, "tauIalphaextras/taualpha.runs", "_alpha.g.t2.",
    "ratio.alpha", "ratio.g.t2", cexhere, pchhere, "sp1 wins", "bottomright",
    "sp2 wins", "topleft")



# Three things vary
plot.paramdiffs.onepanel(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_tauIP.rstar",
    "ratio.tauIP.t2", "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.onepanel(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_tauI.rstar",
    "ratio.tauI", "ratio.rstar", cexhere, pchhere, "??", "bottomleft", "??", "topright")
plot.paramdiffs.onepanel(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_alpha.rstar",
    "ratio.alpha", "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.onepanel(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP.t2", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_tauIP.rstar",
    "ratio.tauIP.t2", "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.twopanel(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_alpha.rstar",
    "ratio.alpha", "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP.t2", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")


# compare these to their equivalents 
plot.paramdiffs.fixedxy(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_tauIP.rstar",
    "ratio.tauIP.t2", "ratio.rstar", tauRstar.runs.df, 1, 16, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.fixedxy(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_tauI.rstar",
    "ratio.tauI", "ratio.rstar", tauRstar.runs.df, 1, 16, "??", "bottomleft", "??", "topright")
plot.paramdiffs.fixedxy(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_alpha.rstar",
    "ratio.alpha", "ratio.rstar", alphaRstar.runs.df, 1, 16, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.fixedxy(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP.t2", taualpha.runs.df, 1, 16, "sp1 wins", "bottomright", "sp2 wins", "topleft")


# and color code by third variable
plot.paramdiffs.colorbyzvar(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_tauIP.rstar",
    "ratio.tauIP.t2", "ratio.rstar", "ratio.alpha", "3 traits vary: survived after stat", 1)
plot.paramdiffs.colorbyzvar(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_alpha.rstar",
    "ratio.alpha", "ratio.rstar", "ratio.tauIP.t2", "3 traits vary: survived after stat", 1)
plot.paramdiffs.colorbyzvar(taualphaRstar.runs.df, "threethingsvary/taualphaRstar.runs", "_alpha.tauIP",
     "ratio.alpha", "ratio.tauIP.t2", "ratio.rstar", "3 traits vary: survived after stat", 1)

# Look at the R0 runs: alphaRstarR0.runs.df
# First we check if R0 is lower than R* much.
# Set up so negative numbers mean resource is lower than species R*
alphaRstarR0.runs.df$RstarR0sp1.t1 <- alphaRstarR0.runs.df$R0_median.t1 - alphaRstarR0.runs.df$Rstar1
alphaRstarR0.runs.df$RstarR0sp2.t1 <- alphaRstarR0.runs.df$R0_median.t1 - alphaRstarR0.runs.df$Rstar2
alphaRstarR0.runs.df$RstarR0sp1.t2 <- alphaRstarR0.runs.df$R0_median.t2 - alphaRstarR0.runs.df$Rstar1 
alphaRstarR0.runs.df$RstarR0sp2.t2 <- alphaRstarR0.runs.df$R0_median.t2 - alphaRstarR0.runs.df$Rstar2
par(mfrow=c(2,2))
hist(alphaRstarR0.runs.df$RstarR0sp1.t1)
hist(alphaRstarR0.runs.df$RstarR0sp2.t1)
hist(alphaRstarR0.runs.df$RstarR0sp1.t2)
hist(alphaRstarR0.runs.df$RstarR0sp2.t2)
# Okay, so it seems like all the species R* values are still okay
# ... and for the extreme R0 declines
alphaRstarR0ext.runs.df$RstarR0sp1.t1 <- alphaRstarR0ext.runs.df$R0_median.t1 - alphaRstarR0ext.runs.df$Rstar1
alphaRstarR0ext.runs.df$RstarR0sp2.t1 <- alphaRstarR0ext.runs.df$R0_median.t1 - alphaRstarR0ext.runs.df$Rstar2
alphaRstarR0ext.runs.df$RstarR0sp1.t2 <- alphaRstarR0ext.runs.df$R0_median.t2 - alphaRstarR0ext.runs.df$Rstar1 
alphaRstarR0ext.runs.df$RstarR0sp2.t2 <- alphaRstarR0ext.runs.df$R0_median.t2 - alphaRstarR0ext.runs.df$Rstar2
par(mfrow=c(2,2))
hist(alphaRstarR0ext.runs.df$RstarR0sp1.t1)
hist(alphaRstarR0ext.runs.df$RstarR0sp2.t1)
hist(alphaRstarR0ext.runs.df$RstarR0sp1.t2)
hist(alphaRstarR0ext.runs.df$RstarR0sp2.t2)

###
# Next let's look at who is left after non-stationary period....
# And compare with the other runs
###
alphaRstarR0.calc.df.start <- df.all.long[which(df.all.long$jobID %in% alphaRstarR0.runs),]
alphaRstarR0.calc.df <- subset(alphaRstarR0.calc.df.start, period==2 & ncoexist>0)
# replace non-coexisting species' alpha and rstar values with NA
alphaRstarR0.calc.df$alpha1[which(alphaRstarR0.calc.df$coexist1==0)] <- NA
alphaRstarR0.calc.df$alpha2[which(alphaRstarR0.calc.df$coexist2==0)] <- NA
alphaRstarR0.calc.df$Rstar1[which(alphaRstarR0.calc.df$coexist1==0)] <- NA
alphaRstarR0.calc.df$Rstar2[which(alphaRstarR0.calc.df$coexist2==0)] <- NA
# repeat above for other comparable runs (without declining R0)
alphaRstar.calc.df.start <- df.all.long[which(df.all.long$jobID %in% alphaRstar.runs),]
alphaRstar.calc.df <- subset(alphaRstar.calc.df.start, period==2 & ncoexist>0)
alphaRstar.calc.df$alpha1[which(alphaRstar.calc.df$coexist1==0)] <- NA
alphaRstar.calc.df$alpha2[which(alphaRstar.calc.df$coexist2==0)] <- NA
alphaRstar.calc.df$Rstar1[which(alphaRstar.calc.df$coexist1==0)] <- NA
alphaRstar.calc.df$Rstar2[which(alphaRstar.calc.df$coexist2==0)] <- NA
# and the extrreme R0 declines
alphaRstarR0ext.calc.df.start <- df.all.long[which(df.all.long$jobID %in% alphaRstarR0ext.runs),]
alphaRstarR0ext.calc.df <- subset(alphaRstarR0ext.calc.df.start, period==2 & ncoexist>0)
# replace non-coexisting species' alpha and rstar values with NA
alphaRstarR0ext.calc.df$alpha1[which(alphaRstarR0ext.calc.df$coexist1==0)] <- NA
alphaRstarR0ext.calc.df$alpha2[which(alphaRstarR0ext.calc.df$coexist2==0)] <- NA
alphaRstarR0ext.calc.df$Rstar1[which(alphaRstarR0ext.calc.df$coexist1==0)] <- NA
alphaRstarR0ext.calc.df$Rstar2[which(alphaRstarR0ext.calc.df$coexist2==0)] <- NA

mean(c(alphaRstar.calc.df$Rstar1, alphaRstar.calc.df$Rstar2), na.rm=TRUE)
mean(c(alphaRstarR0.calc.df$Rstar1, alphaRstarR0.calc.df$Rstar2), na.rm=TRUE)
mean(c(alphaRstarR0ext.calc.df$Rstar1, alphaRstarR0ext.calc.df$Rstar2), na.rm=TRUE)
mean(c(alphaRstar.calc.df$alpha1, alphaRstar.calc.df$alpha2), na.rm=TRUE)
mean(c(alphaRstarR0.calc.df$alpha1, alphaRstarR0.calc.df$alpha2), na.rm=TRUE)
mean(c(alphaRstarR0ext.calc.df$alpha1, alphaRstarR0ext.calc.df$alpha2), na.rm=TRUE)

# We think tracking will be favored less as it makes the tracker use seedbank each year and thus
# it will 'blow through its seedbank' (Megan's words) in all those low R0 years
plot.paramdiffs.onepanel(alphaRstarR0.runs.df, "decliningR0/alphaRstarR0.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(alphaRstarR0.runs.df, "decliningR0/alphaRstarR0.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.manypanel.bfin(alphaRstarR0.runs.df, "decliningR0/alphaRstarR0.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft", colpalettehere)
plot.paramdiffs.onesp.bfin(alphaRstarR0.runs.df, "decliningR0/alphaRstarR0.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft", colpalettehere)
# and for the extreme ...
plot.paramdiffs.onepanel(alphaRstarR0ext.runs.df, "decliningR0/alphaRstarR0ext.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(alphaRstarR0ext.runs.df, "decliningR0/alphaRstarR0ext.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.manypanel.bfin(alphaRstarR0ext.runs.df, "decliningR0/alphaRstarR0ext.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft", colpalettehere)
plot.paramdiffs.onesp.bfin(alphaRstarR0ext.runs.df, "decliningR0/alphaRstarR0ext.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft", colpalettehere)
plot.paramdiffs.tworuntypes(alphaRstar.runs.df, alphaRstarR0ext.runs.df,
    "alphaRstar.runs", "alphaRstarR0ext.runs", "decliningR0/alphaRstar", "_comp.alpha.rstar", "ratio.alpha", "ratio.rstar", 
    cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")


stop(print("stopping here..."))


#################################
### Plots for the manuscript ###
#################################

library(RColorBrewer)
coexist.mscolz = brewer.pal(4, "Spectral")

symbolz <- c(8, 2, 0, 16)
leg.txt <- c("both species extirpated", "species 1 persists", "species 2 persists", "both species persist")

plot.paramdiffs.manuscript <- function(df, runname, figname, colname.x, colname.y, cex, pch,
        corner1.text, corner1.pos, corner2.text, corner2.pos){
    pdf(paste("graphs/modelruns/manuscript/", figname, ".pdf", sep=""),
       width=5, height=8)
    df0 <- subset(df, ncoexist.t2==0)
    df1 <- subset(df, ncoexist.t2==1)
    df2 <- subset(df, ncoexist.t2==2)
    df1.sp1 <- subset(df1, coexist1.t2==1)
    df1.sp2 <- subset(df1, coexist2.t2==1)
    plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
       ylab=colname.y, main="")
    abline(v=1, col="lightgray")
    abline(h=1, col="lightgray")
    fig_label(text=corner1.text, region="plot", pos=corner1.pos, cex=0.75)
    fig_label(text=corner2.text, region="plot", pos=corner2.pos, cex=0.75)
    points(df2[[colname.x]], unlist(df2[colname.y]),
        col=coexist.mscolz [3], pch=pch, cex=cex)
    points(df0[[colname.x]], unlist(df0[colname.y]),
        col=coexist.mscolz [1], pch=pch, cex=cex)
    points(df1.sp1[[colname.x]], unlist(df1.sp1[colname.y]),
        col=coexist.mscolz [2], pch=pch, cex=cex)
    points(df1.sp2[[colname.x]], unlist(df1.sp2[colname.y]),
        col=coexist.mscolz [2], pch=pch, cex=cex)
    legend("topright", leg.txt, pch=pch, col=coexist3col, bty="n")
    dev.off()
}

plot.paramdiffs.manuscript(alphaRstar.runs.df, "alphaRstar.runs", "alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, symbolz, "sp1 wins", "bottomleft", "sp2 wins", "topright")

plot.paramdiffs.manuscript(tauRstar.runs.df, "tauRstar.runs", "tauI.rstar", "ratio.tauI",
    "ratio.rstar", cexhere, symbolz, "sp1 wins", "bottomleft", "sp2 wins", "topright")


#################################
## Pull some of the R0 runs #####
## and check that they decline ##
#################################
if(lookdeeplyatR0runs){
folderID <- "933723"
samplerun <-  read.table(paste("output/OtherOut/envt/", folderID, "/EnvtParms_", folderID,
    "-1.txt", sep=""), header=TRUE) # comment.char = "", 
df.all <- data.frame(matrix(ncol=length(colnames(samplerun)), nrow=0))
file.names <- dir(paste("output/OtherOut/envt/", folderID, sep=""), pattern =".txt")
runse1 <- getfiles.envt(folderID, file.names, colnameshere) # this is SLOW
colnames(runse1) <- colnames(samplerun)
runse1$taskrunID <- paste(runse1$arrayID, runse1$runID, sep="-")

runshere <- unique(runse1$taskrunID)

# Yes, they decline!
par(mfrow=c(3,3))
for (i in 1:9){
    subby <- subset(runse1, taskrunID==unique(runse1$taskrunID)[i])
    plot(subby$R0~c(1:nrow(subby)))
}

# But tauP declines a little more 
quartz()
par(mfrow=c(3,3))
for (i in 1:9){
    subby <- subset(runse1, taskrunID==unique(runse1$taskrunID)[i])
    plot(subby$tauP~c(1:nrow(subby)))
}
}
