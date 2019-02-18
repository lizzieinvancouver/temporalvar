### Started 13 February 2018 ###
### By Lizzie ###

## Analyzing the runs data for the storage effect model! ##
## Last updated August 2018 ##

######################
### To do items!!! ###
######################
# https://www.theanalysisfactor.com/r-tutorial-part-12/
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

runshaveheader <- TRUE

# cheap loop over the files for now
runz2 <- c("858179", "858221", "858241", "858262") # "858282" R0 varies

# c("62025212", "62025233", "62025263", "62025274") # no header row

runz <- c("51803287", "51803320", "51803342",  "51803375",
    "51893537", "51893598", "51893656", "51893711", 
    "51995069", "51995121", "51995125", "51995137",
    "52031904", "52031950", "52031996") # missing 52031833 which should vary everything

# Remember to update below under plotting-related formatting ... 
# 1, 5, 9 are varying everything
# 2, 6, 10, 13 are NOT varying tauI
# 3, 7, 11, 14 are NOT varying tracking
# 4, 8, 12, 15 keeps R* the same across species pairs

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
runnow <- runz2

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
    "coexist1", "coexist2", "taskrunID", "ratio.tauIP", "tauIP1_mean",
    "tauIP2_mean", "diff.bfinslopes", "slopeBfin1", "slopeBfin2", 
    "minslopeBfin"))
df.all.plot <- merge(df.all.coexist1, df.all.t2, by=c("jobID", "taskID", "runID", "taskrunID"),
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
tauRstar.check <-  df.all[which(df.all$jobID %in% runz2[c(3)]),] # tauI and Rstar tradeoff runs
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
if(FALSE){
tauRstar.runs <- runz[c(3,7,11,14)] # NOT varying tracking: tauI and Rstar tradeoff
alphaRstar.runs <- runz[c(2,6,10,13)] # NOT varying tauI: tracking and Rstar tradeoff
taualpha.runs <- runz[c(4,8,12,15)] # keeps R* the same across species pairs
taualphaRstar.runs <- runz[c(1,5,9)] # varying everything (tauI, alpha, Rstar)
}

tauRstar.runs <- runz2[c(3)] # NOT varying tracking: tauI and Rstar tradeoff
alphaRstar.runs <- runz2[c(2)] # NOT varying tauI: tracking and Rstar tradeoff
taualpha.runs <- runz2[c(4)] # keeps R* the same across species pairs
taualphaRstar.runs <- runz2[c(1)] # varying everything (tauI, alpha, Rstar)

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

if(FALSE){ # why do I have this here again? 
tauRstar.runs <- runz[c(3,7,11,14)] # NOT varying tracking: tauI and Rstar tradeoff
alphaRstar.runs <- runz[c(2,6,10,13)] # NOT varying tauI: tracking and Rstar tradeoff
taualpha.runs <- runz[c(4,8,12,15)] # keeps R* the same across species pairs
taualphaRstar.runs <- runz[c(1,5,9)] # varying everything (tauI, alpha, Rstar)
}

## Do some counts
# How many species left after non-stat? First: (# spp left at end of non-stat)/(# of spp left at end of stat)
sum(tauRstar.runs.df$ncoexist.t2)/sum(tauRstar.stat.runs.df$ncoexist) 
sum(alphaRstar.runs.df$ncoexist.t2)/sum(alphaRstar.stat.runs.df$ncoexist)
sum(taualpha.runs.df$ncoexist.t2)/sum(taualpha.stat.runs.df$ncoexist)
sum(taualphaRstar.runs.df$ncoexist.t2)/sum(taualphaRstar.stat.runs.df$ncoexist)
# Second: (# spp left at end of non-stat)/(# of co-existing spp left at end of stat)
sum(tauRstar.runs.df$ncoexist.t2)/sum(tauRstar.runs.df$ncoexist.t1) 
sum(alphaRstar.runs.df$ncoexist.t2)/sum(alphaRstar.runs.df$ncoexist.t1)
sum(taualpha.runs.df$ncoexist.t2)/sum(taualpha.runs.df$ncoexist.t1)
sum(taualphaRstar.runs.df$ncoexist.t2)/sum(taualphaRstar.runs.df$ncoexist.t1)
# What is average alpha before and after stat?
get.mean.alphavalues(df.all[which(df.all$jobID %in% alphaRstar.runs),])
get.mean.alphavalues.ns(df.all[which(df.all$jobID %in% alphaRstar.runs),])
get.mean.alphavalues(df.all[which(df.all$jobID %in% taualpha.runs),])
get.mean.alphavalues.ns(df.all[which(df.all$jobID %in% taualpha.runs),])
get.mean.alphavalues(df.all[which(df.all$jobID %in% taualphaRstar.runs),])
get.mean.alphavalues.ns(df.all[which(df.all$jobID %in% taualphaRstar.runs),])


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
library(RColorBrewer)
cols = brewer.pal(4, "RdBu")
# Define colour pallete
pal = colorRampPalette(c("blue", "red"))
# Use the following line with RColorBrewer
colpalettehere = colorRampPalette(cols)

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

plot.paramdiffs.onepanel(alphaRstar.runs.df, "alphaRstar.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(alphaRstar.runs.df, "alphaRstar.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.manypanel.bfin(alphaRstar.runs.df, "alphaRstar.runs", "_alpha.rstar", "ratio.alpha",
    "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft", colpalettehere)

plot.paramdiffs.onepanel(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP.t1.",
    "ratio.alpha", "ratio.tauIP.t1", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.onepanel(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP.t2", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP.t2", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP.t1.",
    "ratio.alpha", "ratio.tauIP.t1", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualpha.runs.df, "taualpha.runs", "_alpha.tauI",
    "ratio.alpha", "ratio.tauI", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel.fixedxy(taualpha.runs.df, "taualpha.runs", "_alpha.tauI",
    "ratio.alpha", "ratio.tauI", cexhere, pchhere, c(0,3), c(0,5),
    "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel.fixedxy(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP.t2.",
    "ratio.alpha", "ratio.tauIP.t2", cexhere, pchhere, c(0,3), c(0,5),
    "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel.fixedxy(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP.t1", cexhere, pchhere, c(0,3), c(0,5),
    "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.manypanel.bfin(taualpha.runs.df, "taualpha.runs", "_alpha.tauI",
    "ratio.alpha", "ratio.tauI", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft",
    colpalettehere)
# some extras for this tricky run ....
taualpha.runs.df$ratio.tauIP.t1t2 <- (taualpha.runs.df$tauIP1_mean.t1-taualpha.runs.df$tauIP1_mean.t2)/
    (taualpha.runs.df$tauIP2_mean.t1-taualpha.runs.df$tauIP2_mean.t2)
plot.paramdiffs.twopanel(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP.t1t2",
    "ratio.alpha", "ratio.tauIP.t1t2", cexhere, pchhere, "sp1 wins", "bottomright",
    "sp2 wins", "topleft")
plot.paramdiffs.twopanel.fixedxy(taualpha.runs.df, "taualpha.runs", "_alpha.tauIP.t1t2",
    "ratio.alpha", "ratio.tauIP.t1t2", cexhere, pchhere, c(0,3), c(-20,20),
    "sp1 wins", "bottomright", "sp2 wins", "topleft")


# Three things vary
plot.paramdiffs.onepanel(taualphaRstar.runs.df, "taualphaRstar.runs", "_tauIP.rstar",
    "ratio.tauIP.t2", "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.onepanel(taualphaRstar.runs.df, "taualphaRstar.runs", "_tauI.rstar",
    "ratio.tauI", "ratio.rstar", cexhere, pchhere, "??", "bottomleft", "??", "topright")
plot.paramdiffs.onepanel(taualphaRstar.runs.df, "taualphaRstar.runs", "_alpha.rstar",
    "ratio.alpha", "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.onepanel(taualphaRstar.runs.df, "taualphaRstar.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP.t2", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualphaRstar.runs.df, "taualphaRstar.runs", "_tauIP.rstar",
    "ratio.tauIP.t2", "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomleft", "sp2 wins", "topright")
plot.paramdiffs.twopanel(taualphaRstar.runs.df, "taualphaRstar.runs", "_alpha.rstar",
    "ratio.alpha", "ratio.rstar", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.twopanel(taualphaRstar.runs.df, "taualphaRstar.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP.t2", cexhere, pchhere, "sp1 wins", "bottomright", "sp2 wins", "topleft")


# compare these to their equivalents 
plot.paramdiffs.fixedxy(taualphaRstar.runs.df, "taualphaRstar.runs", "_tauIP.rstar",
    "ratio.tauIP.t2", "ratio.rstar", tauRstar.runs.df, 1, 16, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.fixedxy(taualphaRstar.runs.df, "taualphaRstar.runs", "_tauI.rstar",
    "ratio.tauI", "ratio.rstar", tauRstar.runs.df, 1, 16, "??", "bottomleft", "??", "topright")
plot.paramdiffs.fixedxy(taualphaRstar.runs.df, "taualphaRstar.runs", "_alpha.rstar",
    "ratio.alpha", "ratio.rstar", alphaRstar.runs.df, 1, 16, "sp1 wins", "bottomright", "sp2 wins", "topleft")
plot.paramdiffs.fixedxy(taualphaRstar.runs.df, "taualphaRstar.runs", "_alpha.tauIP",
    "ratio.alpha", "ratio.tauIP.t2", taualpha.runs.df, 1, 16, "sp1 wins", "bottomright", "sp2 wins", "topleft")


# and color code by third variable
plot.paramdiffs.colorbyzvar(taualphaRstar.runs.df, "taualphaRstar.runs", "_tauIP.rstar",
    "ratio.tauIP.t2", "ratio.rstar", "ratio.alpha", "3 traits vary: survived after stat", 1)
plot.paramdiffs.colorbyzvar(taualphaRstar.runs.df, "taualphaRstar.runs", "_alpha.rstar",
    "ratio.alpha", "ratio.rstar", "ratio.tauIP.t2", "3 traits vary: survived after stat", 1)
plot.paramdiffs.colorbyzvar(taualphaRstar.runs.df, "taualphaRstar.runs", "_alpha.tauIP",
     "ratio.alpha", "ratio.tauIP.t2", "ratio.rstar", "3 traits vary: survived after stat", 1)


stop(print("stopping here..."))
