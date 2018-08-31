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

# cheap loop over the files for now
runz <- c("51803287", "51803320", "51803342",  "51803375",
    "51893537", "51893598", "51893656", "51893711", 
    "51995069", "51995121", "51995125", "51995137",
    "52031904", "52031950", "52031996") # missing 52031833 which should vary everything

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
runnow <- runz

## Setup for pasting runs together into one df (not sure we permanently want this but I want it now)
folderID <- runnow[1] 
samplerun <-  read.table(paste("output/SummaryFiles/", folderID, "/SummaryOut_", folderID,
    "-1.txt", sep=""), header=TRUE)
df.all <- data.frame(matrix(ncol=length(colnames(samplerun)), nrow=0))
colnames(df.all) <- colnames(samplerun)

for(folderIDhere in c(1:length(runnow))){
    
folderID <- runnow[folderIDhere] # folderID <- 41801534
samplerun <-  read.table(paste("output/SummaryFiles/", folderID, "/SummaryOut_", folderID,
    "-1.txt", sep=""), header=TRUE)
# filenamestart <- c(paste("SummaryOut_", folderID, "-", sep=""))
file.names <- dir(paste("output/SummaryFiles/", folderID, sep=""), pattern =".txt")
colnameshere <- colnames(samplerun)
runs1 <- getfiles(folderID, file.names, colnameshere)
runs1$taskrunID <- paste(runs1$taskID, runs1$runID, sep="-")
 
##
## Data formatting to compare species pairs
##
df <- makediffs(runs1)
df.coexist <- subset(df, ncoexist==2)
print(paste("the current folder ID is", folderID, "the total rows are:", nrow(df), 
    "and the total coexisting rows are:", nrow(df.coexist), sep=" "))

# altogether now! pasting runs together into one df
df.all <- rbind(df.all, df)

##############################
## Set tauP for the graphs! ##
##############################
if(FALSE){
source("sourcefiles/analyses/tauP.R")
dfhere <- taualphaRstar.runs.df # tauRstar.runs.df

df2 <- subset(dfhere, ncoexist.t2==2) #
plot(nhere, tauPfin, type="l")
lines(nhere, tauP, type="l")
lines(density(c(dfhere$tauI1, dfhere$tauI2)), col="lightblue")
lines(density(c(df2$tauI1, df2$tauI2)), col="red")

quartz()
plot(nhere, tauPfin, type="l")
lines(nhere, tauP, type="l")
lines(density(pmax(dfhere$tauI1, dfhere$tauI2)), col="lightblue")
lines(density(pmax(df2$tauI1, df2$tauI2)), col="red")

quartz()
plot(nhere, tauPfin, type="l")
lines(nhere, tauP, type="l")
lines(density(pmin(dfhere$tauI1, dfhere$tauI2)), col="lightblue")
lines(density(pmin(df2$tauI1, df2$tauI2)), col="red")

# ggplot(tauP.plot, aes(x=tauP, fill=when)) + geom_density(alpha=0.25)
}

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
df.all.coexist1 <- subset(df.all, ncoexist==2 & period==1)
df.all.t2 <- subset(df.all, period==2)
df.all.t2 <- subset(df.all.t2, select=c("jobID", "taskID", "runID", "ncoexist", "taskrunID"))
df.all.plot <- merge(df.all.coexist1, df.all.t2, by=c("jobID", "taskID", "runID", "taskrunID"),
    all.x=TRUE, all.y=FALSE, suffixes=c(".t1", ".t2"))


##
## Data formatting to compare species in histograms
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

###############
## Plotting! ##
###############
coexist3col <- add.alpha(c("firebrick", "dodgerblue", "seagreen"), alpha=0.4)
# col2rgb helps here ...
leg.txt <- c("poof", "1 left", "2 survive")

## histograms DEAL WITH!!! 
breaknum <- 20
hist(c(df.long.exist$tauI1, df.long.exist$tau2), xlim=c(0,1), ylim=c(0,30),
     breaks=breaknum, col=coexist3col[3], main="", xlab="number")
par(new=TRUE)
hist(c(df.long.noexist$tauI1, df.long.noexist$tau2), xlim=c(0,1), ylim=c(0,30),
     breaks=breaknum, col=coexist3col[1], main="", xlab="", ylab="") 

plot.histograms(df.long.exist, df.long.noexist, "tauI", "tauI1", "tauI2",
    coexist3col, seq(from=0, to=1, by=0.05), c(0,1), c(0,40))

###
# ggplot(df.plot, aes(x=ratio.tauIP, color=as.factor(ncoexist.t2), fill=as.factor(ncoexist.t2))) +
#    geom_histogram(alpha=0.5, position="identity")

# ggplot() + geom_density(data=tauP.plot, aes(x=tauP), alpha=0.25) +
#    geom_histogram(data=df.plot, aes(x=tauI1, color=as.factor(ncoexist.t2), fill=as.factor(ncoexist.t2)))
###

# 3 and 7 are NOT varying tracking
tauRstar.runs <- runz[c(3,7, 11)] # tauI and Rstar tradeoff
tauRstar.runs.df <- df.all.plot[which(df.all.plot$jobID %in% tauRstar.runs),]

plot.paramdiffs.onepanel(tauRstar.runs.df, "tauRstar.runs", "tauIP.rstar", "ratio.tauIP",
    "ratio.rstar")
plot.paramdiffs.twopanel(tauRstar.runs.df, "tauRstar.runs", "tauIP.rstar", "ratio.tauIP",
    "ratio.rstar")

# 2 and 6 are NOT varying tauI
alphaRstar.runs <- runz[c(2,6,10)] # tracking and Rstar tradeoff
alphaRstar.runs.df <- df.all.plot[which(df.all.plot$jobID %in% alphaRstar.runs),]

# alphaRstar.run.dfs$newratio.alpha <- 1/alphaRstar.runs.df$ratio.alpha

plot.paramdiffs.onepanel(alphaRstar.runs.df, "alphaRstar.runs", "alpha.rstar", "ratio.alpha",
    "ratio.rstar")
plot.paramdiffs.twopanel(alphaRstar.runs.df, "alphaRstar.runs", "alpha.rstar", "ratio.alpha",
    "ratio.rstar")


# 4 and 8 keeps R* the same across species pairs
taualpha.runs <- runz[c(4, 8)] 
taualpha.runs.df <- df.all.plot[which(df.all.plot$jobID %in% taualpha.runs),]

plot.paramdiffs.onepanel(taualpha.runs.df, "taualpha.runs", "alpha.tauIP",
    "ratio.alpha", "ratio.tauIP")
plot.paramdiffs.twopanel(taualpha.runs.df, "taualpha.runs", "alpha.tauIP",
    "ratio.alpha", "ratio.tauIP")


# 1 and 5 are varying everything
taualphaRstar.runs <- runz[c(1,5, 9)] 
taualphaRstar.runs.df <- df.all.plot[which(df.all.plot$jobID %in% taualphaRstar.runs),]

plot.paramdiffs.onepanel(taualphaRstar.runs.df, "taualphaRstar.runs", "tauIP.rstar",
    "ratio.tauIP", "ratio.rstar")
plot.paramdiffs.onepanel(taualphaRstar.runs.df, "taualphaRstar.runs", "alpha.rstar",
    "ratio.alpha", "ratio.rstar")
plot.paramdiffs.onepanel(taualphaRstar.runs.df, "taualphaRstar.runs", "alpha.tauIP",
    "ratio.alpha", "ratio.tauIP")
plot.paramdiffs.twopanel(taualphaRstar.runs.df, "taualphaRstar.runs", "tauIP.rstar",
    "ratio.tauIP", "ratio.rstar")
plot.paramdiffs.twopanel(taualphaRstar.runs.df, "taualphaRstar.runs", "alpha.rstar",
    "ratio.alpha", "ratio.rstar")
plot.paramdiffs.twopanel(taualphaRstar.runs.df, "taualphaRstar.runs", "alpha.tauIP",
    "ratio.alpha", "ratio.tauIP")

# compare these to their equivalents 
plot.paramdiffs.fixedxy(taualphaRstar.runs.df, "taualphaRstar.runs", "tauIP.rstar",
    "ratio.tauIP", "ratio.rstar", tauRstar.runs.df)
plot.paramdiffs.fixedxy(taualphaRstar.runs.df, "taualphaRstar.runs", "alpha.rstar",
    "ratio.alpha", "ratio.rstar", alphaRstar.runs.df)
plot.paramdiffs.fixedxy(taualphaRstar.runs.df, "taualphaRstar.runs", "alpha.tauIP",
    "ratio.tauIP", "ratio.rstar", taualpha.runs.df)


# and color code by third variable
plot.paramdiffs.colorbyzvar(taualphaRstar.runs.df, "taualphaRstar.runs", "tauIP.rstar",
    "ratio.tauIP", "ratio.rstar", "ratio.alpha", "3 traits vary: survived after stat", 1)
plot.paramdiffs.colorbyzvar(taualphaRstar.runs.df, "taualphaRstar.runs", "alpha.rstar",
    "ratio.alpha", "ratio.rstar", "ratio.tauIP", "3 traits vary: survived after stat", 1)
plot.paramdiffs.colorbyzvar(taualphaRstar.runs.df, "taualphaRstar.runs", "alpha.tauIP",
     "ratio.alpha", "ratio.tauIP", "ratio.rstar", "3 traits vary: survived after stat", 1)


stop(print("stopping here..."))
