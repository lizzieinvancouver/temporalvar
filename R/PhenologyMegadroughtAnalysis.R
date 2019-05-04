### Started 30 December 2018 ###
### By Lizzie ###

## Analyzing the runs data for the megadroughts runs! ##
## Based heavily on PhenologyModelAnalysis.R ##

######################
### To do items!!! ###
######################

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## libraries
library(ggplot2)
require(plyr); require(dplyr); require(tidyr) # data formatting

## set working directory
setwd("~/Documents/git/projects/temporalvar/R")

source("sourcefiles/analyses/runanalysisfxsmega.R")

runshaveheader <- TRUE

# cheap loop over the files for now
runz <- c("9014727", "9014765", "9014815")

# c("1009960", "1009991", "1010022") # for these: Runs come in sets of 3, first is low tauI, second is mod tauI, and third is hi tauI
# c("1008362", "1008404", "1008425") # Zero coexistence!


#########################################
## Do some data reading and formatting ##
#########################################
runnow <- runz

## Setup for pasting runs together into one df (not sure we permanently want this but I want it now)
folderID <- runnow[1]
samplerun <-  read.table(paste("output/megadroughts/SummaryFiles/",folderID,"/SummaryOut_", folderID,
    "-1.txt", sep=""), header=TRUE) # comment.char = "", 
df.all <- data.frame(matrix(ncol=length(colnames(samplerun)), nrow=0))
colnames(df.all) <- colnames(samplerun)

for(folderIDhere in c(1:length(runnow))){
    
folderID <- runnow[folderIDhere] # folderID <- 41801534
if(runshaveheader){
    # Do we need the first line here?
    samplerun <-  read.table(paste("output/megadroughts/SummaryFiles/", folderID, "/SummaryOut_", folderID,
        "-1.txt", sep=""),  header=TRUE)
    file.names <- dir(paste("output/megadroughts/SummaryFiles/", folderID, sep=""), pattern =".txt")
    colnameshere <- colnames(samplerun)
    runs1 <- getfiles(folderID, file.names, colnameshere)
    runs1$taskrunID <- paste(runs1$taskID, runs1$runID, sep="-")
}

if(!runshaveheader){
    file.names <- dir(paste("output/megadroughts/SummaryFiles/", folderID, sep=""), pattern =".txt")
    colnameshere <- colnames(samplerun)
    runs1 <- getfiles(folderID, file.names, colnameshere)
    runs1$taskrunID <- paste(runs1$taskID, runs1$runID, sep="-")
}
 
 
##
## Data formatting to compare species pairs
##
df <- makediffs(runs1)
df.coexist <- subset(df, ncoexist==2)
print(paste("the current folder ID is", folderID, "the total rows are:", nrow(df), 
    "and the total coexisting rows are:", nrow(df.coexist), sep=" "))

# altogether now! pasting runs together into one df
df.all <- rbind(df.all, df)
    
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
df.all.t2 <- subset(df.all.t2, select=c("jobID", "taskID", "runID", "ncoexist",
    "coexist1", "coexist2",  "R0_mean", "R0_median", "R0_autocor", "taskrunID", "ratio.s", "ratio.phi"))
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
## Group the runs by what type they are 
##
lowtauI.runs <- runz[seq(1, length(runz), 3)]
modtauI.runs <- runz[seq(2, length(runz), 3)]
hitauI.runs <- runz[seq(3, length(runz), 3)]

## Grab just the wet period
df.all.stat <- subset(df.all, period==1)
lowtauI.stat.runs.df <- df.all.stat[which(df.all.stat$jobID %in% lowtauI.runs),]
modtauI.stat.runs.df <- df.all.stat[which(df.all.stat$jobID %in% modtauI.runs),]
hitauI.stat.runs.df <- df.all.stat[which(df.all.stat$jobID %in% hitauI.runs),]

## Grab the wet coexisting runs and all from the dry period
lowtauI.runs.df <- df.all.plot[which(df.all.plot$jobID %in% lowtauI.runs),]
modtauI.runs.df <- df.all.plot[which(df.all.plot$jobID %in% modtauI.runs),]
hitauI.runs.df <- df.all.plot[which(df.all.plot$jobID %in% hitauI.runs),]

# Check out s and phi range
par(mfrow=c(2,2))
hist(df.all$s1)
hist(df.all$s2)
hist(df.all$phi1)
hist(df.all$phi2)

par(mfrow=c(2,2))
hist(df.all$s1)
hist(df.all.coexist1$s2)
hist(df.all$phi1)
hist(df.all.coexist1$phi2)

# Check out tauI
mean(lowtauI.stat.runs.df$tauI1)
mean(modtauI.stat.runs.df$tauI1)
mean(hitauI.stat.runs.df$tauI1)

# Check other stuff

per1 <- subset(df.all, period==1)
per2 <- subset(df.all, period==2)

par(mfrow=c(2,2))
hist(per1$R0_median)
hist(per2$R0_median)
hist(per1$R0_autocor)
hist(per2$R0_autocor)

mean(per1$R0_median)
mean(per2$R0_median)
mean(per1$R0_autocor)
mean(per2$R0_autocor)

median(per1$R0_median)
median(per2$R0_median)
median(per1$R0_autocor)
median(per2$R0_autocor) # seems odd that megadrought period is negatively autocorrelated.



# Some plotting needs ... (probably more than we need)
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


## FIX MAKE DIFFS!

plot.paramdiffs.twopanel(df.all.stat, "newruns.stat", "_stau", "ratio.s.t1",
    "ratio.tau.t1", cexhere, pchhere, "? wins", "bottomleft", "? wins", "topright")

plot.paramdiffs.twopanel(df.all.plot, "newruns", "_stau", "ratio.s.t1",
    "ratio.tauI", cexhere, pchhere, "? wins", "bottomleft", "? wins", "topright")

plot.paramdiffs.twopanel(df.all.plot, "newruns", "_stau", "ratio.s",
    "ratio.tauIP.t1", cexhere, pchhere, "? wins", "bottomleft", "? wins", "topright")

if(FALSE){
plot.paramdiffs.twopanel(lowtauI.runs.df, "lowtauI.runs", "_sphi", "ratio.s.t1",
    "ratio.phi.t1", cexhere, pchhere, "? wins", "bottomleft", "? wins", "topright")
plot.paramdiffs.twopanel(modtauI.runs.df, "modtauI.runs", "_sphi", "ratio.s.t1",
    "ratio.phi.t1", cexhere, pchhere, "? wins", "bottomleft", "? wins", "topright")
plot.paramdiffs.twopanel(hitauI.runs.df, "hitauI.runs", "_sphi", "ratio.s.t1",
    "ratio.phi.t1", cexhere, pchhere, "? wins", "bottomleft", "? wins", "topright")
}
