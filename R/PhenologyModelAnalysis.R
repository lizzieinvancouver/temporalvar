### Started 13 February 2018 ###
### By Lizzie ###

## Analyzing the runs data for the storage effect model! ##
## Last updated mid March 2018 ##

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

## flags for what to do
runbout <- FALSE

# cheap loop over the files for now
runs3 <- c("51745826", "51745827", "51745828", "51745829") # new set as of 27 Aug 2018, including "51744988", "51744989"

# runinfo <- read.table("Table_of_RunParms.txt", skip=1, header=TRUE)

# varying:
# tauI, Rstar <- 41801567 <- 0 coexist
# alpha, Rstar <- 41801561 <- 1 coexist
# tauI, alpha <- 41801534 <- 179/4000 rows

# tauI, alpha with 10 rstar <- 41801557 <- 27/3205
# alpha with 10 rstar <- 41801558 <- 31/3436
# tauI with 10 rstar <- 41801559 <- 25/6189

# alpha, tauI, Rstar <- 41801560 <- 25/1885

# Rstar <-  41801578 <- 1/346
# alpha <- 41801579 <- 7/178
# alpha  <- 41801535 <- 219/4221
# tauI  <- 41801556 <- 83/6239

##############################
## Notes on other runs .... ##
##############################
# runs1 <- c("36426477", "36511349", "36511352", "36691943","36511384","36691954", "36691955", "36691956") # this is the set from Feb 2018
# runs2 <- c(paste("417251", 25:32, sep="")) # new set as of 16 Apr 2018
# EMPTY: "41725133", "41774659", "41775246", "41776446", "41779263","41779441"
# Note that , "41801581" ends at 33 (not 40)
# not coexisting after p2: "41801580"

#########################################
## Do some data reading and formatting ##
#########################################
runnow <- runs3 

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


##############################
## Set tauP for the graphs! ##
##############################
if(FALSE){
source("sourcefiles/analyses/tauP.R")
ggplot(tauP.plot, aes(x=tauP, fill=when)) + geom_density(alpha=0.25)
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

##
## Data formatting to compare species in histograms
## Get runs with coexist=2 in period 1 (stat), keep data only in period 2 (ns)
##
df.t2.wp <- subset(df, period==2)
# get just the period 2 runs that had coexistence in period 1
df.long <- df.t2.wp[which(df.t2.wp$taskrunID %in% unique(df.coexist1$taskrunID)),]
df.long.exist <- subset(df.long, coexist1==1 | coexist2==1)
df.long.noexist <- subset(df.long, coexist1==0 | coexist2==0)
}

###############
## Plotting! ##
###############
coexist3col <- add.alpha(c("firebrick", "dodgerblue", "seagreen"), alpha=0.4)
# col2rgb helps here ...
leg.txt <- c("2 survive", "1 left", "poof")

plot.paramdiffs(df.plot, "tauIP.alpha", "ratio.tauIP", "ratio.alpha")
plot.paramdiffs(df.plot, "tauI.alpha", "ratio.tauI", "ratio.alpha")

plot.paramdiffs(df.plot, "tauI.alpha.diff", "diff.tauIP", "diff.alpha")
plot.paramdiffs(df.plot, "tauI.alpha.diff", "diff.tauI", "diff.alpha")


## histograms
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


stop(print("stopping here..."))
