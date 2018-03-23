### Started 13 February 2018 ###
### By Lizzie ###

## Analyzing the runs data for the storage effect model! ##
## Last updated mid March 2018 ##

######################
### To do items!!! ###
######################
# Add to plots (see plotting section and notes at top of it)
# Make gi plots (need EnvtParms) and some Bfin plots
# tauIP (realized distance from tauI versus tauP) is more than tauIhat (so delete all the code for tauIhat)
# Make bout plots
######################

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## libraries
library(ggplot2)
require(plyr); require(dplyr); require(tidyr) # data formatting

## set working directory
setwd("~/Documents/git/projects/temporalvar/R")

source("sourcefiles/multiplot.R") # used in plot.params
source("sourcefiles/runanalysisfxs.R")

## flags for what to do
runbout <- FALSE

# cheap loop over the files for now
sruns <- c("36426477", "36511349","36691943", "36691954", "36691955")
nsruns <- c("36511352", "36511384", "36691956")

#####################################################
## Loop over all folders with stationary-only runs ##
#####################################################
for(folderIDhere in c(1:length(sruns))){
    
folderID <- sruns[folderIDhere]
samplerun <-  read.table(paste("output/", folderID, "/SummaryOut_", folderID,
    "-1.txt", sep=""), header=TRUE)
filenamestart <- c(paste("SummaryOut_", folderID, "-", sep=""))
colnameshere <- colnames(samplerun)
numhere <- c(1:20)

runs1 <- getfiles(folderID, filenamestart, numhere, colnameshere)
runs1$taskrunID <- paste(runs1$taskID, runs1$runID, sep="-")

    
##
## Data formatting 
##
df <- makediffs(runs1)
df.coexist <- subset(df, ncoexist==2)
print(paste("the current folder ID is", folderID, "the total rows are:", nrow(df), 
    "and the total coexisting rows are:", nrow(df.coexist), sep=" "))
 
##
## Get a list of X coexisting and not coexisting runs 
##
howmanytoget <- 10
df.nocoexist <- subset(df, ncoexist==0)
subsample <- rbind(df.coexist[1:howmanytoget,], df.nocoexist[1:howmanytoget,])
write.csv(subsample, file=paste("output/", folderID, "/", "subsample", ".csv", sep=""))

##
## Plotting 
##
## Look at the differences in parameters between the two species
# and color dots by coexisting yay or nay
plot.params(df, df.coexist, "alpha", "alpha1", "alpha2")
plot.params(df, df.coexist, "tauI", "tauIPini1_pre", "tauIPini2_pre")
plot.params(df, df.coexist, "germ", "g1mean_pre", "g2mean_pre")
plot.params(df, df.coexist, "c", "c1", "c2")
plot.params(df, df.coexist, "rstar", "Rstar1", "Rstar2")
    
## Look at parameter differences between species
plot.paramdiffs(df, df.coexist, "rstar_vs_tauI", "diff.tauIPini", "diff.rstar")
plot.paramdiffs(df, df.coexist, "rstar_vs_alpha", "diff.alpha", "diff.rstar")
plot.paramdiffs(df, df.coexist, "rstar_vs_tauIPini", "diff.tauIPini", "diff.rstar")
plot.paramdiffs(df, df.coexist, "g_vs_tauIPini", "diff.tauIPini", "diff.gmean")
plot.paramdiffs(df, df.coexist, "g_vs_alpha", "diff.alpha", "diff.gmean")
plot.paramdiffs(df, df.coexist, "g_vs_c", "diff.c", "diff.gmean")
plot.paramdiffs(df, df.coexist, "rstar_vs_tauI_ratios", "ratio.tauIini", "ratio.rstar")
plot.paramdiffs(df, df.coexist, "rstar_vs_alpha_ratios", "ratio.alpha", "ratio.rstar")
}



####################################################
## Loop over all folders with non-stationary runs ##
####################################################
for(folderIDhere in c(1:length(nsruns))){
    
folderID <- nsruns[folderIDhere]
samplerun <-  read.table(paste("output/", folderID, "/SummaryOut_", folderID,
    "-1.txt", sep=""), header=TRUE)
filenamestart <- c(paste("SummaryOut_", folderID, "-", sep=""))
colnameshere <- colnames(samplerun)
numhere <- c(1:20)

runs1 <- getfiles(folderID, filenamestart, numhere, colnameshere)
runs1$taskrunID <- paste(runs1$taskID, runs1$runID, sep="-")
    
##
## Data formatting 
##
df <- makediffs.ns(runs1)
df.coexist <- subset(df, ncoexist==2)
print(paste("the current folder ID is", folderID, "the total rows are:", nrow(df), 
    "and the total coexisting rows are:", nrow(df.coexist), sep=" "))
 
##
## Get a list of X coexisting and not coexisting runs 
##
howmanytoget <- 10
df.nocoexist <- subset(df, ncoexist==0)
subsample <- rbind(df.coexist[1:howmanytoget,], df.nocoexist[1:howmanytoget,])
write.csv(subsample, file=paste("output/", folderID, "/", "subsample", ".csv", sep=""))

##
## Plotting 
##
## Look at the differences in parameters between the two species
# and color dots by coexisting yay or nay
plot.params(df, df.coexist, "alpha", "alpha1", "alpha2")
plot.params(df, df.coexist, "tauI", "tauIPini1_pre", "tauIPini2_pre")
plot.params(df, df.coexist, "tauI", "tauIPns1_pre", "tauIPns2_pre")
plot.params(df, df.coexist, "tauI", "tauIPfin1_pre", "tauIPfin2_pre")
plot.params(df, df.coexist, "germ", "g1mean_pre", "g2mean_pre")
plot.params(df, df.coexist, "c", "c1", "c2")
plot.params(df, df.coexist, "rstar", "Rstar1", "Rstar2")
    
## Look at parameter differences between species
    ## Need to add in tauI ns and tauI fin!! ##
plot.paramdiffs(df, df.coexist, "rstar_vs_alpha", "diff.alpha", "diff.rstar")
plot.paramdiffs(df, df.coexist, "rstar_vs_tauIPini", "diff.tauIPini_pre", "diff.rstar")
plot.paramdiffs(df, df.coexist, "g_vs_tauIPini", "diff.tauIPini_pre", "diff.gmean")
plot.paramdiffs(df, df.coexist, "g_vs_alpha", "diff.alpha", "diff.gmean")
plot.paramdiffs(df, df.coexist, "g_vs_c", "diff.c", "diff.gmean")
plot.paramdiffs(df, df.coexist, "rstar_vs_tauI_ratios", "ratio.tauIini_pre", "ratio.rstar")
plot.paramdiffs(df, df.coexist, "rstar_vs_alpha_ratios", "ratio.alpha", "ratio.rstar")
}





stop(print("stopping here..."))
