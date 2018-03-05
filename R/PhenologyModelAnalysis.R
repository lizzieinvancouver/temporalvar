### Started 13 February 2018 ###
### By Lizzie ###

## Analyzing the runs data for the storage effect model! ##
## Last updated late Feb 2018 ##
## Happy not-Valentine's Day! ##

######################
### To do items!!! ###
######################
# Change the filename structure (write into folders of run name IDs?)
# Add to plots (see plotting section and notes at top of it)
######################

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## libraries
library(ggplot2)
require(plyr); require(dplyr); require(tidyr) # data formatting

## set working directory
setwd("~/Documents/git/projects/temporalvar/R")

source("sourcefiles/multiplot.R")

## flags for what to do
runbfin <- FALSE
runbout <- FALSE

# cheap loop over the files for now
currentruns <- c("36426477", "36511349", "36511352", "36511384", "36691943",
    "36691954", "36691955", "36691956")
folderID <- "36511349"
samplerun <-  read.table(paste("output/", folderID, "/SummaryOut_", folderID,
    "-1.txt", sep=""), header=TRUE)
filenamestart <- c(paste("SummaryOut_", folderID, "-", sep=""))
colnameshere <- colnames(samplerun)
numhere <- c(1:20)

getfiles <- function(folderID, filenamestart, numhere, colnameshere){
    filepack <- lapply(numhere, function(numhere) {
    filename <- paste("output/", folderID, "/", filenamestart, numhere, ".txt", sep="")
    dat <- read.table(filename, skip=1)
    names(dat) <- colnameshere
    return(data.frame(dat))
    })
    datahere <- do.call("rbind", filepack)
}

runs1 <- getfiles(folderID, filenamestart, numhere, colnameshere)
runs1$taskrunID <- paste(runs1$taskID, runs1$runID, sep="-")


if(runbfin){
samplerunbfin <-  read.table(paste("output/", folderID, "/BfinN_", folderID,
    "-1.txt", sep=""), header=TRUE)
filenamestartBfin <- c(paste("BfinN_", folderID, "-", sep=""))
colnamesherebfin <- colnames(samplerunbfin)

runs1bfin <- getfiles(folderID, filenamestartBfin, numhere, colnamesherebfin)
runs1bfin$taskrunID <- paste(runs1bfin$taskID, runs1bfin$runID, sep="-")

}

#####################
## Data formatting ##
#####################
makediffs <- function(df){
    dathere <- df
    dathere$diff.alpha <- dathere$alpha1-dathere$alpha2
    dathere$diff.c <- dathere$c1-dathere$c2
    dathere$diff.rstar <- dathere$Rstar1-dathere$Rstar2
    dathere$diff.gmean <- dathere$g1mean_pre-dathere$g2mean_pre
    dathere$diff.tauIPini_pre <- dathere$tauIPini1_pre-dathere$tauIPini2_pre
    dathere$diff.tauIPns_pre <- dathere$tauIPns1_pre-dathere$tauIPns2_pre
    dathere$ratio.rstar <-  dathere$Rstar1/dathere$Rstar2
    dathere$ratio.tauIini_pre <- dathere$tauIPini1_pre/dathere$tauIPini2_pre
    dathere$ratio.tauIPns_pre <-dathere$tauIPns1_pre/dathere$tauIPns2_pre
    dathere$ratio.alpha <- dathere$alpha1/dathere$alpha2
    return(dathere)
    }


df <- makediffs(runs1)
df.coexist <- subset(df, ncoexist==2)

########################################################
## Get a list of X coexisting and not coexisting runs ##
########################################################
howmanytoget <- 10
df.nocoexist <- subset(df, ncoexist==0)
subsample <- rbind(df.coexist[1:howmanytoget,], df.nocoexist[1:howmanytoget,])
write.csv(subsample, file=paste("output/", folderID, "/", "subsample", ".csv", sep=""))

#################
## Plotting ##
#################

## Things I have not yet done yet for plotting ...
# Plots for Bfin, see PhenologyModelAnalysisADD.R 
# Plots for Bout (but see code to read it in below) #

## Code to read in Bout files (no figures for this yet) ##
if(runbfin){
getBoutfiles <- function(folderID, boutfilenamestart){
    boutfilenamestart <- 
    filename <- paste("output/", folderID, "/", "Bout/", "Bout*", ".txt", sep="")
    dat <- lapply(Sys.glob(filename), function(i) read.table(i, header=TRUE))
    # datahere <- do.call("rbind", dat)
    return(dat)
}

goober <- getBoutfiles(folderID, "Bout_36511349")
}


## Look at the differences in parameters between the two species
# and color dots by coexisting yay or nay

plot.params <- function(df, df.coexist, colname, colname.sp1, colname.sp2){
    pdf(paste("graphs/modelruns/params/Track_varR_2spp_", folderID, colname, "_compare.pdf", sep=""),
        width=9, height=4)
    plot.allruns <- ggplot(data=df, aes(x=df[colname.sp1], y=df[colname.sp2],
        colour=ncoexist)) +
        geom_point() +
        labs(x = colname.sp1, y=colname.sp2)
    plot.coexistruns <- ggplot(data=df.coexist, aes(x=df.coexist[colname.sp1], 
        y=df.coexist[colname.sp2], colour=ncoexist)) +
        geom_point() +
        labs(x = colname.sp1, y=colname.sp2)
    multiplot(plot.allruns, plot.coexistruns, cols=2)
    dev.off()
}


plot.params(df, df.coexist, "alpha", "alpha1", "alpha2")
plot.params(df, df.coexist, "tauI", "tauIPini1_pre", "tauIPini2_pre")
# plot.params(df, df.coexist, "tauIPini", "tauIPini.sp1", "tauIPini.sp2")
# plot.params(df, df.coexist, "tauIhat", "tauIhat.sp1", "tauIhat.sp2")
plot.params(df, df.coexist, "germ", "g1mean_pre", "g2mean_pre")
plot.params(df, df.coexist, "c", "c1", "c2")
plot.params(df, df.coexist, "rstar", "Rstar1", "Rstar2")

## Look at parameter differences between species

plot.paramdiffs <- function(df, df.coexist, figname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/Track_varR_2spp_", folderID, figname, ".pdf", sep=""),
        width=9, height=4)
    plot.allruns <- ggplot(data=df, aes(x=df[colname.x], y=df[colname.y],
        colour=ncoexist)) +
        geom_point() +
        labs(x = colname.x, y=colname.y)
    plot.coexistruns <- ggplot(data=df.coexist, aes(x=df.coexist[colname.x], 
        y=df.coexist[colname.y], colour=ncoexist)) +
        geom_point() +
        labs(x = colname.x, y=colname.y)
    multiplot(plot.allruns, plot.coexistruns, cols=2)
    dev.off()
}

plot.paramdiffs(df, df.coexist, "rstar_vs_tauI", "diff.tauIPini_pre", "diff.rstar")
plot.paramdiffs(df, df.coexist, "rstar_vs_alpha", "diff.alpha", "diff.rstar")
# plot.paramdiffs(df, df.coexist, "rstar_vs_tauIhat", "diff.tauIhat", "diff.rstar")
# plot.paramdiffs(df, df.coexist, "rstar_vs_tauIPini", "diff.tauIPini", "diff.rstar")
plot.paramdiffs(df, df.coexist, "g_vs_tauIPini", "diff.tauIPini_pre", "diff.gmean")
# plot.paramdiffs(df, df.coexist, "gi_vs_alpha", "diff.alpha", "diff.gi")
plot.paramdiffs(df, df.coexist, "g_vs_c", "diff.c", "diff.gmean")

plot.paramdiffs(df, df.coexist, "rstar_vs_tauI_ratios", "ratio.tauIini_pre", "ratio.rstar")
plot.paramdiffs(df, df.coexist, "rstar_vs_alpha_ratios", "ratio.alpha", "ratio.rstar")
# plot.paramdiffs(df, df.coexist, "rstar_vs_tauIhat_ratios", "tauIhat.ratio", "ratio.rstar")


