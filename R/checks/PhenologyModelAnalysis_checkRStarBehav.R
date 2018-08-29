## Started 29 August 2018 ## 
## By Lizzie ##

## Check dynamics when R* is reached ## 
## Based somewhat off PhenologyModelAnalysis and # PhenologyModelAnalysis_intraann.R (in the archive) ##

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

## Not that it matters but I might as well pick a coexisting run!
# 51803342-1-56 has two coexisting ...
summjob <- read.table("output/SummaryFiles/51803342/SummaryOut_51803342-1.txt", header=TRUE)
# Get Bout, even just getting one is slow!
folderID <- 51803342
numherebout <- "51803342-1"
whichrun <- "56"
bouthere <-  read.table(paste("output/Bout/", folderID, "/Bout_", numherebout, "-",
    whichrun, ".txt", sep=""), header=TRUE)
table(bouthere$yr)

par(mfrow=c(3,1))

yrhere <- 139 # ends before timeout (1001)
# yrhere <- 140 # times out example
subbout <- subset(bouthere, yr==yrhere)
summrun <- subset(summjob, runID==whichrun)

subbout[nrow(subbout),]
summrun

plot(R~time, data=subbout)
abline(0, summrun$Rstar1[1], col="red") # lower Rstar
abline(0, summrun$Rstar2[1], col="blue") 

plot(B1~time, data=subbout, col="red", type="l")
plot(B2~time, data=subbout, col="blue", type="l")

# all years
yrends <- aggregate(bouthere[c("B1", "B2")], bouthere["yr"], FUN=max)
plot(B1~yr, data=yrends, col="red", type="l")
plot(B2~yr, data=yrends, col="blue", type="l")
