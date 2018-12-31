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
runz <- c("58137170", "58137225", "58137269", "58378327") # not the complete list


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
    "coexist1", "coexist2", "taskrunID", "ratio.tauIP"))
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
