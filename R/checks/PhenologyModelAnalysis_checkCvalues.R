## Started 28 August 2018 ## 
## By Lizzie ##

## Check which c and R* values seem to lead to coexistence ## 
## Based off PhenologyModelAnalysis) ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## libraries
library(ggplot2)
require(plyr); require(dplyr); require(tidyr) # data formatting

## set working directory
setwd("~/Documents/git/projects/temporalvar/R")

## Get the f(x)s (extracted from sourcefiles/analyses/runanalysisfxs.R)
## Just grabbing them here in case they change someday ...

getfiles <- function(folderID, file.names, colnameshere){
    # numhere <- as.numeric(gsub("^[^-]*-([^.]+).*", "\\1", file.names))
    filepack <- lapply(file.names, function(file.names) {
    filename <- paste("output/SummaryFiles/", folderID, "/", file.names, sep="")
    dat <- read.table(filename, skip=1)
    names(dat) <- colnameshere
    return(data.frame(dat))
    })
    datahere <- do.call("rbind", filepack)
}


makediffs <- function(df){
    dathere <- df
    dathere$diff.c <-  dathere$c1-dathere$c2
    dathere$diff.rstar <-  dathere$Rstar1-dathere$Rstar2
    dathere$diff.tauI <-  dathere$tauI1-dathere$tauI2
    dathere$diff.tauIP <- dathere$tauIP1_mean-dathere$tauIP2_mean
    dathere$diff.alpha <-dathere$alpha1-dathere$alpha2
    dathere$ratio.c <-  dathere$c1/dathere$c2
    dathere$ratio.rstar <-  dathere$Rstar1/dathere$Rstar2
    dathere$ratio.tauI <-  dathere$tauI1/dathere$tauI2
    dathere$ratio.tauIP <- dathere$tauIP1_mean/dathere$tauIP2_mean
    dathere$ratio.alpha <-dathere$alpha1/dathere$alpha2
    return(dathere)
    }

##
## And now do some real work
## 
runs3 <- c("51803375", "51803287", "51803320", "51803342") 


runnow <- runs3

## Set up to paste together runs in one DF for now 
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

## Paste together runs in one DF in the loop!
df.all <- rbind(df.all, df)
    } 
plot(itertime~diff.rstar, data=df.all) # there is no relationship
summary(lm(itertime~diff.rstar, data=df.all))
ggplot(df.all, aes(y=itertime, x=diff.rstar)) + geom_point(aes(color=as.factor(ncoexist)))
ggplot(subset(df.all, ncoexist==2), aes(y=itertime, x=diff.rstar)) + geom_point() # there is no relationship
ggplot(df.all, aes(y=itertime, x=diff.rstar)) + geom_point(aes(color=as.factor(period))) # relationship!

df.p1 <- subset(df.all, period==1) # assumes all runs start with stationary,
# may miss some stationary periods but a good start
## runs with just rstar varying ...
# 51803342 has no tracking varying but does vary tauI so basically all the runs have something other than Rstar varying ...
# But all have Rstar varying! So we'll just use all?
# can code someday, get runparms and look for: varRstar2==NA and varRstar1==1 
plot(diff.rstar~ncoexist, df.p1)
rstar.summ.all <-
      ddply(df.all, c("ncoexist"), summarise,
      mean.diff.rstar = mean(diff.rstar),
      min.diff.rstar = min(diff.rstar),
      max.diff.rstar = max(diff.rstar),
      mean.diff.c = mean(diff.c),
      min.diff.c = min(diff.c),
      max.diff.c = max(diff.c),
      mean.Rstar1 = mean(Rstar1),
      mean.Rstar2 = mean(Rstar2),
      min.Rstar1 = min(Rstar1),
      min.Rstar2 = min(Rstar2), 
      max.Rstar1 = max(Rstar1),
      max.Rstar2 = max(Rstar2),
      mean.c1 = mean(c1),
      mean.c2 = mean(c2),
      min.c1 = min(c1),
      min.c2 = min(c2), 
      max.c1 = max(c1),
      max.c2 = max(c2))

rstar.summ.p1 <-
      ddply(df.p1, c("ncoexist"), summarise,
      mean.diff.rstar = mean(diff.rstar),
      min.diff.rstar = min(diff.rstar),
      max.diff.rstar = max(diff.rstar),
      mean.diff.c = mean(diff.c),
      min.diff.c = min(diff.c),
      max.diff.c = max(diff.c),
      mean.Rstar1 = mean(Rstar1),
      mean.Rstar2 = mean(Rstar2),
      min.Rstar1 = min(Rstar1),
      min.Rstar2 = min(Rstar2), 
      max.Rstar1 = max(Rstar1),
      max.Rstar2 = max(Rstar2),
      mean.c1 = mean(c1),
      mean.c2 = mean(c2),
      min.c1 = min(c1),
      min.c2 = min(c2), 
      max.c1 = max(c1),
      max.c2 = max(c2))


require(gridExtra)

cplot1 <- ggplot(df.all, aes(x=as.factor(ncoexist), y=diff.c)) + geom_boxplot(aes(fill=as.factor(period)))
cplot2 <- ggplot(df.all, aes(x=as.factor(ncoexist), y=c1)) + geom_boxplot(aes(fill=as.factor(period)))
cplot3 <- ggplot(df.all, aes(x=as.factor(ncoexist), y=c2)) + geom_boxplot(aes(fill=as.factor(period)))
rsplot1 <- ggplot(df.all, aes(x=as.factor(ncoexist), y=diff.rstar)) + geom_boxplot(aes(fill=as.factor(period)))
rsplot2 <- ggplot(df.all, aes(x=as.factor(ncoexist), y=Rstar1)) + geom_boxplot(aes(fill=as.factor(period)))
rsplot3 <- ggplot(df.all, aes(x=as.factor(ncoexist), y=Rstar2)) + geom_boxplot(aes(fill=as.factor(period)))

grid.arrange(cplot1, cplot2, cplot3, rsplot1, rsplot2, rsplot3, ncol=3)

