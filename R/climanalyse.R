## Started 17 August 2019 ##
## Home alone on a gray Saturday in August ##

## Code to look at climate shifts in timing ... ##

# housekeeping
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/temporalvar/R") 
} else setwd("~/Documents/git/boopboop")

# Libraries
library(ggplot2)

# What to do ...
runMarseille <- FALSE
runSantaBarbara <- FALSE
runSanDiego <- FALSE
runMunich <- TRUE
runBoston <- TRUE
runButte <- TRUE

# f(xs)
# This f(x) is not fast, but it works
makestartdatecounter <- function(dater, monthcol, daycol, countercol, whatmon, whatday){
    dater[[countercol]][which(dater[[monthcol]]==whatmon & dater[[daycol]]==whatday)] <- 0
    for(j in c(2:nrow(dater))){
        if(is.na(dater[[countercol]][j])){
        dater[[countercol]][j] <- dater[[countercol]][j-1]+1
        }
     }
    return(dater)
}

addwhen <- function(df){
    # add when
    df[["when"]][which(df[["year"]]<1980)] <- "before 1980"
    df[["when"]][which(df[["year"]]>1979)] <- "after 1980"
    return(df) 
    }


makethreshold.data <- function(dater, temp.col, thresh){ 
    ifelse(dater[[temp.col]]>thresh,
       (dater[[temp.col]]-thresh), 0)
  }

## this f(x) adds up gdd
## requires data ordered by doy (I do this in the loop below)
## this f(x) returns the value while treating NA as zeroes
# needstartdate is when you require data to start that year, otherwise it returns NA
# for example, 5 means you need data that starts before 5 January to actually count
makegdd.data.skipNA <- function(dater, gdd.col, doy.col, startdate, needstartdate){
     saveme <- c()
     for(i in 1:nrow(dater)){
     # start the counter with the starting doy of the data ...
     # j <- dater[[doy.col]][1]
     # deal with cases where the data start after Jan 1
     if (dater[[doy.col]][1]>needstartdate) saveme[i] <- NA
     else
     # deal with cases where the entire column is NA
     if (sum(is.na(dater[[gdd.col]]))==length(dater[[gdd.col]])) saveme[i] <- NA
     else
     # deal with cases before startdate
     if (dater[[doy.col]][i]<startdate) saveme[i] <- NA
     else
     # okay, finally calculate the GDD
     if (dater[[doy.col]][i]==startdate) saveme[i] <- (dater[[gdd.col]][i])
     else
     # if a cell is NA, just add 0 instead of the cell
     if (is.na(dater[[gdd.col]][i])) saveme[i] <- (0+saveme[i-1])
     else
     saveme[i] <- (dater[[gdd.col]][i]+saveme[i-1])
 }
 return(saveme)
}


if(runMarseille){
mont <- read.delim("climdata/marseille_precip_ed.txt", sep="\t", header=FALSE)
names(mont) <- c("year", "month", "day", "precip")

hist(mont$precip, breaks=1000, xlim=c(0,20))
mont.nonzerorain <- subset(mont, precip>0)
hist(mont.nonzerorain$precip, breaks=100)
hist(mont.nonzerorain$precip, breaks=100, xlim=c(0,20))

mont.annualprecip <- aggregate(mont["precip"], mont["year"], FUN=sum)
plot(density(mont.annualprecip$precip))
}


if(runSantaBarbara){
sb <- read.delim("climdata/santabarbara_precip_ed.txt", sep="\t", header=FALSE)
names(sb) <- c("year", "month", "day", "precip")

hist(sb$precip, breaks=1000, xlim=c(0,20))
sb.nonzerorain <- subset(sb, precip>0)
hist(sb.nonzerorain$precip, breaks=100)
hist(sb.nonzerorain$precip, breaks=100, xlim=c(0,20))

sb.annualprecip <- aggregate(sb["precip"], sb["year"], FUN=sum)
plot(density(sb.annualprecip$precip))
}


if(runSanDiego){
## Get the data
sd <- read.delim("climdata/sandiego_precip_ed.txt", sep="\t", header=FALSE)
names(sd) <- c("year", "month", "day", "precip")
    
sd.annualprecip <- aggregate(sd["precip"], sd["year"], FUN=sum)
plot(density(sd.annualprecip$precip))


# Add in some critical missing dates (onxw nwwsws)
# sd <- sd[order(sd$year, sd$month, sd$day),]

sd$seasonday <- NA
sd <- makestartdatecounter(sd, "month", "day", "seasonday", 8, 31)
subset(sd, month==8 & day==31 & seasonday>400) # check for missing years

hist(sd$precip, breaks=1000, xlim=c(0,20))
sd.nonzerorain <- subset(sd, precip>0)
hist(sd.nonzerorain$precip, breaks=100)
hist(sd.nonzerorain$precip, breaks=100, xlim=c(0,20))

## Hmm, let's just try 5 mm days
sd.firstrain <- sd[1,]
rainydays <- subset(sd, precip>10)
for(j in unique(rainydays$year)){
    subby <- subset(rainydays, year==j)
    minday <- subby[which(subby$seasonday==min(subby$seasonday, na.rm=TRUE)),]
    sd.firstrain <- rbind(sd.firstrain, minday)
   }
sd.firstrain <- sd.firstrain[-1,]

sd.firstrain$when <- NA
sduse <- addwhen(sd.firstrain)
nrow(subset(sduse, when=="before 1980"))
nrow(subset(sduse, when=="after 1980"))
summary(lm(seasonday~when, sduse))
summary(lm(precip~when, sduse))

# check that all years are represented...
nrow(sd.firstrain)
length(unique(sd$year))

hist(sd.firstrain$seasonday, breaks=30)
plot(density(sd.firstrain$seasonday))

ggplot(sduse, aes(x=seasonday, fill=when)) + geom_density(alpha=0.4)
ggplot(sduse, aes(x=seasonday, fill=when)) + geom_histogram(alpha=0.4, aes(y = ..density..), position = 'identity')

ggplot(sduse, aes(x=precip, fill=when)) + geom_density(alpha=0.4)
ggplot(sduse, aes(x=precip, fill=when)) + geom_histogram(alpha=0.4, aes(y = ..density..), position = 'identity')

sd$date <- as.Date(paste(sd$year, sd$month, sd$day, sep="-"), format="%Y-%m-%d")
plot(precip~date, data=sd)
}


if(runMunich){
## Get the data
mun <- read.delim("climdata/munichmin_ed.txt", sep="\t", header=FALSE)
names(mun) <- c("year", "month", "day", "tmin")

mun$date <- as.Date(paste(mun$year, mun$month, mun$day, sep="-"), format="%Y-%m-%d")
mun$doy <- as.numeric(format(mun$date, "%j"))

mun$tthres <- makethreshold.data(mun, "tmin", 0)
mun$gdd <- makegdd.data.skipNA(mun, "tthres", "doy", 1, 2)

## Day 200 GDD is reached
munwarm <- mun[1,]
munwarmall <- subset(mun, gdd>200)
for(j in unique(munwarmall$year)){
    subby <- subset(munwarmall, year==j)
    minday <- subby[which(subby$doy==min(subby$doy, na.rm=TRUE)),]
    munwarm <- rbind(munwarm, minday)
   }
munwarm <- munwarm[-1,]

munwarm$when <- NA
munuse <- addwhen(munwarm)
nrow(subset(munuse, when=="before 1980"))
nrow(subset(munuse, when=="after 1980"))
summary(lm(doy~when, munuse))

hist(munuse$doy, breaks=30)
plot(density(munuse$doy))

ggplot(munuse, aes(x=doy, fill=when)) + geom_density(alpha=0.4)
ggplot(munuse, aes(x=doy, fill=when)) + geom_histogram(alpha=0.4, aes(y = ..density..), position = 'identity')

# Get only 40 first years
munusesm <- subset(munuse, year<1919 | year>1979)
nrow(subset(munusesm, when=="before 1980"))
nrow(subset(munusesm, when=="after 1980"))
munwhen <- ggplot(munusesm, aes(x=doy, fill=when)) + geom_density(alpha=0.4)    
}


if(runBoston){
## Get the data
bos <- read.delim("climdata/loganmin_ed.txt", sep="\t", header=FALSE)
names(bos) <- c("year", "month", "day", "tmin")

bos$date <- as.Date(paste(bos$year, bos$month, bos$day, sep="-"), format="%Y-%m-%d")
bos$doy <- as.numeric(format(bos$date, "%j"))

bos$tthres <- makethreshold.data(bos, "tmin", 0)
bos$gdd <- makegdd.data.skipNA(bos, "tthres", "doy", 1, 2)

## Day 200 GDD is reached
boswarm <- bos[1,]
boswarmall <- subset(bos, gdd>200)
for(j in unique(boswarmall$year)){
    subby <- subset(boswarmall, year==j)
    minday <- subby[which(subby$doy==min(subby$doy, na.rm=TRUE)),]
    boswarm <- rbind(boswarm, minday)
   }
boswarm <- boswarm[-1,]

boswarm$when <- NA
bosuse <- addwhen(boswarm)
bosusesm <- subset(bosuse, year<1976 | year>1979)

nrow(subset(bosusesm, when=="before 1980"))
nrow(subset(bosusesm, when=="after 1980"))
summary(lm(doy~when, bosusesm))

hist(bosusesm$doy, breaks=30)
plot(density(bosusesm$doy))

boswhen <- ggplot(bosusesm, aes(x=doy, fill=when)) + geom_density(alpha=0.4)
ggplot(bosusesm, aes(x=doy, fill=when)) + geom_histogram(alpha=0.4, aes(y = ..density..), position = 'identity')
}



if(runButte){
## Get the data
butte <- read.delim("climdata/crestedbuttesnowdepth_ed.txt", sep="\t", header=FALSE)
names(butte) <- c("year", "month", "day", "snowdepth") # I assume in mm

butte$date <- as.Date(paste(butte$year, butte$month, butte$day, sep="-"), format="%Y-%m-%d")
butte$doy <- as.numeric(format(butte$date, "%j"))

plot(snowdepth~doy, data=butte) # maybe do max snowdepth in February-April?

## Max snow in March or mean...
buttesnow <- butte[1,]
buttesnow$mean <- NA
buttesnowy <- subset(butte, month==3)
for(j in unique(buttesnowy$year)){
    subby <- subset(buttesnowy, year==j)
    maxday <- subby[which(subby$snowdepth==max(subby$snowdepth, na.rm=TRUE)),]
    maxday <- maxday[1,]
    maxday$mean <- mean(subby$snowdepth, na.rm=TRUE)
    buttesnow <- rbind(buttesnow, maxday)
    
   }
buttesnow <- buttesnow[-1,]

buttesnow$when <- NA
butteuse <- addwhen(buttesnow)
nrow(subset(butteuse, when=="before 1980"))
nrow(subset(butteuse, when=="after 1980"))
summary(lm(snowdepth~when, butteuse))
summary(lm(mean~when, butteuse))

hist(butteuse$snowdepth, breaks=30)
plot(density(butteuse$snowdepth))

ggplot(butteuse, aes(x=snowdepth, fill=when)) + geom_density(alpha=0.4)
ggplot(butteuse, aes(x=snowdepth, fill=when)) + geom_histogram(alpha=0.4, aes(y = ..density..), position = 'identity')

# Get only 40 first years
butteusesm <- subset(butteuse, year<1954 | year>1979)
nrow(subset(butteusesm, when=="before 1980"))
nrow(subset(butteusesm, when=="after 1980"))
buttemaxsnow <- ggplot(butteusesm, aes(x=snowdepth, fill=when)) + geom_density(alpha=0.4)
ggplot(butteusesm, aes(x=mean, fill=when)) + geom_density(alpha=0.4)

summary(lm(snowdepth~when, butteusesm))
summary(lm(mean~when, butteusesm))


## Now do date of no snow in spring ...
if(FALSE){ # takes first no snow day, better (below, not commented out) takes first day of no snow, followed by 9 days of no snow
buttenowsnowday <- butte[1,]
buttenowsnowdaydat <- subset(butte, month>3 & month<8 & snowdepth<5)
for(j in unique(buttenowsnowdaydat$year)){
    subby <- subset(buttenowsnowdaydat, year==j)
    minday <- subby[which(subby$doy==min(subby$doy, na.rm=TRUE)),]
    buttenowsnowday <- rbind(buttenowsnowday, minday)
   }
buttenowsnowdayorig <- buttenowsnowday[-1,]
}
    
buttenowsnowday <- butte[1,]
buttenowsnowdaydat <- subset(butte, month>2 & month<8)
for(j in unique(buttenowsnowdaydat$year)){ # j <- 1981
    subby <- subset(buttenowsnowdaydat, year==j)
    subby$snow10d <- NA
    for(i in c(10:nrow(subby))){
        subby$snow10d[(i-9)] <- sum(subby$snowdepth[(i-9):i])
    }
    subbyzero <- subset(subby, snow10d==0)
    subbyzero$snow10d <- NULL
    minday <- subbyzero[which(subbyzero$doy==min(subbyzero$doy, na.rm=TRUE)),]
    buttenowsnowday <- rbind(buttenowsnowday, minday)
}
buttenowsnowday <- buttenowsnowday[-1,]
buttenowsnowday <- subset(buttenowsnowday, year!=1909) # first year is incomplete data!


buttenowsnowday$when <- NA
butteuse2 <- addwhen(buttenowsnowday)
nrow(subset(butteuse2, when=="before 1980"))
nrow(subset(butteuse2, when=="after 1980"))
summary(lm(doy~when, butteuse2))

hist(butteuse2$doy, breaks=30)
plot(density(butteuse2$doy))

ggplot(butteuse2, aes(x=doy, fill=when)) + geom_density(alpha=0.4)
ggplot(butteuse2, aes(x=doy, fill=when)) + geom_histogram(alpha=0.4, aes(y = ..density..), position = 'identity')

# Get only 40 first years
butteuse2sm <- subset(butteuse2, year<1951 | year>1979)
nrow(subset(butteuse2sm, when=="before 1980"))
nrow(subset(butteuse2sm, when=="after 1980"))
buttewhen40yrs <- ggplot(butteuse2sm, aes(x=doy, fill=when)) + geom_density(alpha=0.4)

summary(lm(doy~when, butteuse2sm))

plot(snowdepth~date, subset(butte, year==1990)) # dates look accurate ...
subset(butteuse2sm, year==1990)
    
# some tricky years ... 
plot(snowdepth~date, subset(butte, year==1934))
plot(snowdepth~date, subset(butte, year==1981))
subset(butteuse2sm, year==1934|year==1981) 
}

if(runMunich==TRUE & runBoston==TRUE & runButte==TRUE){
print("I think you need to run the below separately for some reason.")
require(cowplot)
pdf(paste("graphs/otherdat/climdata.pdf", sep=""), width = 16, height = 10)
plot_grid(boswhen, munwhen, buttewhen40yrs, buttemaxsnow,
    labels = c(' (a) Boston: Day of 200 GDD', '(b) Munich: Day of 200 GDD', '(c) Butte: Snowfree day', '(d) Butte: Max snow depth in March'), ncol=2)
dev.off()
}
