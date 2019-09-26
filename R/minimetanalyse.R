### Started 19 August 2019 ###
### By Lizzie ###

## Reviewing mini meta-analysis of tracking x traits papers ##

######################
### To do items!!! ###
######################
# 
######################

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## set your working directory
setwd("~/Documents/git/projects/temporalvar/R")

###########################################
## Examine accept/reject #s and reasons ##
###########################################

## get part of the data
dall <- read.csv("minimeta/accept_or_reject.csv", header=TRUE)
dall[220:235,]
d <- dall[1:231,]

# trends over time
timetrends <- aggregate(d["Author"], d["Publishing.Year"], FUN=length)
names(timetrends)[names(timetrends)=="Author"] <- "n"

pdf("graphs/otherdat/papersovertime.pdf", width=7.5, height=5)
plot(n~Publishing.Year, data=timetrends)
lines(n~Publishing.Year, data=timetrends)
dev.off()

sum(timetrends$n)
recentpapers <- subset(timetrends, Publishing.Year>2010)
sum(recentpapers$n)

# Why were things rejected?
table(d$Accept.Reject) # this has 68 a, but below shows 69 -- WHY?
table(d$why)

# From the above I combined:
# model + theoretical model + theory/model - no data
# no phenology change measured + no trait or phenological change measured + no treat or phenological change measured
# meta-analysis + review - no data + no data
table(d$accept.reject_2) # this round is single species studies
table(d$why_3) # these are added to above....
# Traits not associated with phenological shifts - focus on herbivory ... need more info from Kelley
46+45+32+10+8+6+2 # studies removed

# Above was just quick look-see, now I try to extract info I want
whogotin <- data.frame(why1=d$why, why2=d$accept.reject_2, why3=d$why_3)
whogotin$why2[whogotin$why2=="a"] <- ""
whogotin$why2[whogotin$why2=="r"] <- "single species"

# Combine these 3 columns ...
whogotin$allreasons <- paste(whogotin$why1, whogotin$why2, whogotin$why3, sep="")
# Great, and now group answers
whogotin$allreasons.simple <- whogotin$allreasons
sort(unique(whogotin$allreasons.simple))

# No phenology data 
whogotin$allreasons.simple[whogotin$allreasons.simple=="Did not measure any phenological tracking"] <- "no phentracking"
whogotin$allreasons.simple[whogotin$allreasons.simple=="No phenological shift data"] <- "no phentracking"
whogotin$allreasons.simple[whogotin$allreasons.simple=="no phenology change measured"] <- "no phentracking"
whogotin$allreasons.simple[whogotin$allreasons.simple=="no measured phenology"] <- "no phentracking"
whogotin$allreasons.simple[whogotin$allreasons.simple=="no phenology data"] <- "no phentracking"
whogotin$allreasons.simple[whogotin$allreasons.simple=="no phenological tracking data"] <- "no phentracking"
# Note below in caption -- no traits or phen is under phen *AND* we excluded one grazing study
whogotin$allreasons.simple[whogotin$allreasons.simple=="no trait or phenological change measured"] <- "no phentracking"
whogotin$allreasons.simple[whogotin$allreasons.simple=="no trat or phenological change measured"] <- "no phentracking"
whogotin$allreasons.simple[whogotin$allreasons.simple=="no phenological shift"] <- "no phentracking"
whogotin$allreasons.simple[whogotin$allreasons.simple=="no phenological tracking"] <- "no phentracking"
whogotin$allreasons.simple[whogotin$allreasons.simple=="No phenological tracking"] <- "no phentracking"
whogotin$allreasons.simple[whogotin$allreasons.simple=="No phenological traking data (artificially forced delayed phenology instead of looking for a response to and environmental variable)"] <- "no phentracking"


whogotin$allreasons.simple[whogotin$allreasons.simple=="Phenological response to grazing treatmentand not to a natural phenomena"] <- "no phentracking"
# Theory/model papers (no data)
whogotin$allreasons.simple[whogotin$allreasons.simple=="Just a model, no data"] <- "model/theory"
whogotin$allreasons.simple[whogotin$allreasons.simple=="model"] <- "model/theory"
whogotin$allreasons.simple[whogotin$allreasons.simple=="theoretical model"] <- "model/theory"
whogotin$allreasons.simple[whogotin$allreasons.simple=="theory/model - no data"] <- "model/theory"
# Reviews
whogotin$allreasons.simple[whogotin$allreasons.simple=="meta-analysis"] <- "review"
whogotin$allreasons.simple[whogotin$allreasons.simple=="review - no data"] <- "review"
# No traits
whogotin$allreasons.simple[whogotin$allreasons.simple=="no trait data"] <- "no trait measured or analysed"
whogotin$allreasons.simple[whogotin$allreasons.simple=="No trait data or mention"] <- "no trait measured or analysed"
whogotin$allreasons.simple[whogotin$allreasons.simple=="no traitanalysis"] <- "no trait measured or analysed"


whogotin$allreasons.simple[whogotin$allreasons.simple=="no trait  analysis"] <- "no trait measured or analysed"
whogotin$allreasons.simple[whogotin$allreasons.simple=="Doesn\xd5t have traits, just geographic origin"] <- "no trait measured or analysed"
whogotin$allreasons.simple[whogotin$allreasons.simple=="Traits not associated with phenological shifts - focus on herbivory"] <- "no trait measured or analysed"
whogotin$allreasons.simple[whogotin$allreasons.simple=="No clear trait to extract"] <- "no trait measured or analysed"
whogotin$allreasons.simple[whogotin$allreasons.simple=="No trait"] <- "no trait measured or analysed"
whogotin$allreasons.simple[whogotin$allreasons.simple=="no trait measured"] <- "no trait measured or analysed"
# Single species
whogotin$allreasons.simple[whogotin$allreasons.simple=="one species - maize"] <- "single species"
whogotin$allreasons.simple[grep("1 species", whogotin$allreasons.simple)] <- "single species"

table(whogotin$allreasons.simple)

# Questions for Kelley:
# (1) What does no data mean (compared to phenological or not trait data)?
# (2) Any differences between 'model' versus 'theory/model - no data' or 'theoretical model'?

# So many studies excluded!
nrow(whogotin)-33 # or 27? Seems to be that below ...


#############################
## Examine extracted data! ##
#############################
datall <- read.csv("minimeta/phentrack_meta.csv", header=TRUE)
dat <- datall[1:207,]
dim(dat)

length(unique(dat$paperID)) # 28 studies (or 27?)
unique(dat$study_level) # safety-check

table(dat$taxongroup_studied) # at least 22 are angiosperms, 12 are lepidoptera, missing seem to be plants also
plantshere <- subset(dat, taxongroup_studied=="angiosperms"|taxongroup_studied=="grasses and forbs"|
    taxongroup_studied=="grasses, forbs, and legumes"|taxongroup_studied=="poaceae")
unique(plantshere$paperID) # 17 papers
butterflies <- subset(dat, taxongroup_studied=="Lepidoptera" | taxongroup_studied=="lepidopterans")
unique(butterflies$paperID) # 4 paper
birdz <- subset(dat, taxongroup_studied=="passerine birds")
unique(birdz$paperID) # 3 papers

linkedstudiez <- subset(dat, link_trackingandtrait_yesno=="yes")
unique(linkedstudiez$paperID)

plantsherelinked <- subset(plantshere, link_trackingandtrait_yesno=="yes")
butterflieslinked <- subset(butterflies, link_trackingandtrait_yesno=="yes")
birdzlinked <- subset(birdz, link_trackingandtrait_yesno=="yes")

unique(plantsherelinked$paperID)
unique(butterflieslinked$paperID)
unique(birdzlinked$paperID)


table(dat$phenophase) # mostly flowering (generally early) or 10th percentile collection date or a little fruiting and leadout
table(dat$track_what) # temperature and precipitation

# Clean up traits ...
table(dat$trait)
dat$trait.simple <- dat$trait
dat$trait.simple[dat$trait.simple=="early/late flowering"] <- "early/late plant phenophase"
dat$trait.simple[dat$trait.simple=="early/late"] <- "early/late plant phenophase"
dat$trait.simple[dat$trait.simple=="earyness"] <- "early/late plant phenophase"
dat$trait.simple[dat$trait.simple=="Host breadth"] <- "niche breadth"
dat$trait.simple[dat$trait.simple=="Number all host plants"] <- "niche breadth"
dat$trait.simple[dat$trait.simple=="Niche breadth score"] <- "niche breadth"
dat$trait.simple[dat$trait.simple=="Number core host plants"] <- "niche breadth"
dat$trait.simple[dat$trait.simple=="polination syndrome (wind or animal)"] <- "wind/insect polination"
# dat$trait.simple[dat$trait.simple==""] <- ""
table(dat$trait.simple)

traits1 <- subset(dat, trait.simple=="early/late plant phenophase" | trait.simple=="rooting depth & earlyness")
length(unique(traits1$paperID))

traits2 <- subset(dat, trait.simple=="nativeness")
length(unique(traits2$paperID))

traits3 <- subset(dat, trait.simple=="early/late migrating" | trait.simple=="Mean of min and max forewing span"|
   trait.simple=="Mobility score" | trait.simple=="Max number generations"| trait.simple=="niche breadth" |
   trait.simple=="overwintering stage")
length(unique(traits3$paperID))

traits4 <- subset(dat, trait.simple=="wind/insect polination")
length(unique(traits4$paperID))

#
dat$link_trackingandtrait_rsq
median(dat$link_trackingandtrait_rsq, na.rm=TRUE) # Damn, that is low

# table we eventually want to build:
supptable <- data.frame(taxa=character(), trait=character(), phenophase=character(),
    nstudiespapers=numeric(), linked=numeric(), notlinked=numeric(), nottested=numeric())
