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
plot(n~Publishing.Year, data=timetrends, xlab="Publication year", ylab="n papers")
lines(n~Publishing.Year, data=timetrends)
dev.off()

sum(timetrends$n)
recentpapers <- subset(timetrends, Publishing.Year>2010)
sum(recentpapers$n)

# Why were things rejected?
table(d$Accept.Reject) 
table(d$why)

# Now I try to extract info I want
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
subset(d, why=="" & accept.reject_2=="" & why_3=="")
# 'Check, not all data possibly entered' was entered



#############################
## Examine extracted data! ##
#############################
datall <- read.csv("minimeta/phentrack_meta.csv", header=TRUE)
dat <- datall[1:210,1:32]
dim(dat)

length(unique(dat$paperID)) # 30 studies 
unique(dat$study_level) # safety-check

## clean up taxa
table(dat$taxongroup_studied)
dat$taxaclean <- dat$taxongroup_studied
dat$taxaclean[grep("epidopter", dat$taxongroup_studied)] <- "Lepidoptera"
dat$taxaclean[grep("grasses", dat$taxongroup_studied)] <- "plants"
dat$taxaclean[grep("poaceae", dat$taxongroup_studied)] <- "plants"
dat$taxaclean[grep("angiosperms", dat$taxongroup_studied)] <- "plants"
dat$taxaclean[grep("warblers", dat$taxongroup_studied)] <- "passerine birds"

table(dat$taxaclean)
# Papers per taxa grouping... 
plantshere <- subset(dat, taxaclean=="plants")
unique(plantshere$paperID) # 20 papers!
butterflies <- subset(dat, taxaclean=="Lepidoptera")
unique(butterflies$paperID) # 4 papers
birdz <- subset(dat, taxaclean=="passerine birds")
unique(birdz$paperID) # 4 papers
aphidz <- subset(dat, taxaclean=="Aphidoidea")
unique(aphidz$paperID) # 1 paper
plankton <- subset(dat, taxaclean=="plankton")
unique(plankton$paperID) # 1 paper


# linking traits and tracking ....
table(dat$link_trackingandtrait_yesno)
linkedstudiez <- subset(dat, link_trackingandtrait_yesno=="yes")
unique(linkedstudiez$paperID)

plantsherelinked <- subset(plantshere, link_trackingandtrait_yesno=="yes")
butterflieslinked <- subset(butterflies, link_trackingandtrait_yesno=="yes")
birdzlinked <- subset(birdz, link_trackingandtrait_yesno=="yes")
aphizlinked <- subset(aphidz, link_trackingandtrait_yesno=="yes")
planktonlinked <- subset(plankton, link_trackingandtrait_yesno=="yes")

unique(plantsherelinked$paperID)
unique(butterflieslinked$paperID)
unique(birdzlinked$paperID)
unique(aphizlinked$paperID) # didnottry
unique(planktonlinked$paperID) # didnottry

# clean up phenophases
table(dat$phenophase) 
dat$phenophase.simple <- dat$phenophase
dat$phenophase.simple[dat$phenophase=="10th percentile collection date"] <- "appearance/collection date"
dat$phenophase.simple[dat$phenophase=="date of first appearance"] <- "appearance/collection date"
dat$phenophase.simple[dat$phenophase=="first emergance date"] <- "appearance/collection date"
dat$phenophase.simple[dat$phenophase=="first appearance"] <- "appearance/collection date"
dat$phenophase.simple[dat$phenophase=="flight season timing"] <- "flight timing"
dat$phenophase.simple[dat$phenophase=="first flight"] <- "flight timing"
dat$phenophase.simple[dat$phenophase=="spring passing"] <- "flight timing"

dat$phenophase.simple[dat$phenophase=="leaf emergance"] <- "budbreak/leafing"
dat$phenophase.simple[dat$phenophase=="leafout"] <- "budbreak/leafing"
dat$phenophase.simple[dat$phenophase=="leaf unfolding"] <- "budbreak/leafing"
dat$phenophase.simple[dat$phenophase=="bud break"] <- "budbreak/leafing"
dat$phenophase.simple[dat$phenophase=="budbreak"] <- "budbreak/leafing"
dat$phenophase.simple[dat$phenophase=="timing of flowering"] <- "flowering"
dat$phenophase.simple[dat$phenophase=="days to first flower"] <- "flowering"
dat$phenophase.simple[dat$phenophase=="days to last flower"] <- "flowering"
dat$phenophase.simple[dat$phenophase=="first flowering"] <- "flowering"
dat$phenophase.simple[dat$phenophase=="flowering"] <- "flowering"
dat$phenophase.simple[dat$phenophase=="peak flowering"] <- "flowering"
dat$phenophase.simple[dat$phenophase=="flowering onset"] <- "flowering"
dat$phenophase.simple[dat$phenophase=="flowering time"] <- "flowering"
dat$phenophase.simple[dat$phenophase=="flowering timing"] <- "flowering"
dat$phenophase.simple[dat$phenophase=="first flowering date"] <- "flowering"

dat$phenophase.simple[dat$phenophase=="fruting"] <- "fruiting"
dat$phenophase.simple[dat$phenophase=="fruiting time"] <- "fruiting"
dat$phenophase.simple[dat$phenophase=="days to first fruit"] <- "fruiting"
dat$phenophase.simple[dat$phenophase=="timing of fruit coloring"] <- "fruiting"
dat$phenophase.simple[dat$phenophase=="flowering period"] <- "flowering length"
dat$phenophase.simple[dat$phenophase=="flowering duration"] <- "flowering length"
dat$phenophase.simple[dat$phenophase=="flowring duration"] <- "flowering length"
dat$phenophase.simple[dat$phenophase=="lay date"] <- "breeding time"
table(dat$phenophase.simple)

# clean up trackwhat ...
table(dat$track_what) 
dat$trackwhat.simple <- dat$track_what
dat$trackwhat.simple[grep("temperature", dat$track_what)] <- "temperature"
dat$trackwhat.simple[grep("NAO", dat$track_what)] <- "climate mode"
dat$trackwhat.simple[dat$track_what=="monthly rainfall"] <- "precipitation"
table(dat$trackwhat.simple)
# Munson2017 is the missing one and it's because they did 15 climate metrics mainly temp and precip
dat$trackwhat.simple[dat$paperID=="Munson2017"][1] <- "temperature"
dat$trackwhat.simple[dat$paperID=="Munson2017"][2] <- "precipitation"
table(dat$trackwhat.simple)
100*(table(dat$trackwhat.simple)/sum(table(dat$trackwhat.simple))) # percents

# how many measured multiple climate metrics?
multiclimmetrics1 <- aggregate(dat["track_what"], dat[c("paperID", "trackwhat.simple")], FUN=length)
multiclimmetrics <- aggregate(multiclimmetrics1["trackwhat.simple"], multiclimmetrics1["paperID"],
    FUN=length)
subset(multiclimmetrics, trackwhat.simple>1) # 4 paper measured 2 metrics

## Before we clean traits and phenophases ...
## How many did each study measure? (Not exact since for ones with MANY values we didn't always enter all of them)
# f(x) for mode
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
# phenophase
phenmetrics1 <- aggregate(dat["track_what"], dat[c("paperID", "phenophase")], FUN=length)
phenmetrics <- aggregate(phenmetrics1["phenophase"], phenmetrics1["paperID"],
    FUN=length)
getmode(phenmetrics$phenophase)
median(phenmetrics$phenophase)
# traits
traitmetrics1 <- aggregate(dat["track_what"], dat[c("paperID", "trait")], FUN=length)
traitmetrics <- aggregate(traitmetrics1["trait"], traitmetrics1["paperID"],
    FUN=length)
getmode(traitmetrics$trait)
median(traitmetrics$trait)

#####
# Clean up traits ...
# Note that minnimum spring temperature is '10th quantile of minimum temperature in native range'
table(dat$trait)
dat$trait.simple <- dat$trait
dat$trait.simple[grep("root", dat$trait.simple)] <- "root traits"
dat$trait.simple[dat$trait.simple=="number of tillers"] <- "root traits"
table(dat$trait)

dat$trait.simple[dat$trait.simple=="early/late flowering"] <- "early/late phenophase"
dat$trait.simple[dat$trait.simple=="early/late"] <- "early/late phenophase"
dat$trait.simple[dat$trait.simple=="earyness"] <- "early/late phenophase"
dat$trait.simple[dat$trait.simple=="earliness"] <- "early/late phenophase"
dat$trait.simple[dat$trait.simple=="early/late migrating"] <- "early/late phenophase"
dat$trait.simple[dat$trait.simple=="timing of flight season"] <- "early/late phenophase"
dat$trait.simple[dat$trait.simple=="breeding season"] <- "early/late phenophase" # CHECK!
dat$trait.simple[dat$trait.simple=="baseline date first appearance"] <- "early/late phenophase"

dat$trait.simple[(dat$trait.simple=="percentage national 10-km grid cells occupied X U.K. latitudinal extent")] <- "range traits"
dat$trait.simple[(dat$trait.simple=="minnimum spring temperature")] <- "range traits"

dat$trait.simple[grep("bread", dat$trait.simple)] <- "niche breadth"
dat$trait.simple[dat$trait.simple=="Host breadth"] <- "niche breadth"
dat$trait.simple[dat$trait.simple=="Number all host plants"] <- "niche breadth"
dat$trait.simple[dat$trait.simple=="Niche breadth score"] <- "niche breadth"
dat$trait.simple[dat$trait.simple=="Number core host plants"] <- "niche breadth"
dat$trait.simple[dat$trait.simple=="larval host breadth"] <- "niche breadth"
dat$trait.simple[grep("host plants", dat$trait.simple)] <- "niche breadth"

dat$trait.simple[grep("obilit", dat$trait.simple)] <- "mobility"
dat$trait.simple[grep("dispresal ability", dat$trait.simple)] <- "mobility"

dat$trait.simple[grep("wing", dat$trait.simple)] <- "wing size"

dat$trait.simple[grep("growth speed", dat$trait.simple)] <- "growth speed (plankton)"

dat$trait.simple[grep("leaf size", dat$trait.simple)] <- "leaf/shoot size"
dat$trait.simple[grep("leaf area", dat$trait.simple)] <- "leaf/shoot size"
dat$trait.simple[grep("shoot length", dat$trait.simple)] <- "leaf/shoot size"

dat$trait.simple[grep("leaf long", dat$trait.simple)] <- "leaf longevity"
dat$trait.simple[grep("leaf nitrogen", dat$trait.simple)] <- "leaf nitrogen"
dat$trait.simple[grep("woody", dat$trait.simple)] <- "woody/herbaceous"
dat$trait.simple[grep("overwinter", dat$trait.simple)] <- "overwintering"

dat$trait.simple[dat$trait.simple=="photosynthetic rate"] <- "other leaf traits"
dat$trait.simple[dat$trait.simple=="C3/C4"] <- "other leaf traits"
dat$trait.simple[dat$trait.simple=="shade tolerance"] <- "other leaf traits"
dat$trait.simple[dat$trait.simple=="wood anatomy"] <- "other leaf traits"
dat$trait.simple[dat$trait.simple=="phylogenetic group"] <- "other leaf traits"
dat$trait.simple[dat$trait.simple=="wind/insect pollination"] <- "other leaf traits"
dat$trait.simple[dat$trait.simple=="life form"] <- "other leaf traits"
dat$trait.simple[dat$trait.simple=="annual/perrenial"] <- "other leaf traits"
dat$trait.simple[dat$trait.simple=="canopy width"] <- "other leaf traits"
dat$trait.simple[dat$trait.simple=="fruit type"] <- "other leaf traits"
dat$trait.simple[dat$trait.simple=="growth habit"] <- "other leaf traits"

dat$trait.simple[dat$trait.simple=="leaf water & nitrogen content"] <- "other leaf traits"
dat$trait.simple[dat$trait.simple=="leaf toughness"] <- "other leaf traits"
dat$trait.simple[dat$trait.simple=="leaf tannin content"] <- "other leaf traits"

dat$trait.simple[dat$trait.simple=="# flowers"] <- "number of flowers/fruits"
dat$trait.simple[dat$trait.simple=="# fruits"] <- "number of flowers/fruits"
dat$trait.simple[grep("seed", dat$trait.simple)] <- "offspring traits (weight/size/number)"
dat$trait.simple[grep("pods", dat$trait.simple)] <- "offspring traits (weight/size/number)"
dat$trait.simple[dat$trait.simple=="number of broods"] <- "offspring traits (weight/size/number)"
dat$trait.simple[dat$trait.simple=="number of broods"] <- "offspring traits (weight/size/number)"

dat$trait.simple[grep("polin", dat$trait.simple)] <- "seed weight/size/number"

dat$trait.simple[dat$trait.simple=="deciduois/evergreen"] <- "other plant traits"
dat$trait.simple[dat$trait.simple=="Deciduous/evergreen"] <- "other plant traits"
dat$trait.simple[dat$trait.simple=="flower longevity"] <- "other plant traits"
dat$trait.simple[dat$trait.simple=="flower type"] <- "other plant traits"
dat$trait.simple[dat$trait.simple=="stomatal conductance"] <- "other plant traits"
dat$trait.simple[dat$trait.simple=="species type"] <- "other plant traits"
dat$trait.simple[dat$trait.simple=="mean isotherm"] <- "other plant traits"
dat$trait.simple[dat$trait.simple=="mean spring variability"] <- "other plant traits"
dat$trait.simple[dat$trait.simple=="polination syndrome (wind or animal)"] <- "other plant traits"
dat$trait.simple[dat$trait.simple=="monoecious/heteroecious"] <- "other plant traits"


dat$trait.simple[dat$trait.simple=="mass"] <- "other bird traits"
dat$trait.simple[dat$trait.simple=="annual brood number"] <- "other bird traits"
dat$trait.simple[dat$trait.simple=="brain mass"] <- "other bird traits"

dat$trait.simple[dat$trait.simple=="food type"] <- "other Lepidopteran traits"
dat$trait.simple[dat$trait.simple=="average length of flight season"] <- "other Lepidopteran traits"
dat$trait.simple[dat$trait.simple=="wing size"] <- "other Lepidopteran traits"
dat$trait.simple[dat$trait.simple=="voltinism"] <- "other Lepidopteran traits"
dat$trait.simple[dat$trait.simple=="Max number generations"] <- "other Lepidopteran traits"
dat$trait.simple[dat$trait.simple=="average number of generations"] <- "other Lepidopteran traits"

dat$trait.simple[dat$trait.simple=="diet"] <- "diet traits"
dat$trait.simple[dat$trait.simple=="diet breadth"] <- "diet traits"
dat$trait.simple[dat$trait.simple=="food specialization"] <- "diet traits"
dat$trait.simple[grep("abitat", dat$trait.simple)] <- "habitat traits"
dat$trait.simple[grep("migrat*", dat$trait.simple)] <- "migration traits"
dat$trait.simple[grep("range", dat$trait.simple)] <- "range traits"

table(dat$trait.simple)

# Simplify even more ...
dat$phenophase.supersimple <- dat$phenophase.simple
dat$phenophase.supersimple[dat$phenophase.simple=="peak flowering"] <- "flowering/fruiting"
dat$phenophase.supersimple[grep("fruit", dat$phenophase.simple)] <- "flowering/fruiting"
dat$phenophase.supersimple[grep("flower", dat$phenophase.simple)]  <- "flowering/fruiting"
dat$phenophase.supersimple[dat$phenophase.simple=="median emergance date"]  <- "last/median emergence dates"
dat$phenophase.supersimple[dat$phenophase.simple=="last emergance date"]  <- "last/median emergence dates"

table(dat$phenophase.supersimple)

# try to summarize
dattouse <- subset(dat, link_trackingandtrait_yesno=="yes" | link_trackingandtrait_yesno=="no")
table(dattouse$trackwhat.simple)

didnottrydat <- subset(dat, link_trackingandtrait_yesno=="didnottry")
unique(didnottrydat$paperID)
length(unique(didnottrydat$paperID))
length(unique(dattouse$paperID)) # good! These numbers add up.
intersect(didnottrydat$paperID, dattouse$paperID) # two studies linked some traits, but did not try on others

datlinked <- subset(dat, link_trackingandtrait_yesno=="yes")
datforsummlinked <- subset(datlinked, select=c("paperID", "taxaclean", "phenophase.supersimple",
    "trait.simple"))
datforsummlinked.nodups <- datforsummlinked[!duplicated(datforsummlinked), ]

summlinked <- aggregate(datforsummlinked.nodups["paperID"], datforsummlinked.nodups[c("taxaclean",
    "phenophase.supersimple", "trait.simple")], FUN=length)
names(summlinked)[names(summlinked)=="paperID"] <- "n papers linked"

datnotlinked <- subset(dat, link_trackingandtrait_yesno=="no")
datforsummnotlinked <- subset(datnotlinked, select=c("paperID", "taxaclean", "phenophase.supersimple",
    "trait.simple"))
datforsummnotlinked.nodups <- datforsummnotlinked[!duplicated(datforsummnotlinked), ]

summnotlinked <- aggregate(datforsummnotlinked.nodups["paperID"], datforsummnotlinked.nodups[c("taxaclean",
    "phenophase.supersimple", "trait.simple")], FUN=length)
names(summnotlinked)[names(summnotlinked)=="paperID"] <- "n papers not linked"

tableforpaper <- merge(summlinked, summnotlinked, by=c("taxaclean", "phenophase.supersimple",
    "trait.simple"), all.x=TRUE, all.y=TRUE)

sum(tableforpaper[,4], na.rm=TRUE)
sum(tableforpaper[,5], na.rm=TRUE)


tableforpaper <- tableforpaper[-c(1:2),]
names(tableforpaper) <- c("Taxa", "Phenophase", "Trait", "n linked", "n not linked")

# What traits are most common? And are they linked?
traitz <- aggregate(tableforpaper[c("n linked", "n not linked")],
    tableforpaper["Trait"], FUN=sum, na.rm=TRUE)
# traitz <- traitz[sort(traitz$Trait),]
checkearlylate <- subset(tableforpaper, Trait=="early/late phenophase")
sum(checkearlylate[,4], na.rm=TRUE)
sum(checkearlylate[,5], na.rm=TRUE)
# which studies?
checkstudies2 <- dat[which(dat$phenophase.simple %in% checkearlylate$Phenophase),]
checkstudies1 <-  checkstudies2[which(checkstudies2$trait.simple %in% checkearlylate$Trait),]
checkstudies <- subset(checkstudies1, link_trackingandtrait_yesno=="yes")
##

## Look at some specific plant traits ...
native <- dat[grep("nativeness", dat$trait.simple),] # basically one study found a link and one did not 
otherleaf <- dat[grep("other leaf", dat$trait.simple),] # mixed!

rootz1 <- dat[grep("root", dat$trait),]
rootz2 <- dat[grep("tiller", dat$trait),] # this study did not try to link traits and tracking 
rootz <- rbind(rootz1, rootz2)

# high SLA is thin leaves that draw down resources quickly but do not last long ... 
sla <- dat[grep("SLA", dat$trait.simple),] # but wait, all the SLA studies did not try to link ...
