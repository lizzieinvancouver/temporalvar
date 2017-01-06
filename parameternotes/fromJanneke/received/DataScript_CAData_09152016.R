### Data Exploration for Lizzie
## Exploring R*, relative abundance, seed production for 11 California Annuals

## Set home directory
setwd("C:/Users/janneke/Dropbox/Research/CA_projects/Midland_main/Lizzie")

## Read in data
PlotDat <- read.csv("PlotDat.csv", header = TRUE)
SpeciesDat <- read.csv("SpeciesDat.csv", header = TRUE)
InflDat <- read.csv("InflDat.csv", header=TRUE)
SeedDat <- read.csv("SeedDat.csv", header=TRUE)

### Examine correlation between R* and relative abundance in competition, ungrazed plots
#Extract R* per species
rstarH20 <- tapply(PlotDat[,4], as.factor(PlotDat[,3]), mean)
rstarPAR <- tapply(PlotDat[,5], as.factor(PlotDat[,3]), mean)
rstarP <- tapply(PlotDat[,6], as.factor(PlotDat[,3]), mean)
rstarN <- tapply(PlotDat[,7], as.factor(PlotDat[,3]), mean)
Rstarall <- cbind(rstarH20,rstarPAR, rstarP, rstarN)

#Remove Rstar values for mixtures, bare plots, and one species that did not germinate
Rstarall <- Rstarall[-c(2,4,8,11,12),]

#Calculate relative abundance of each species in mixture plots (Above ground biomass)
MixtureDat <- SpeciesDat[SpeciesDat[,3]=="All",]
totbio <- tapply(MixtureDat[,5],as.factor(MixtureDat[,2]),sum)
MixtureDat <- merge(MixtureDat, totbio, by.x="Plot", by.y="row.names")
relabun <- MixtureDat[,5]/MixtureDat[,8]
MixtureDat <- cbind(MixtureDat, relabun)
dimnames(MixtureDat)[[2]][8:9] <- c("Total Biomass","RAb")

#Calculate the average relative abundance
mrelabun <- tapply(MixtureDat[,9], as.factor(MixtureDat[,4]), mean)

#Merge
RstarBio <- merge(Rstarall,mrelabun, by="row.names")
dimnames(RstarBio)[[2]][c(1,6)] <- c("Spp","RA")

#plot
X11(width=6, height=6)
par(mfrow=c(2,2), omi=c(0,0,0,0), mai=c(0.5,0.5,0.4,0.4), tck=-0.02, mgp=c(1.2,0.4,0))

for(i in 2:5){
  plot(RstarBio[,i],RstarBio[,6], xlab=dimnames(RstarBio)[[2]][i], cex=1.5, ylab="Relative Abundance", pch=21, bg="grey")
  test <- cor.test(RstarBio[,i],RstarBio[,6], alternative="less")
  mtext(paste("r=",signif(test$estimate,3), sep="", collapse=NULL), side=3, adj=0.95, line=-1, cex=0.75)
  mtext(paste("p=",signif(test$p.value,3), sep="", collapse=NULL), side=3, adj=0.95, line=-2, cex=0.75)
}


##Calucate biomass of aboveground, belowground, seeds in monoculture
#Subset monocultures
MonoDat <- PlotDat[PlotDat[,3]!="All"&PlotDat[,3]!="Exotic"&PlotDat[,3]!="Native"
                      &PlotDat[,3]!="Bare"&PlotDat[,3]!="Mm",] #only monocultures that had significant germination

#Scale aboveground / below ground biomass to m2 basis
MonoDat[,8] <- MonoDat[,8] / (0.5*0.1) #clipstrip of 0.5 x 0.1 meters
MonoDat[,9] <- MonoDat[,9] / (2 * pi * 0.025^2) #2 soil cores, 5 cm diameter

#Add seed set
Sdwght <- rep(NA, times=dim(MonoDat)[1])

for(i in 1:dim(MonoDat)[1]){
  plt <- MonoDat[i,2]; sp <- as.character(MonoDat[i,3])
  inflm2 <- SpeciesDat[SpeciesDat[,2]==plt,6] / SpeciesDat[SpeciesDat[,2]==plt,7] #inflorescence per m2
  inflwht <- mean(InflDat[InflDat[,2]==plt & InflDat[,3]==sp,7])
  if(inflm2==0 &is.na(inflwht)==TRUE){inflwht <- 0}
  if(inflm2>0 &is.na(inflwht)==TRUE){inflwht <- mean(InflDat[InflDat[,3]==sp,7])}
  
  Sdwght[i] <- inflm2*inflwht
  
  if(sp=="Ab"|sp=="Am"|sp=="Cc"|sp=="Cp"){
    seedunitinfl <- mean(na.omit(InflDat[InflDat[,2]==plt & InflDat[,3]==sp,5]))
    seednounit <- mean(na.omit(InflDat[InflDat[,2]==plt & InflDat[,3]==sp,6]))
    if(is.na(seedunitinfl)==TRUE){seedunitinfl <- mean(na.omit(InflDat[InflDat[,3]==sp,5]))}
    if(is.na(seednounit)==TRUE){seednounit <- mean(na.omit(InflDat[InflDat[,3]==sp,6]))}
    seedweight <- SeedDat[SeedDat[,1]==sp,2]
    Sdwght[i] <- inflm2 * seedunitinfl * seednounit * seedweight
  }
}

#Merge
MonoDat <- cbind(MonoDat, Sdwght)

#averages
AbvBio <- tapply(MonoDat[,8], as.factor(MonoDat[,3]), mean); AbvBio <- AbvBio[is.na(AbvBio)==FALSE]
BlwBio <- tapply(MonoDat[is.na(MonoDat[,9])==FALSE,9], as.factor(MonoDat[is.na(MonoDat[,9])==FALSE,3]), mean); BlwBio <- BlwBio[is.na(BlwBio)==FALSE] 
SdBio <- tapply(MonoDat[,10], as.factor(MonoDat[,3]), mean); SdBio <- SdBio[is.na(SdBio)==FALSE]


#Plot of biomass per m2 in monoculture
X11(width=5, height=8)
par(mfrow=c(3,1), omi=c(0,0,0,0), mai=c(0.5,0.5,0.4,0.4), tck=-0.02, mgp=c(1.2,0.4,0))
barplot(AbvBio[order(AbvBio,decreasing=TRUE)], cex.names=0.75, ylab="g / m2")
title("Aboveground Biomass")
barplot(BlwBio[order(BlwBio,decreasing=TRUE)], cex.names=0.75, ylab="g / m2")
title("Belowground Biomass")
barplot(SdBio[order(SdBio,decreasing=TRUE)], cex.names=0.75, ylab="g / m2")
title("Seed Biomass")


##calculate relative abundance as above ground biomass and seed in mixtures
#first calculate total seed biomass per m2
Sdwght <- rep(NA, times=dim(MixtureDat)[1])

for(i in 1:dim(MixtureDat)[1]){
  plt <- MixtureDat[i,1]; sp <- as.character(MixtureDat[i,4])
  inflm2 <- MixtureDat[i,6] / MixtureDat[i,7] #inflorescence per m2
  inflwht <- mean(InflDat[InflDat[,2]==plt & InflDat[,3]==sp,7])
  if(inflm2==0 &is.na(inflwht)==TRUE){inflwht <- 0}
  Sdwght[i] <- inflm2*inflwht
  
  if(sp=="Ab"|sp=="Am"|sp=="Cc"|sp=="Cp"){
    seedunitinfl <- mean(na.omit(InflDat[InflDat[,2]==plt & InflDat[,3]==sp,5]))
    seednounit <- mean(na.omit(InflDat[InflDat[,2]==plt & InflDat[,3]==sp,6]))
       if(is.na(seedunitinfl)==TRUE){seedunitinfl <- mean(na.omit(InflDat[InflDat[,3]==sp,5]))}
       if(is.na(seednounit)==TRUE){seednounit <- mean(na.omit(InflDat[InflDat[,3]==sp,6]))}
    seedweight <- SeedDat[SeedDat[,1]==sp,2]
    Sdwght[i] <- inflm2 * seedunitinfl * seednounit * seedweight
  }
}

#Merge
MixtureDat <- cbind(MixtureDat, Sdwght)
totsbio <- tapply(MixtureDat[,10],as.factor(MixtureDat[,1]),sum)
MixtureDat <- merge(MixtureDat, totsbio, by.x="Plot", by.y="row.names")
relabuns <- MixtureDat[,10]/MixtureDat[,11]
MixtureDat <- cbind(MixtureDat, relabuns)
dimnames(MixtureDat)[[2]][11:12] <- c("Total Seed Biomass", "RAs")

##Average relativ abundance: Biomass & seeds
RAbio <- tapply(MixtureDat[,9], as.factor(MixtureDat[,4]), mean) 
RAseed <- tapply(MixtureDat[,12], as.factor(MixtureDat[,4]), mean)

#Plot of relative abundance in mixture
X11(width=6, height=8)
par(mfrow=c(2,1), omi=c(0,0,0,0), mai=c(0.5,0.5,0.4,0.4), tck=-0.02, mgp=c(1.2,0.4,0))
barplot(RAbio[order(RAbio,decreasing=TRUE)], cex.names=0.75, ylab="Relative Abundance")
title("Aboveground Biomass")
barplot(RAseed[order(RAseed,decreasing=TRUE)], cex.names=0.75, ylab="Relative Abundance")
title("Seed production")
