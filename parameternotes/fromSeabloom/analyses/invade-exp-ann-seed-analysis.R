# READ Annual Grass Seed DATA FROM RECIPROCAL INVASION EXPERIMENT CONDUCTED AT SEDGWICK
# IN SANTA YNEZ CALIFORNIA. 
# There are three experiment which were all done in different parts of the same field
# and share a core treatment set: Starting Community (Annual or Perennial) and 
# Seed Addition (Annual seed added to perennial plots and Perennial seed added 
# to annual plots).
#
# Example publications:
# Seabloom, E. W., W. S. Harpole, O. J. Reichman, and D. Tilman. 2003. 
# Invasion, competitive dominance, and resource use by exotic and native 
# California grassland species. Proceedings of the National Academy of Sciences 
# of the United States of America 100:13384-13389.
#
# Seabloom, E. W., and S. A. Richards. 2003. Multiple stable equilibria in 
# grasslands mediated by herbivore population dynamics and foraging behavior. 
# Ecology 84:2891–2904.
#
# Seabloom, E. W., O. N. Bjornstad, B. M. Bolker, and O. J. Reichman. 2005. 
# The spatial signature of environmental heterogeneity, dispersal, and competition 
# in successional grasslands. Ecological Monographs 75:199-214.
#
# Everard, K., E. W. Seabloom, W. S. Harpole, and C. de Mazancourt. 2010. 
# Plant Water Use Affects Competition for Nitrogen: Why Drought Favors Invasive 
# Species in California. American Naturalist 175:85-97.

## Janueary 2017, edits and calculations by Lizzie ##

# Clear all existing data
rm(list=ls())
# Close graphics devices
graphics.off()
library(plyr)
library(ggplot2)
library(nlme)

## Set home directory
setwd("~/Documents/git/projects/temporalvar/parameternotes/fromSeabloom/analyses")

panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    #r <- abs(cor(x, y))
    r <-(cor(x, y, use='complete.obs'))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    #text(0.5, 0.5, txt, cex = cex * r)
    text(0.5, 0.5, txt, cex = 1.5)
}

# Output plot (quadrat) level data data
df <- read.csv(file = "input/invade-exp-data-output-annual-seeds.csv")

# Log some data
df$plntmass.lg <- log(df$plntmass)
df$seedmass.lg <- log(df$seedmass)
df$sdms.plt.lg <- log(df$sdms.plt)
df$sds.plt.lg <- log(df$sds.plt)
df$sdratio.lg <- log(df$sdratio)

# Make some factors
df$year.f <- as.factor(as.character(df$year))

# Relabel Plant Community
df$plant.txt[df$plant == "A"] <- "Annual"
df$plant.txt[df$plant == "P"] <- "Perennial"

# Plot Raw Data Correlations
pairs(df[c("plntmass", "seedmass", "sds.plt", "sdms.plt", "sdratio")],
 upper.panel=panel.cor, lower.panel=panel.smooth,
 diag.panel=panel.hist, cex.labels = 1.1, font.labels=2)

# Plot logged data correlations
pairs(df[c("plntmass.lg", "seedmass.lg", "sds.plt.lg", "sdms.plt.lg", "sdratio.lg")],
 upper.panel=panel.cor, lower.panel=panel.smooth,
 diag.panel=panel.hist, cex.labels = 1.1, font.labels=2)

# Do some quick plotting of data
# Pull out unseeded annual nutrient plots
nut.df <- with(df, df[exp == "NUTRIENT" & seed == 0 & plant == "A",])

qplot(trt, sds.plt.lg, data=nut.df, xlab="Nutrient Treatment", ylab="log(Seeds per Plant)") + 
# facet_wrap(~taxa, scales="fixed", nrow=1) +
facet_grid(year.f~taxa, scales="fixed") +
geom_boxplot()

qplot(trt, sdms.plt.lg, data=nut.df, xlab="Nutrient Treatment", ylab="log(Seed Mass per Plant)") + 
# facet_wrap(~taxa, scales="fixed", nrow=1) +
facet_grid(year.f~taxa, scales="fixed") +
geom_boxplot()

# Pull out unseeded annual watering plots
wat.df <- with(df, df[exp == "WATER" & seed == 0 & plant == "A",])

qplot(trt, sds.plt.lg, data=wat.df, xlab="Watering Treatment", ylab="log(Seed Mass per Plant)") + 
# facet_wrap(~taxa, scales="fixed", nrow=1) +
facet_grid(year.f~taxa, scales="fixed") +
geom_boxplot()

qplot(trt, sdms.plt.lg, data=wat.df, xlab="Watering Treatment", ylab="log(Seeds per Plant)") + 
# facet_wrap(~taxa, scales="fixed", nrow=1) +
facet_grid(year.f~taxa, scales="fixed") +
geom_boxplot()

# Compare plants in annual and perennial plots with no other treatment
avp.df <- with(df, df[trt == "CONTROL" & seed == 0,])

qplot(plant.txt, sds.plt.lg, data=avp.df, xlab="Plant Community", ylab="log(Seeds per Plant)") + 
facet_grid(year.f~taxa, scales="fixed") +
geom_boxplot()

qplot(plant.txt, sdms.plt.lg, data=avp.df, xlab="Plant Community", ylab="log(Seed Mass per Plant)") + 
facet_grid(year.f~taxa, scales="fixed") +
geom_boxplot()

qplot(plant.txt, seedmass.lg, data=avp.df, xlab="Plant Community", ylab="log(Seed Mass)") + 
facet_grid(year.f~taxa, scales="fixed") +
geom_boxplot()

# Quick regression of nutrient experiment
# Number of Seeds per Plant
lme.sds.nut <- lme(sds.plt.lg~(taxa*trt*plant.txt*seed), 
 random=~1|block/plot/year, 
 data=df[df$exp == "NUTRIENT",], 
 method="REML")
summary(lme.sds.nut)
anova(lme.sds.nut)

# Mass of Seeds per Plant
lme.sdms.nut <- lme(sdms.plt.lg~(taxa*trt*plant.txt*seed), 
 random=~1|block/plot/year, 
 data=df[df$exp == "NUTRIENT",], 
 method="REML")
summary(lme.sdms.nut)
anova(lme.sdms.nut)

# Mass of a Seed
lme.seed.nut <- lme(seedmass.lg~(taxa*trt*plant.txt*seed), 
 random=~1|block/plot/year, 
 data=df[df$exp == "NUTRIENT",], 
 method="REML")
summary(lme.seed.nut)
anova(lme.seed.nut)


## Some plots by Lizzie

# Notes to self:
# leaving in unseeded annual watering plot

# Biomass per individual (Bi)
ggplot(df,
     aes(x=trt, y=plntmass, color=taxa)) +
    # scale_x_discrete(name="Temperature") + scale_y_continuous(name="Days to BB") +
     facet_wrap(~year, nrow=3) + 
     geom_boxplot()

# Per seed biomass (b0)
ggplot(df,
     aes(x=trt, y=seedmass, color=taxa)) +
    # scale_x_discrete(name="Temperature") + scale_y_continuous(name="Days to BB") +
     facet_wrap(~year, nrow=3) + 
     geom_boxplot()

# Biomass to seed conversion (like phi, but without overwintering)


# some calculations ...

## from invade-exp-file-descriptions.txt:
# plntmass	Mean mass of a single plant (g)
# seedmass	Mean mass of a single seed (g)
# sds.plt	Mean number of seeds per plant
# sdms.plt	Mean mass of seeds per plant (g)

# leaving in unseeded annual watering plot
# but remove the seeded plots
df.control <- subset(df, seed==0)
df.control$biomasstoseed <- (df.control$sds.plt/df.control$plntmass)
# sanity check that I understand the data ...
seedmass.calc <- df.control$seedmass * df.control$sds.plt
plot(df.control$sdms.plt~seedmass.calc)

param.summary <-
    ddply(df.control, c("taxa", "trt", "plant.txt", "year"), summarise,
        mean.seedmass = mean(seedmass), sd = sd(seedmass),
          # Mean mass of a single plant (g)
        mean.plntmass = mean(plntmass), sd = sd(plntmass),
          # Mean mass of a single seed (g)
        mean.sds.plt = mean(sds.plt), sd = sd(sds.plt),
          # Mean number of seeds per plant (we don't really need this)
        mean.sdms.plt = mean(sdms.plt), sd = sd(sdms.plt),
          # Mean mass of seeds per plant (g)
        mean.biomasstoseed = mean(biomasstoseed), sd = sd(biomasstoseed))
          # biomass to seed number conversion

# get min max for parameters spreadsheet
param.minmax <-
      ddply(param.summary, c("taxa"), summarise,
      minofmean.seedmass = min(mean.seedmass), maxofmean.seedmass=max(mean.seedmass),
          meanofmean.seedmass=mean(mean.seedmass),
      minofmean.plntmass = min(mean.plntmass), maxofmean.plntmass=max(mean.plntmass),
          meanofmean.plntmass=mean(mean.plntmass),
      minofmean.sds.plt = min(mean.sds.plt), maxofmean.sds.plt=max(mean.sds.plt),
          meanofmean.sds.plt=mean(mean.sds.plt),
      minofmean.mean.sdms.plt = min(mean.sdms.plt),
          maxofmean.mean.sdms.plt=max(mean.sdms.plt),
          meanofmean.sdms.plt=mean(mean.sdms.plt),
      minofmean.biomasstoseed = min(mean.biomasstoseed),
          maxofmean.biomasstoseed=max(mean.biomasstoseed),
          meanofmean.biomasstoseed=mean(mean.biomasstoseed))

### Hmm the biomasstoseed values are now way different than HillRisLambers, need to check both!

# and some plots of those calculations
ggplot(param.summary, aes(x=year,y=mean.seedmass, fill=trt))+
    facet_wrap(~taxa, nrow=3) +
    geom_boxplot()

ggplot(param.summary, aes(x=year,y=mean.seedmass, fill=plant.txt))+
    facet_wrap(~taxa, nrow=3) +
    geom_boxplot()

ggplot(param.summary, aes(x=year,y=mean.biomasstoseed, fill=trt))+
    facet_wrap(~taxa, nrow=3) +
    geom_boxplot()
