### Started 28 Jan 2015 ###
### By Lizzie & Megan ###

## Simpler model for within year dynamics (L-V within year, with variable between year germination##

# safety feature(s)
setwd(getwd()) # Lizzie: setwd("~/Documents/git/temporalvar/R")
options(stringsAsFactors=FALSE)

# packages
library(ggplot2)
library(deSolve)

set.seed(2)

#Number of species
nsp = 2  #when nsp=2, tauI is assigned known values from chesson 2004

#Number of years to run
nyrs <- 100

#Define arrays
#interannual dynamics set-up (R0 is in getEnvt.R)
N0 <- rep(10,nsp)          # initial number of seeds (per meter square?)
N <- matrix(rep(0), nyrs, nsp) # number of seeds by yr and spp
N[1,] <- N0  #initialize
#rcrt <- matrix(rep(0),nyrs,nsp) # recruitment in year y
#rcrt0 <- matrix(rep(0),nyrs,nsp) # recruitment WO competition in year y

#Define environment parameters & germination
p <- 2  #first parameter for beta distribution of tau
q <- 2  #second parameter for beta distribution of tau
tauP <- rbeta(nyrs, p, q) # change once not doing stationary+nonstationary run
gmax <-  c(.32,.12)#rep(0.5,nsp)   # max germination fraction
h <-  rep(100,nsp)      # max rate of germination decrease following pulseif (nsp==2) tauI <- c(0.35, 0.4) else tauI <-runif(nsp,0.1, 0.9)  # time of max germ for sp i
g <- gmax*exp(-h*(matrix(rep(tauP,nsp),nrow=length(tauP),ncol=nsp)-matrix(rep(tauI,nyrs),ncol=nsp,nrow=nyrs,2))^2)  #germination fraction in year y


#Define plant parameters
alpha <- matrix(data = c(.122,.0017,.0132,.0018),nrow=nsp, ncol=nsp)  #interaction coefficients; the ith column is the competition experienced by sp i
lambda <- c(1993,23)    #seeds produced per germinant
s <-  c(.02,.31)      # seedbank survival for non-germinating seeds


for (y in c(1:(nyrs-1))){
  C <- rep(NA,nsp)
  for (i in c(1:nsp)){
    C[i] <- sum(alpha[,i]*g*N[y,])
  }
  N[y+1,] <- N[y,]*s*(1-g[y,]) + lambda*g[y,]*N[y,]/(1+C)
  N[y+1,] <- N[y+1,]*(N[y+1,]>ext)  #if density does not exceed ext, set to zero
}

#most basic plot
plot(c(1:nyrs),N[,1],type="l")
lines(c(1:nyrs),N[,2],col="red")



# source("sourcefiles/plotNyears.R")  #plots dynamics of seedbank abundance over years
# #source("sourcefiles/plotBinSeason.R")  #plot within season dynamics of biomass & R for a subset of years
# source("sourcefiles/plotBwCnoC.R")  #plot within season biomass resource dynamics w and wo competition
# source("sourcefiles/plotBinSeason_Lizzie.R")
# 
# 
# dev.new(width=7, height=6)
# 
# # (2) using ggplot, which really is good for this sort of thing
# tau.df <- data.frame(coexisted=Bfin[max(y),]>0, tauI=tauI)
# tau.df$coexisted[tau.df$coexisted==TRUE] <- "coexisted"
# tau.df$coexisted[tau.df$coexisted==FALSE] <- "doomed"
# tauP.df <- data.frame(coexisted=rep("tauP"), tauP=tauP)
# 
# if (nonsta>0){
#   ggplot(tau.df, aes(tauI, fill = coexisted)) + geom_histogram(alpha=0.5) +
#     geom_density(data=tauP.df[1:(nyrs-nonsta),], aes(tauP),  alpha = 0.2) +
#     geom_density(data=tauP.df[(nonsta+1):nyrs,], aes(tauP),  alpha = 0.4) +
#     labs(title=paste(sum(Bfin[max(y),]>0), "out of", nsp, "coexisted", sep=" "))
# }else {
#   ggplot(tau.df, aes(tauI, fill = coexisted)) + geom_histogram(alpha=0.5) +
#     geom_density(data=tauP.df[1:nyrs,], aes(tauP),  alpha = 0.2) +
#     labs(title=paste(sum(Bfin[max(y),]>0), "out of", nsp, "coexisted", sep=" "))
# }
# 
# 
# print(c("The number of coexisting species are",sum(Bfin[max(y),]>0),"out of",nsp))
