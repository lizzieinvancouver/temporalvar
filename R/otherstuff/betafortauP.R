### Started 10 March 2016 ###
### By Lizzie ###

## Lizzie learns about beta distribution ##


## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## some setup
library(e1071) # skewmness parameter

betaskew <- function(alp, bet){
    return((2*(bet-alp)*(sqrt(alp+bet+1)))/((alp+bet+2)*(sqrt(alp*bet))))
}
    

## here we go
# vary both!
brange <- c(1:100)
alpharange <- c(0.5, 1, 2, 3, 5, 10, 15, 20, 30, 50, 100, 150) 

par(mfrow=c(3,4))

for (j in c(1:length(alpharange))){
    meanz <- c()
    skewz <- c()
    alphahere <- alpharange[j]
    for(i in c(1:length(brange))){
        betame <- rbeta(30000, alphahere, brange[i])
        meanz[i] <- mean(betame)
        skewz[i] <- skewness(betame)
    }
    plot(meanz~brange, main=paste("alpha is", alphahere), ylim=c(-2,1))
    points(skewz~brange, col="dodgerblue")
}    


# vary beta only
setalpha <- 5
brange <- c(1:200)

quartz()
meanz <- c()
skewz <- c()
betaskewz <- c()

for(i in c(1:length(brange))){
   betame <- rbeta(30000, setalpha, brange[i])
   meanz[i] <- mean(betame)
   skewz[i] <- skewness(betame)
   betaskewz[i] <- betaskew(setalpha, brange[i])
    }
plot(meanz~brange, main=paste("alpha is", alphahere), ylim=c(-2,1))
points(skewz~brange, col="dodgerblue")
points(betaskewz~brange, col="red")


# options close to 0.8125 without too much skew
hist(rbeta(3000,100,22), xlim=c(0,1))
hist(rbeta(3000,50,11), xlim=c(0,1))

# options close to 0.1875 without too much skew
hist(rbeta(3000,3,156), xlim=c(0,1))
hist(rbeta(3000,5,22), xlim=c(0,1))
hist(rbeta(3000,10,43)) # opposite is rbeta(3000,10,43)

# some means
mean(rbeta(3000,11,50))
mean(rbeta(3000,50,11))
mean(rbeta(3000,12,51))
mean(rbeta(3000,51,12))

# calculate how closer we're getting to goal of 0.125 variation on either side of mean
qbeta(0.99, 10,43)-qbeta(0.5, 10,43)
qbeta(0.01, 10,43)-qbeta(0.5, 10,43)

qbeta(0.99, 50, 11)-qbeta(0.5, 50,11)
qbeta(0.01, 50, 11)-qbeta(0.5, 50,11)

qbeta(0.99, 51, 12)-qbeta(0.5, 51,12)
qbeta(0.01, 51, 12)-qbeta(0.5, 51,12)

# overlay histograms
h1 <- rbeta(3000,12,51)
h2 <- rbeta(3000,51,12)

h3 <- rbeta(3000,11,50)
h4 <- rbeta(3000,50,11)

hist(h1, col=rgb(1,0,0,0.5), xlim=c(0,1))
hist(h2, col=rgb(0,0,1,0.5), add=TRUE)

quartz() # WTF?
hist(h3, col=rgb(1,0,0,0.5), xlim=c(0,1))
hist(h4, col=rgb(0,0,1,0.5), add=TRUE)

