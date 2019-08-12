## Started 12 August 2019 ##
## By Lizzie ##

## Working on idealized plotting ## 

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)


## Plotting the ideal tauP distributions ...
x <- seq(0, 1, length = 10000)
plot.start.tauP <- dbeta(x, 10, 10)
plot.end.tauP <- dbeta(x, 5, 15)
plot(plot.start.tauP, type="l", ylab="", xlab="", yaxt="n", ylim=c(0, 4.25))
lines(plot.end.tauP, ylab="", xlab="", yaxt="n", col="blue")
mean(rbeta(10000, 10, 10))
mean(rbeta(10000, 5, 15))

## Next! 
taui.spa <- 0.65
taui.spb <- 0.99
alpha <- 0.5
tauP <- 0.5
tauihat.spb <- alpha*tauP+(1-alpha)*taui.spb

# arrows()
