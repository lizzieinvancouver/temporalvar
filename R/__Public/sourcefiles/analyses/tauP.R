## Mainly taken from getEnvt.R ##
## Adjusted to return the tauP stat and tauP nonsta ##

nhere <- seq(0,1, by=0.0001)
#Resource parameters
# stationary
p <- 10  #first parameter for beta distribution of tau
q <- 10  #second parameter for beta distribution of tau
tauP <- dbeta(nhere, p, q) 
# nonstat
pfin <- 5
qfin <- 15 
tauPfin <- dbeta(nhere, pfin, qfin)

tauP.df <- data.frame(tauP=tauP, when="stat")
tauPfin.df <- data.frame(tauP=tauPfin, when="nonstatfin")

tauP.plot <- rbind(tauP.df, tauPfin.df)

if(FALSE){
plot(density(tauP)) # plots the results
lines(density(tauPfin), col="red")
ggplot(tauP.plot, aes(x=tauP, fill=when)) + geom_density(alpha=0.25)
}

# plot(seq(0,1, by=0.0011), dbeta(seq(0,1, by=0.0011), 5,15))
