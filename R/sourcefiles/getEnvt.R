#Resource parameters
#time varying
mu <- log(2)  #mean of resource distribution
sigma <- 0.2  #sd of resource distribution
R0 <- rlnorm(nyrs, mu, sigma) # intial R in a season
p <- 10  #first parameter for beta distribution of tau
q <- 10  #second parameter for beta distribution of tau
tauP <- rbeta(nyrs, p, q) # change once not doing stationary+nonstationary run

#constant (for now)
eps <- 1              # evaporative stress 

#nonstationary tauP includes nonstationary period of length nonsta[1] and stationary "final" period of length nonsta[2]
if (nonsta[1] > 0) {
  pfin <- 5
  qfin <- 15   #we should think about what is the appropriate value for qfin??
  pns <- seq(p, pfin, length.out=nonsta[1])
  qns <- seq(q, qfin, length.out=nonsta[1])
  tauPns <- rbeta(nonsta[1], pns, qns) #get tau during period of nonstationarity
  tauPfin <- rbeta(nonsta[2],pfin, qfin) #get tau during period after nonstationarity when q = qfin
  tauP <- c(tauP, tauPns,tauPfin)
  plot(tauP~c(1:(nyrs+sum(nonsta))))
  nyrs <- nyrs + sum(nonsta)  #when a nonstationary period is added to a stationary 
}
