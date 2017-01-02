#Resource parameters
#time varying
mu <- log(2)  #mean of resource distribution
sigma <- 0.2  #sd of resource distribution
R0 <- rlnorm(nyrs, mu, sigma) # intial R in a season
p <- 10  #first parameter for beta distribution of tau
q <- 10  #second parameter for beta distribution of tau
tauP <- rbeta(nonsta[1], p, q) # change once not doing stationary+nonstationary run

#constant (for now)
eps <- 1              # evaporative stress 

#nonstationary tauP includes nonstationary period of length nonsta[1] and stationary "final" period of length nonsta[2]
if (sum(nonsta[2:3]) > 0) {
  pfin <- 5
  qfin <- 15   #we should think about what is the appropriate value for qfin??
  pns <- seq(p, pfin, length.out=nonsta[2])
  qns <- seq(q, qfin, length.out=nonsta[2])
  tauPns <- rbeta(nonsta[2], pns, qns) #get tau during period of nonstationarity
  tauPfin <- rbeta(nonsta[3],pfin, qfin) #get tau during period after nonstationarity when q = qfin
  tauP <- c(tauP, tauPns,tauPfin)
  #plot(tauP~c(1:nyrs))
}

envtvars <- as.data.frame(cbind(R0,tauP, eps))