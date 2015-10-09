#Resource parameters
#time varying
mu <- log(2)  #mean of resource distribution
sigma <- 0.2  #sd of resource distribution
R0 <- rlnorm(nyrs, mu, sigma) # intial R in a season
p <- 2  #first parameter for beta distribution of tau
q <- 2  #second parameter for beta distribution of tau
tauP <- rbeta(nyrs, p, q) # change once not doing stationary+nonstationary run

#constant (for now)
eps <- 1              # evaporative stress 

#nonstationary tauP
if (nonsta > 0) {
  qns <- seq(2, 20, length.out=nonsta)
  tauPns <- rbeta(nonsta, p, qns) # yes, it takes a vector! Yay!
  tauP <- tauPns
}
