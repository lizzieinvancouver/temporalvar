#Define the species characteristics for this run
#CONSIDER:
#   may want to vary initial density, N0
#   Megan reparameterized tau_g to be a Poisson(start day + chill-dependent delay 

nsp <- 2  #this is a 2-species model

#create arrays for within and between year dynamics
N <- matrix(rep(0), nyrs, nsp)    # number of seeds prior to winter
N[1,] <- rep(10,nsp)             # initial density of seeds
B0 <- matrix(rep(0),nyrs,nsp)     # initial biomass in year y
B <- list()                       # B contains the within season dynamics   

#converting from within-year to between-year dynamics
s <-  rep(0.8,nsp)      # seedbank survival overwinter
b <-  rep(1,nsp)        # biomass of seedling per seed
phi <- rep(0.05,nsp)    # conversion of end-of-season plant biomass to seeds

#within-year competition parameters
#  note: vary Rstar using c, but don't allow R* to be neg (c>m*u)
a <-  rep(0.2,nsp)                    # slope of species uptake rate with increasing R
u <-  rep(10,nsp)                    # inverse of the max uptake rate
theta <- rep(1,nsp)                  # nonlinearity in resource response
m <-  rep(0.005,nsp)                 # mortality
c <- runif(nsp,max(m*u),3*max(m*u))  # conversion of resource to biomass
Rstar <- (m/(a*(c-m*u)))^(1/theta)

#germination timing tau_g (describes days of delay as a function of weeks of chilling)
#   tau_g is days of germ delay; avg min delay = 2; max delay ~15, depends on xi 
tau_start <- 2              # average start day, poisson
xi_tau <- runif(nsp,0,4)    # delay sensitivity to chill (up to 4 days per week of chill)
tau_delay <- xi_tau%*%t(xi) # delay depends on chill
tau_g <-rpois(nsp,tau_start+tau_delay)

#germination fraction g (describes germination rate as a function of chilling)
#   g increases at rate gamma_g from gmin to gmax, where gmin, gmax, and gamma_g are species-specific
gmin.zero <- 0.1                                   #prob that g0 is a true zero
gmin <- ifelse(runif(nsp,0,1)<gmin.zero,0,1)*runif(nsp,0,1) #min germination with g0.zero true zeros
gmax <- runif(nsp,gmin,rep(1,nsp))                 #max germination ranges from g0 to 1
gamma_g <- c(runif(nsp,0.2,1.25))  #species-specific rate of decline from max to min germ
g <- t(gmin + t(1 - exp(xi %*% t(-gamma_g)))*(gmax-gmin))
