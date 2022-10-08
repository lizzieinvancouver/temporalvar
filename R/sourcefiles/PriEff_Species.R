#Define the species characteristics for this run
#CONSIDER:
#   may want to vary initial density, N0
#   Megan reparameterized tau_g to be a Poisson(start day + chill-dependent delay 

nsp <- 2  #this is a 2-species model

#create arrays for within and between year dynamics
N <- matrix(rep(0), nyrs, nsp)   # number of seeds prior to winter
N0 <- rep(100,nsp)             # initial density of seeds
B0 <- matrix(rep(0),nyrs,nsp)    # initial biomass in year y
Bfin <- matrix(rep(0),nyrs,nsp)  # end of season biomass in year y
Bout <- list()                   # holds within season dynamics for each year
ext <- 0.00001                   # extinction threshold for density

#converting from within-year to between-year dynamics
s <-  rep(0.8,nsp)      # seedbank survival overwinter
b <-  rep(1,nsp)        # biomass of seedling per seed
phi <- rep(0.5,nsp)    # conversion of end-of-season plant biomass to seeds

#within-year competition parameters
#  note: vary Rstar using c, but don't allow R* to be neg (c>m*u)
a <-  rep(0.8,nsp)                    # slope of species uptake rate with increasing R
u <-  rep(5,nsp)                    # inverse of the max uptake rate
theta <- rep(1,nsp)                  # nonlinearity in resource response
m <-  rep(0.005,nsp)                 # mortality
c <- runif(nsp,max(m*u),3*max(m*u))  # conversion of resource to biomass
Rstar <- (m/(a*(c-m*u)))^(1/theta)


#germination timing tau_g (describes days of delay as a function of weeks of chilling)
#   tau_g is days of germ delay; avg min delay = 2; max delay ~15, depends on xi 
tau_start <- 2              # average start day, poisson
xi_tau <- runif(nsp,0,3)    # delay sensitivity to chill (up to 3 days per week of chill)
tau_delay <- t(xi_tau%*%t(xi)) # delay depends on chill
tau_g <-rpois(nsp*nyrs,tau_start+tau_delay)
dim(tau_g) <- dim(tau_delay)
#  tau_spr is a (days x spp) matrix that spreads the germination over the season
#          it sums to 1 for each species in each year (first pass has a 1 on one day, 0 elsewhere)
tau_spr <- list()
for (yr in c(1:nyrs)){
  ts1 <- data.frame(x=seq(0,days,by=dt),y = rep(0,days+1))
  ts2 <- data.frame(x=seq(0,days,by=dt),y = rep(0,days+1))
  ts1[ts1$x==tau_g[yr,1],2] = 1
  ts2[ts2$x==tau_g[yr,2],2] = 1
  tau_spr[[yr]] <- list(as.data.frame(ts1),as.data.frame(ts2))
}

#germination fraction g (describes germination rate as a function of chilling)
#   g increases at rate gamma_g from gmin to gmax, where gmin, gmax, and gamma_g are species-specific
#   CONSIDER - these values may need more thought - not a lot of variability bt yrs
gmin.zero <- 0.4                            #prob that g0 is a true zero
gmin <- ifelse(rbeta(nsp, 1,8)<gmin.zero,0,1)*rbeta(nsp, 1,8) #min germination with g0.zero true zeros
gmax <- rbeta(nsp, 8,1)
gmax<-ifelse(gmax-gmin<0,rbeta(nsp, 8,1),gmax)#max germination ranges from g0 to 1

gamma_g <- c(runif(nsp,.05,.5))  #species-specific rate of decline from max to min germ

#note that g is the proportion of seeds germinating in this year (mulit by N0 to get total num seeds)
g <- t(gmin + t(1 - exp(xi %*% t(-gamma_g)))*(gmax-gmin)) 
