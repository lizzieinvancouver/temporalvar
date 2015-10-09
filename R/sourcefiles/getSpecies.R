#getSpecies characteristics

#survivorship
s <-  rep(0.8,nsp)      # seedbank survival overwinter

#germination
b <-  rep(1,nsp)          # biomass of seedling
gmax <-  rep(0.5,nsp)     # max germination fraction
h <-  rep(100,nsp)             # max rate of germination decrease following pulse
phi <- rep(0.05,nsp)     # conversion of end-of-season plant biomass to seeds
#germination: tau I and alpha below
if (nsp==2) tauIhat <- c(0.35, 0.4) else tauIhat <-runif(nsp,0.1, 0.9)  # time of max germ for sp i
tauIhat <- matrix(rep(tauIhat),nyrs,nsp, byrow = TRUE)
alpha <- runif(nsp,0.3, 0.99)
# alpha <- rep(0, nsp) # turns off alpha
alpha <- matrix(rep(alpha),nyrs,nsp, byrow = TRUE)
tauI <- alpha*tauP+(1-alpha)*tauIhat 

g <- gmax*exp(-h*(matrix(rep(tauP,nsp),nrow=length(tauP),ncol=nsp)-tauI)^2)  #germination fraction in year y

#competition
a <-  rep(20,nsp)        # slope of species uptake rate with increasing R
u <-  rep(1,nsp)          # inverse of the max uptake rate
c <-  rep(12,nsp)        # conversion of resource to biomass
# c <- runif(nsp, 2, 20) # this leads to much lower coexistence when I checked it on 6 Oct 2015
theta <- rep(1,nsp)      #nonlinearity in resource response
m <-  rep(0.05,nsp)     # mortality
Rstar <- (m/(a*(c-m*u)))^(1/theta)

#concatenate all the parms to save
crossyrsvars <- as.data.frame(cbind(b, s, a, u, c, m, gmax, h, phi, theta, tauI=colMeans(tauI[,]), Rstar))
