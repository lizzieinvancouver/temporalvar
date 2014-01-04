getSpecies <- function(nsp)

b <-  rep(1,nsp)          # biomass of seedling
s <-  rep(0.8,nsp)      # seedbank survival overwinter
a <-  rep(20,nsp)        # slope of species uptake rate with increasing R
u <-  rep(1,nsp)          # inverse of the max uptake rate
c <-  rep(12,nsp)        # conversion of resource to biomass
m <-  rep(0.05,nsp)     # mortality
gmax <-  rep(0.5,nsp)     # max germination fraction
h <-  rep(100,nsp)             # max rate of germination decrease following pulse
phi <- rep(0.05,nsp)     # conversion of end-of-season plant biomass to seeds
tauI <- runif(nsp,0.1, 0.9)    # time of max germ for sp i

return(b,s,a,u,c,m,gmax,h,phi,tauI)