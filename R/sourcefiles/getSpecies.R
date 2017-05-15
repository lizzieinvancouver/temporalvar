#getSpecies characteristics

#survivorship
s <-  rep(0.8,nsp)      # seedbank survival overwinter

#germination
b <-  rep(1,nsp)          # biomass of seedling
gmax <-  rep(0.5,nsp)     # max germination fraction
h <-  rep(100,nsp)             # max rate of germination decrease following pulse
phi <- rep(0.05,nsp)     # conversion of end-of-season plant biomass to seeds
#germination: tau I and alpha below
tauI <-runif(nsp,0.1,0.9)  # time of max germ for sp i
alpha <- rep(0,nsp)
#add tracking with alpha to create tauIhat
if (tracking > 0) {
  alpha <- runif(nsp,0.3, 0.99)
  tauIhat <- matrix(rep(alpha),nyrs,nsp, byrow = TRUE)*tauP+matrix((1-alpha)*tauI, nyrs, nsp, byrow = TRUE)
} else {
  tauIhat <- matrix(rep(tauI),nyrs,nsp, byrow = TRUE)
} 

#effective tauI (=tauIhat) for initial, nonstationary, and final periods 
tauIPini <- colMeans(abs(tauP[1:nonsta[1]] - tauIhat[1:nonsta[1],]))
if (nonsta[2]>0) {
  tauIPns <- colMeans(abs(tauP[(nonsta[1]+1):(nonsta[1]+nonsta[2])] - tauIhat[(nonsta[1]+1):(nonsta[1]+nonsta[2]),]))
} else {  
  tauIPns <-  rep(NA,nsp)
}
if (nonsta[3]>0){
  tauIPfin <- colMeans(abs(tauP[(nonsta[1]+nonsta[2]+1):nyrs] - tauIhat[(nonsta[1]+nonsta[2]+1):nyrs,]))
}else{
  tauIPfin <- rep(NA,nsp)
}

g <- gmax*exp(-h*(matrix(rep(tauP,nsp),nrow=length(tauP),ncol=nsp)-tauIhat)^2)  #germination fraction in year y

#competition
if (varRstar>0) {
  c <- runif(nsp, 2, 20) # this leads to much lower coexistence when I checked it on 6 Oct 2015
} else {
    c <-  rep(12,nsp)        # conversion of resource to biomass
}
a <-  rep(20,nsp)        # slope of species uptake rate with increasing R
u <-  rep(1,nsp)          # inverse of the max uptake rate
theta <- rep(1,nsp)      #nonlinearity in resource response
m <-  rep(0.05,nsp)     # mortality
Rstar <- (m/(a*(c-m*u)))^(1/theta)

#concatenate all the parms to save
sppvars <- as.data.frame(cbind(b, s, phi, a, u, c, m, theta, Rstar, gmax, h, alpha, tauI,tauIPini,tauIPns,tauIPfin))
