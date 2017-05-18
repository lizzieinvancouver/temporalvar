#getSpecies characteristics

#survivorship
s <-  rep(0.8,nsp)      # seedbank survival overwinter

#germination
b <-  rep(1,nsp)          # biomass of seedling
gmax <-  rep(0.5,nsp)     # max germination fraction
h <-  rep(100,nsp)             # max rate of germination decrease following pulse
phi <- rep(0.05,nsp)     # conversion of end-of-season plant biomass to seeds

#germination: tau I and alpha below; tauI is the time of max germ for sp i
if (length(vartauI)>1) {  #if vartauI is a vector, then it is giving particular values for each species
  tauI = vartauI
} else if (vartauI == 0) {  #if vartauI is 0, then give all species the same randomly selected tauI
  tauI <-rep(runif(1,0.1,0.9),nsp)
} else {                     #if vartauI is 1, then give random values for tauI for each species
  tauI <-runif(nsp,0.1,0.9)  
}

#tracking
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
} else {
  tauIPfin <- rep(NA,nsp)
}

g <- gmax*exp(-h*(matrix(rep(tauP,nsp),nrow=length(tauP),ncol=nsp)-tauIhat)^2)  #germination fraction in year y

#competition
# c is conversion of resource to biomass; we vary this to vary R*

if (length(varRstar)>1) {     #if varRstar is vector then it gives c for each species
  c <- varRstar
} else if (varRstar==0) {      #if varRstar==0 then is gives the same (randomly generated) R* for all species in the run
  c <- rep(runif(1, 2, 20),nsp) 
} else if (varRstar == -1){     #if varRstar is -1, then go with old default value for all species
    c <-  rep(12,nsp)
} else                        #otherwise, randomly select c for each species
    c <- runif(nsp,2,20)
}

a <-  rep(20,nsp)        # slope of species uptake rate with increasing R
u <-  rep(1,nsp)          # inverse of the max uptake rate
theta <- rep(1,nsp)      #nonlinearity in resource response
m <-  rep(0.05,nsp)     # mortality
Rstar <- (m/(a*(c-m*u)))^(1/theta)

#concatenate all the parms to save
sppvars <- as.data.frame(cbind(b, s, phi, a, u, c, m, theta, Rstar, gmax, h, alpha, tauI,tauIPini,tauIPns,tauIPfin))
