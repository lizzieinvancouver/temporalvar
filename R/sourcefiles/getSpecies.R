#getSpecies characteristics
#write out Species Parms and Envt Parms

#survivorship
s <-  rep(0.8,nsp)      # seedbank survival overwinter

#germination
b <-  rep(1,nsp)          # biomass of seedling
gmax <-  rep(0.5,nsp)     # max germination fraction
h <-  rep(100,nsp)             # max rate of germination decrease following pulse
phi <- rep(0.05,nsp)     # conversion of end-of-season plant biomass to seeds

#megaDrought - tradeoff phi and surv with correlation rho=-0.5
if (megaD==1) {
  cmat <- matrix(c(1,rho,rho,1), nrow=2, ncol=2) 
  sphi <- draw.d.variate.uniform(no.row=1,d=2,cov.mat=cmat)
  s <- sphi[,1]*(0.95 - 0.05) + 0.05
  phi <- sphi[,2]*(0.1-0.01) + 0.01
}
  

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

if (!is.na(varRstar[2])) {     #if varRstar is vector then it gives c for each species
  c <- varRstar
} else if (varRstar[1]==0) {      #if varRstar==0 then is gives the same (randomly generated) R* for all species in the run
  c <- rep(runif(1, 2, 20),nsp) 
} else if (varRstar[1] == -1){     #if varRstar is -1, then go with old default value for all species
    c <-  rep(12,nsp)
} else {                       #otherwise, randomly select c for each species
    c <- runif(nsp,8,20)
}

a <-  rep(20,nsp)        # slope of species uptake rate with increasing R
u <-  rep(1,nsp)          # inverse of the max uptake rate
theta <- rep(1,nsp)      #nonlinearity in resource response
m <-  rep(0.05,nsp)     # mortality
Rstar <- (m/(a*(c-m*u)))^(1/theta)

#concatenate all the parms to save
sppvars <- as.data.frame(cbind(b, s, phi, a, u, c, m, theta, Rstar, gmax, h, alpha, tauI,tauIPini,tauIPns,tauIPfin))

#write out species parameters at the beginning of each run;
#  include headers if this is the first run
if (j==1) {
  col.names.SpeciesParms <- c("jobID","arrayID","runID",
                            paste0(rep(names(sppvars),each=nsp),
                                   rep(c(1:nsp),times=length(sppvars))))
  col.names.EnvtParms <- c("jobID","arrayID","runID","yr",names(envtvars),
                           paste0(rep("tauIhat",nsp),c(1:nsp)),
                           paste0(rep("g",nsp),c(1:nsp)))
} else {
  col.names.SpeciesParms <- FALSE  
  col.names.EnvtParms <- FALSE
}

write.table(matrix(data=c(as.numeric(jobID[1]),as.numeric(jobID[2]),j,as.matrix(sppvars)),nrow=1),
            file = paste0(OtherOut_loc,"/SpeciesParms",suffix),
            col.names = col.names.SpeciesParms,row.names = FALSE, 
            append=TRUE,sep="\t", quote=FALSE)
write.table(matrix(data=c(rep(as.numeric(jobID[1]),nyrs),rep(as.numeric(jobID[2]),nyrs),rep(j,nyrs),yrs,
                          envtvars$R0,envtvars$tauP,envtvars$eps,tauIhat,g),
                   nrow=nyrs,ncol=(nsp*2+7)),
            file = paste0(OtherOut_loc,"/EnvtParms",suffix),
            col.names =col.names.EnvtParms,row.names = FALSE,
            append=TRUE,sep="\t", quote=FALSE) 

