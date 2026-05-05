# To do...
# tau_{g,i} which varies by species (taug)
# chilling hours threshold (varies by species with some species always going the same no matter the chill hours) -- Lizzie could start with simple linear f(x) -- see below
# calculate \hat(tau_g_i) -- for each species given chilling (ghat, see below)
# g(t) = d/(1 + (t/tau_g)^b_g)    g(t) is germination fraction on day t
#  where d =     f1(xi) = g0(i) + (gmax(i)-g0(i))*(1-exp(-gamma_g(i)*xi(y)))  d is max germ in year y, depends on xi(y)
#        tau_g = f2(xi) = tau_max - (tau_max/xi_tau(i))  tau_max is max delay same for all species
#                                                        xi_tau(i) chill weeks that result in min tau_g for species i
#germination
#timing as a function of xi (=chill)
tau_max <- 30               #max delay in germination, same for all species; think more about value (y-int of tau_g v chill)
xi_tau <- runif(nsp,0,16)   #chill weeks for min delay - species specific (x-intercept of tau_g v chill)
#tau_g is the delya i germination in year(y) defined by the line, where y=tau_g, yint = tau_max, slope=tau_max/xi_tau, and x=chill=xi
tau_g <- tau_max - (tau_max/matrix(rep(xi_tau,each=nyrs),nrow=nyrs, ncol=nsp))*matrix(rep(xi,nsp),nrow=length(xi),ncol=nsp)

#germ frac as a function of xi (=chill)
g0 <- runif(nsp,0,1)      #species-specific min germination, depending on chill
                          ##do we need to have a probability of a true 0 (won't get any real zeros with a runif)
gmax <- c(runif(1,g0[1],1),
          runif(1,g0[2],1)) #species-specific max germination, depending on chill, must be larger than g0
                           #clunky way to code b/c doesn't allow nsp to vary
gamma_g <- c(runif(nsp,0.2,1.25))  #species-specific rate of decline from max to min germ
#d is the amount of germination in year y 
d <- matrix(rep(g0,each=nyrs),nrow=nyrs,ncol=nsp) + matrix(rep(gmax-g0,each=nyrs),nrow=nyrs,ncol=nsp)*(1-exp(matrix(rep(-gamma_g,each=nyrs),nrow=nyrs, ncol=nsp)* matrix(rep(xi,nsp),nrow=length(xi),ncol=nsp)))  #this needs same matrix fussing (or could figure out the use of %*%)

if(FALSE){ ##we did this to check that the maxtriz multiplication  for d (above) was doing what we think
goo<-(1-exp(matrix(rep(-gamma_g,each=nyrs),nrow=nyrs, ncol=nsp)* matrix(rep(xi,nsp),nrow=length(xi),ncol=nsp)))
goo2<-matrix(rep(gmax-g0,each=nyrs),nrow=nyrs,ncol=nsp)

goo3<-matrix(rep(g0,each=nyrs),nrow=nyrs,ncol=nsp)
##NEXT - model of germination g(t) determines time of entry into the model at what level
head(d,1)

g0[1]+goo2[1,1]*goo[1,1]
}

## Lizzie started editing this on 19 July 2021, but definitely not done ...

#getSpecies characteristics
#write out Species Parms and Envt Parms

#survivorship
s <-  rep(0.8,nsp)      # seedbank survival overwinter

#biomass conversion
b <-  rep(1,nsp)          # biomass of seedling
phi <- rep(0.05,nsp)     # conversion of end-of-season plant biomass to seeds

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

#concatenate all the species-specific parms to save (not the time-varying parms) 
sppvars <- as.data.frame(cbind(b, s, phi, a, u, c, m, theta, Rstar, g0, gmax, gamma_g, tau_xi, tau_max))

#write out species parameters at the beginning of each run;
#  include headers if this is the first run
if (j==1) {
  col.names.SpeciesParms <- c("jobID","arrayID","runID",
                            paste0(rep(names(sppvars),each=nsp),
                                   rep(c(1:nsp),times=length(sppvars))))
  col.names.EnvtParms <- c("jobID","arrayID","runID","yr",names(envtvars),
                           paste0(rep("tau_g",nsp),c(1:nsp)),
                           paste0(rep("d",nsp),c(1:nsp)))
} else {
  col.names.SpeciesParms <- FALSE  
  col.names.EnvtParms <- FALSE
}
#this needs to be cleaned up so that it is writing out the new parameters
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

