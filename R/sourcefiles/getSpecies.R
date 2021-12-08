# To do...
# tau_{g,i} which varies by species (taug)
# chilling hours threshold (varies by species with some species always going the same no matter the chill hours) -- Lizzie could start with simple linear f(x) -- see below
# calculate \hat(tau_g_i) -- for each species given chilling (ghat, see below)
# g(t) = d/(1 + (t/tau_g)^b_g)    g(t) is germination fraction on day t
#  where d =     f1(xi) = g0(i) + (gmax(i)-g0(i))*(1-exp(-gamma_g(i)*xi(y)))  d is max germ in year y, depends on xi(y)
#        tau_g = f2(xi) = tau_max - (tau_max/xi_tau(i))  tau_max is max delay same for all species
#                                                        xi_tau(i) chill weeks that result in min tau_g for species i
#germination
#chilling and timing
tau_max <- 30               #max delay in germination, same for all species; think more about value
xi_tau <- runif(nsp,0,16)   #chill weeks for min delay - species specific 
tau_g <- tau_max - (tau_max/xi_tau)

#chilling and germ frac
g0 <- runif(nsp,0,1)      #species-specific min germination, depending on chill
                          ##do we need to have a probability of a true 0 (won't get any real zeros with a runif)
gmax <- c(runif(1,g0[1],1),
          runif(1,g0[2],1) #species-specific max germination, depending on chill, must be larger than g0
                           #clunky way to code b/c doesn't allow nsp to vary
gamma_g <- c(runif(nsp,0.2,1.25))  #species-specific rate of decline from max to min germ
d <- g0 + (gmax-g0)*(1-exp(-gamma_g * xi[y]))  #this prob won't work - think about vector mult

##NEXT - model of germination g(t) determines time of entry into the model at what level
          
h <-  rep(100,nsp)             # max rate of germination decrease following pulse


#species-specific (ASK MEGAN! about the below onward ...)
# ... as of meeting on 20 July 2021 may need to revisit edits here ... do we have a tau_c still etc.?

tauc <-runif(nsp,0.01,0.2)
taug <- runif(nsp,0.1,0.5)
#okay ... I think ghat needs to be a matrix, similar to tauIhat...
# tauIhat <- matrix(rep(tauI),nyrs,nsp, byrow = TRUE)
ghat <- taug-tauc*xi

g <- gmax*exp(-h*(matrix(rep(tauP,nsp),nrow=length(tauP),ncol=nsp)-tauIhat)^2)  #germination fraction in year y

# Stuff we did not change ... 
# g_max constant for all species (currently) -- leave as is
# s, h, phi -- also leave as is
# Rstar etc. stuff -- also leave as is

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

