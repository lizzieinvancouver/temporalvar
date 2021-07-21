# To do...
# tau_{g,i} which varies by species (taug)
# chilling hours threshold (varies by species with some species always going the same no matter the chill hours) -- Lizzie could start with simple linear f(x) -- see below
# calculate \hat(tau_g_i) -- for each species given chilling (ghat, see below)

# Stuff we did not change ... 
# g_max constant for all species (currently) -- leave as is
# s, h, phi -- also leave as is
# Rstar etc. stuff -- also leave as is

#getSpecies characteristics
#write out Species Parms and Envt Parms

#survivorship
s <-  rep(0.8,nsp)      # seedbank survival overwinter

#germination
b <-  rep(1,nsp)          # biomass of seedling
gmax <-  rep(0.5,nsp)     # max germination fraction
h <-  rep(100,nsp)             # max rate of germination decrease following pulse
phi <- rep(0.05,nsp)     # conversion of end-of-season plant biomass to seeds

#species-specific (ASK MEGAN! about the below onward ...)
tauc <-runif(nsp,0.01,0.2)
taug <- runif(nsp,0.1,0.5)
#okay ... I think ghat needs to be a matrix, similar to tauIhat...
# tauIhat <- matrix(rep(tauI),nyrs,nsp, byrow = TRUE)
ghat <- taug-tauc*xi

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

