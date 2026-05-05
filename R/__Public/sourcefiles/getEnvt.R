## To do ...
# Generate chilling hours (need to decide on distribution) 
# Generate pulse size (we alreaday have this, see R_0) 
# Eventually allow these to covary?
# Generate evapotranspration (could alternatively covary with chilling hours)

# Have done (July 2021, by Lizzie) 
# Removed \tau_p (start times always the same now)
# Removed all the stuff about megadroughts and non-stationarity (I think Dan can do more of an ANOVA set of studies, instead of adding in non-stationarity)
# ... as of meeting on 20 July 2021 may need to revisit edits here


#Resource parameters
#time varying
Rmu <- log(2)  #mean of resource distribution
Rsigma <- 0.2  #sd of resource distribution
R0 <- rlnorm(nyrs, Rmu, Rsigma) # intial R in a season
ximu <- log(2)  #mean of chilling distribution
xisigma <- 0.2  #sd of chilling distribution
xi <- rlnorm(nyrs, ximu, xisigma) # chilling accumuulated before each season
##add a copula for covariance between R0 and xi

#constant (for now) - can use this to extend the timing of the growing season
eps <- 1              # evaporative stress 

envtvars <- as.data.frame(cbind(R0, xi, eps))
