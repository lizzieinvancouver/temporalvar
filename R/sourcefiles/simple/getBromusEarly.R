# get species and community characteristics for run
# taken from Godoy & Levine 2014 (Ecology)

# tauI taken from Chesson et al. 2004
if (nsp==2) tauI <- c(0.35, 0.4) else tauI <-runif(nsp,0.1, 0.9)  # time of max germ for sp i

gmax <-  c(0.32,0.12) #rep(0.5,nsp)   # max germination fraction

#Define plant parameters
alpha <- matrix(data = c(0.122,0.0132,0.0117,0.0018), nrow=nsp, ncol=nsp)  #interaction coefficients; the ith row is the competition experienced by sp i
lambda <- c(1993,23)    #seeds produced per germinant
s <-  c(0.02,0.31)      # seedbank survival for non-germinating seeds
