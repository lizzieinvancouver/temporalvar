# get species and community characteristics for run
# updated from getSimilarSpp.R to make alphas vary (tauI already varies)
# I also changed gmax, lambda and s so they worked for nsp>2 (I checked that they still work n=2 also)

# tauI taken from Chesson et al. 2004 (mu_tauP is 0.5, so tauI closer to tauP is better)
if (nsp==2) tauI <- c(0.42, 0.4) else tauI <-runif(nsp,0.1, 0.9)  # time of max germ for sp i

gmax <-  rep(0.8, nsp) # max germination fraction

#Define plant parameters
#interaction coefficients (alpha matrix): the ith ROW is the competition experienced by sp i
alphadata <- runif(nsp*nsp, 0.1, 0.9)
alpha <- matrix(data = alphadata, nrow=nsp, ncol=nsp) # matrix fills in order (1,1),(2,1),(2,1),(2,2)
lambda <- rep(5, nsp) # seeds produced per germinant
s <- rep(0.8, nsp) # seedbank survival for non-germinating seeds
