# get species and community characteristics for run

# tauI taken from Chesson et al. 2004 (mu_tauP is 0.5, so tauI closer to tauP is better)
if (nsp==2) tauI <- c(0.42, 0.4) else tauI <-runif(nsp,0.1, 0.9)  # time of max germ for sp i

gmax <-  c(.8,.8) #rep(0.5,nsp)   # max germination fraction

#Define plant parameters
#interaction coefficients (alpha matrix): the ith ROW is the competition experienced by sp i
alpha <- matrix(data = c(0.06,0.03,0.03,0.06), nrow=nsp, ncol=nsp)  #matrix fills in order (1,1),(2,1),(2,1),(2,2)
lambda <- c(5,5)    #seeds produced per germinant
s <-  c(0.8,0.8)      # seedbank survival for non-germinating seeds
