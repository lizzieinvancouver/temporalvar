#Environment calculations

#mean, variance, skew, etc of resource pulse (for calulating var(Envt) and mu(Envt))
mu_tauP <- p/(p+q)   #single value
sig2_tauP <- p*q/((p+q)^2*(p+q+1))  #single value
sig2_tauP2 <- sig2_tauP^2 + 2*mu_tauP*sig2_tauP  #var(tauP^2)
skew_tauP <- 2*(q-p)*sqrt(p+q+1)/((p+q+2)*sqrt(p*q))
mu_tauP3 <- skew_tauP*(sig2_tauP)^(3/2) + mu_tauP^3 + 3*mu_tauP*sig2_tauP

#Environment experienced by each species
E <- log(g)  #row = nyrs, col= nsp
#expected value of Envt = <ln g(i,t)>: nrows = 1, ncols = nsp
mean_E <- log(gmax) - h*(sig2_tauP + mu_tauP^2 - 2*tauI*mu_tauP + tauI^2)
#expected value of Envt based on "observed" (actual germination values) rather than derived from distribution)
mean_Ehat <- apply(log(g),2,mean)
#variance of environment: nrows = 1, ncols = nsp
var_E <- h^2*(sig2_tauP2 + 4*tauI^2*sig2_tauP - 4*tauI*(mu_tauP3-mu_tauP*(sig2_tauP + mu_tauP^2)))
var_Ehat <- apply(E,2,var)

#########################################
#Competition Calculations

#C_noi[y,i] is competition experienced by species i when rare (i.e., no intraspecific comp)
C_noi <- matrix(NA, nrow=nyrs,ncol=nsp)

for (y in c(1:nyrs)){
  for (i in c(1:nsp)){
    ind <- seq(1:nsp)[-c(i)]
    a_noi <- sum(alpha[i,ind]*g[y,ind]*N[y,ind])
    C_noi[y,i] <- -log(lambda[i]/(1+ a_noi) -s[i])  #should this have a negative sign? YES
  }
}

