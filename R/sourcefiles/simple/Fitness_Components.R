#Fitness components

#######################
#Mean fitness:  rbar_i_prime/di

meanFit <- rep(NA,nsp)
del_muI <- rep(NA,nsp)   #based on moments of beta distribution
del_muIhat <- rep(NA,nsp)  #based on samples from beta dist
del_vI <- rep(NA,nsp)
del_vIhat <- rep(NA,nsp)
mu <- mean_E - log(1-s)  #(r,c) = (1,nsp)
mu_hat <- mean_Ehat - log(1-s)  #(r,c) = (1,nsp)

for (i in c(1:nsp)){
  ind <- seq(1:nsp)[-c(i)]
  del_muI[i] <- mu[i] -mean(mu[ind])
  del_muIhat[i] <- mu_hat[i] -mean(mu_hat[ind])
  del_vI[i] <- (1/2)*s[i]*var_E[i] - mean((1/2)*s[ind]*var_E[ind])
  del_vIhat[i] <- (1/2)*s[i]*var_Ehat[i] - mean((1/2)*s[ind]*var_Ehat[ind])
}

meanFit <- del_muI + del_vI  #(r,c) = (1,nsp)
meanFit_hat <- del_muIhat + del_vIhat

##########################
#Storage Effect
#   (delI/d) is the difference between individual growth rate as invader 
#   and average community growth rate experienced by invader

#cov(E,C) = cov(log(g(t)),log(lambda/(a+sum(comp(t)))-surv)
storEff <- rep(NA,nsp)
covEC <- rep(NA,nsp)

for (i in c(1:nsp)){
  covEC[i] <- cov(E[,i],C_noi[,i])  #chi_j{-i}
}
#calculate storage effect for each species 
for (i in c(1:nsp)){
  ind <- seq(1:nsp)[-c(i)]
storEff[i] <- mean(s[ind]*covEC[ind]) - s[i]*covEC[i]
}

###########################
#Relative Nonlinearity

#relNL = 1/2 var(C_noi)*(di-avg(d))
relNL <- rep(NA,nsp)

for (i in c(1:nsp)){
  ind <- seq(1:nsp)[-c(i)]
  relNL[i] <- (1/2)*var(C_noi[,i])*(1-s[i]-mean(1-s[ind]))
}