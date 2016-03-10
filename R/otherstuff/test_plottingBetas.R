a = 12
b = 51


x = seq(0,1,.01)
xp = seq(0,1,1/10)
mu <- matrix(NA,1,length(xp))
sig <- matrix(NA,1,length(xp))

plot(x,dbeta(x,a,b),xlim=c(0,1),type="n")

for (i in c(0:10)){
  ahat <- a+(b-a)*i/10
  bhat <- b-(b-a)*i/10
  lines(x,dbeta(x,ahat,bhat),type="l")
  mu[i] <- ahat/(ahat+bhat)
  sig[i] <- ahat*bhat/((ahat+bhat)^2*(ahat+bhat+1))
}

plot(xp,mu,type="p")
plot(xp,sig,type="p")

