#quick plots of within season dynamics

par(mfrow=c(1,2))
plot(Bout[[y]]$R~Bout[[y]]$time, type="l",
     xlab="days",ylab=NA, main="Resource")
plot(Bout[[y]]$B1~Bout[[y]]$time, type="l", ylim=c(0,max(Bout[[y]]$B1,Bout[[y]]$B2)),
     xlab="days",ylab=NA,main="Sp1 & Sp2 Density")
lines(Bout[[y]]$B2~Bout[[y]]$time, type="l",col="blue")
