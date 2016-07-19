#plot tauIini vs Rstar with points coded for coexistence


## Helpful reminders on indexing within lists from Lizzie!
# how to call something (given that you have at least 2 runs) ....
modelruns[[2]][1][[1]]$gmax
# how to call internal list by name
modelruns[[2]]["sppvars"]
# when slicing out bits of vector within a list you have to get beyond the $ 
# this is often something I forget so here are two examples
modelruns[[2]]["sppvars"][[1]]$gmax
modelruns[[2]]["tauP"][[1]][1:10]


gtauIPini <- matrix(0,nruns)

for (j in seq(1, nruns)){
numcoexist[j,1] <- sum(Bfin[(nonsta[1]-2),]>0)
if (nonsta[2]>0) numcoexist[j,2] <- sum(Bfin[(nonsta[1]+nonsta[2]-2),]>0)
if (nonsta[3]>0) numcoexist[j,3] <- sum(Bfin[(nyrs-2),]>0) 

gtauIPini <- matrix(0,nruns,nsp)
gtauIPns <- matrix(0,nruns,nsp)
gtauIPfin <- matrix(0,nruns,nsp)
gRstar <- matrix(0,nruns,nsp)

for (i in seq(1,nruns)){
  gtauIPini[i,] <- modelruns[[i]][[1]]$tauIPini
  gtauIPns[i,] <- modelruns[[i]][[1]]$tauIPns
  gtauIPfin[i,] <- modelruns[[i]][[1]]$tauIPfin  
  gRstar[i,] <- modelruns[[i]][[1]]$Rstar
}
plot(gRstar[,1]-gRstar[,2], gtauIPini[,1]-gtauIPini[,2],col=numcoexist+1,main="initial")
plot(gRstar[,1]-gRstar[,2], gtauIPns[,1]-gtauIPns[,2],col=numcoexist+1, main="nonstationary")
plot(gRstar[,1]-gRstar[,2], gtauIPfin[,1]-gtauIPfin[,2],col=numcoexist+1,main = "final")

##
## plotting what tracking looks like in comparison to tauP
##
nonsta = c(200,200,0)   #number of [1] initial stationary,[2]nonstationary,[3]final nonstationary years
tracking = 1   #tracking in these runs?
varRstar = 1   #flag for variation in Rstar; if 1, then c is drawn randomly, and R* varies
nsp = 2        #Number of species to start in these runs?
nyrs <- sum(nonsta)  # number of yrs to run if nonsta=0 or for initial period if nonsta>0

p <- 10  #first parameter for beta distribution of tau
q <- 10  #second parameter for beta distribution of tau
tauP <- rbeta(nonsta[1], p, q) # change once not doing stationary+nonstationary run

#nonstationary tauP includes nonstationary period of length nonsta[1] and stationary "final" period of length nonsta[2]
if (sum(nonsta[2:3]) > 0) {
  pfin <- 5
  qfin <- 15   #we should think about what is the appropriate value for qfin??
  pns <- seq(p, pfin, length.out=nonsta[2])
  qns <- seq(q, qfin, length.out=nonsta[2])
  tauPns <- rbeta(nonsta[2], pns, qns) #get tau during period of nonstationarity
  tauPfin <- rbeta(nonsta[3],pfin, qfin) #get tau during period after nonstationarity when q = qfin
  tauP <- c(tauP, tauPns,tauPfin)
  #plot(tauP~c(1:nyrs))
}

tauI <-runif(nsp,0.2,0.2)  # time of max germ for sp i
alpha <- rep(0,nsp)
#add tracking with alpha to create tauIhat
if (tracking > 0) {
  alpha <- runif(nsp,0, 0.99) # what we actually do
  alpha <- c(0.1,0.9) # added! make the alphas super different by over-writing random draw (right above)
  tauI <- matrix(rep(alpha),nyrs,nsp, byrow = TRUE)*tauP+matrix((1-alpha)*tauI, nyrs, nsp, byrow = TRUE) # this seems correct
  # tauI <- matrix(rep(alpha),nyrs,nsp, byrow = TRUE)*tauP+(1-alpha)*tauI # this seems wrong matrix work
}

plot(tauP~c(1:nyrs), type="l", ylim=c(-0.2,1.2))
lines(tauI[,1]~c(1:nyrs), col="blue", ylim=c(-0.2,1.2))
lines(tauI[,2]~c(1:nyrs), col="firebrick", ylim=c(-0.2,1.2))
