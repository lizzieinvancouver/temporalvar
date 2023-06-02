#Define the species characteristics for this run

nsp <- 2  #this is a 2-species model


# R* COMPETITION parameters -----------------------------------------------

#create arrays for within and between year dynamics
N <- matrix(rep(0), nyrs, nsp)   # number of seeds prior to winter
N0 <- rep(100,nsp)           # initial density of seeds $THis makes it a one species model
B0 <- matrix(rep(0),nyrs,nsp)    # initial biomass in year y
Bfin <- matrix(rep(0),nyrs,nsp)  # end of season biomass in year y
Bout <- list()                   # holds within season dynamics for each year
ext <- 0.0001                   # extinction threshold for density

#converting from within-year to between-year dynamics
s <-  rep(0.8,nsp)      # seedbank survival overwinter 
b <-  rep(1,nsp)        # biomass of seedling per seed
phi <- rep(0.5,nsp)    # conversion of end-of-season plant biomass to seeds

#within-year competition parameters
#  note: vary Rstar using c, but don't allow R* to be neg (c>m*u)
a <-  rep(0.8,nsp)                    # slope of species uptake rate with increasing R
u <-  rep(5,nsp)                    # inverse of the max uptake rate
theta <- rep(1,nsp)                  # nonlinearity in resource response
m <-  rep(0.005,nsp)                 # mortality
c <- runif(nsp,max(m*u),30*max(m*u))  # conversion of resource to biomass
Rstar <- (m/(a*(c-m*u)))^(1/theta)


# TIMING ------------------------------------------------------------------
#germination timing tau_g (describes days of delay as a function of weeks of chilling)
#   maximum mean delay from chilling is 30 days
#   mean delay decreases exponentially at rate xi_tau with increasing days of chilling, xi
#   cut off max delay at 90 days
# 
xi_tau <- runif(nsp,0,1)    # delay sensitivity to chill - exponential rate of decreased delay from 30
tau_delay <- 30*exp(-t(xi_tau%*%t(xi)))  
tau_g50 <-rpois(nsp*nyrs,tau_delay)  #Day of 50% germination
tau_g50 <- pmin(tau_g50,90)          #max delay is 90 days
dim(tau_g50) <- dim(tau_delay)

#test plots
#histogram of day of 50% germination
par(mfrow=c(3,2))
ax <- seq(0,max(tau_g50),1)
h1 <- hist(tau_g50[,1],breaks=ax,plot= FALSE)
h2 <- hist(tau_g50[,2],breaks=ax,plot= FALSE)
if(FALSE){
plot(h1,col=alpha(1,0.3))
plot(h2,col=alpha(2,0.5),add=TRUE)
#Day of 50% germination versus chilling (xi)

plot(xi,tau_g50[,1], xlim=c(0,ceiling(max(xi))), ylim=c(0,max(tau_g50)),
     xaxs="i",yaxs="i", col=1,
     ylab="Day of 50% germination",xlab="chilling (xi)",
     main="Day of 50% Germ ~ Chilling, 2 spp, nyrs")
points(tau_g50[,2]~xi, col=2)
points(tau_delay[,1]~xi, col=1, pch=20)
points(tau_delay[,2]~xi, col=2, pch=20)
}

# FRACTION ----------------------------------------------------------------
#germination fraction g (describes germination rate as a function of chilling)
#   20% of runs have both spp chill-insensitive w.r.t germination fraction; these have gmax = 0.8
#   in remaining 80% of runs, 20% of species will be insensitive with const gmax that varies between 0.5 and 1
#   remaining species are chilling sensitive w.r.t. germination fraction
#     gmin - germ w 0 weeks chilling; 0 for all chilling-sensitive species 
#     max possible germination for all chilling-sensistive species is 100%, but annual gmax will vary with chill
#     xi_0 is chilling weeks required to exceed 0 germination ~rpois(2) 
#     xi_rng is the number of chilling weeks (range) from xi_0 to xi_100 
#     xi_100 = xi_0 + xi_rngis the number of chilling until  species reaches 100% germination
#     gmax is the total germination fraction in year yr with chilling xi[yr]

g_notxi <- 0.2  #proportion of runs where both species are %g- insensitive to chilling
                 # and proportion of species that are %g-insens in remaining runs
xi_100 <- 0

if (runif(1,0,1)<g_notxi) {
  #25% of runs have both species are chilling insensitive at have gmax=0.8
  gmin = rep(0,nsp)
  gmax = matrix(rep(0.8,times=nsp*nyrs),ncol=nsp)
  #print(gmin)
} else {
  #remaining 75% of runs have a 25% insensitive species and 75% sensitive species
  #25% of species are insensitive and get a fixed min germ between 0.5-1, rest are sensitive and get min germ of 0
  gmin <- as.numeric(runif(nsp,0,1)<g_notxi) * runif(nsp,0.5,1)  #if insensitive to xi, get gmin between 0.5,1
  xi_0 <- rnorm(nsp,2,1)    #weeks of chilling to exceed 0 germination
  while (isFALSE(all(xi_0>0))) xi_0 <- rnorm(nsp,2,1)  #truncates normal to >0
  xi_rng <- rnorm(nsp,6,2) 
  xi_100 <- xi_0 + xi_rng  #chilling weeks required for max germ
  
  #max germination for chilling sensitive species
  gmax <- xi %*% t(1/xi_rng) - matrix(rep(xi_0/xi_rng,nyrs),nrow=nyrs,byrow=TRUE)
  gmax <- pmin(gmax,1)    #max germination fraction is 1
  gmax <- pmax(gmax,0)    #max germination fraction has minimum of 0
  #now, replace gmax with gmin for species that were flagged as chilling insensitive
  #  and replace gmin for insensitive species with 0 so that gmin will be 0 for all species in the dose-response curve
  if(gmin[1]>0) {
    gmax[,1] <- gmin[1]
    gmin[1] <-0
  }
  if(gmin[2]>0) {
    gmax[,2] <- gmin[2]
    gmin[2] <-0
  }
}

#test plots to show relationship between germ and chill
#  with the annual chill values on the x
###close plot for speeding computation maybe 
if(FALSE){
for (i in c(1,nsp)){
  plot(x=seq(0,max(xi),length=10), y= seq(0,1,length=10), 
       type="n", xlim=c(0,40),ylim=c(0,1), xaxs="i",yaxs="i",
       xlab="chilling",ylab=paste("%g for species", i),
       main = paste("%g-Chilling Sensitivity for species", i))
  points(x=xi,y=rep(0,length(xi)),col="black",pch=3,cex=0.8)
  if (var(gmax[,i])==0) {  #if gmax is constant
    abline(a=gmax[1,i],b=0,col="blue") 
  } else { 
    points(x=xi_0[i],y=0,col=i,pch=20,cex=2)
    points(x=xi_100[i],y=1,col=i,pch=20, cex=2)
    abline(b=1/xi_rng[i],a=-xi_0[i]/xi_rng[i],col="blue", lty=1)
  }
}
}

# DISTRIBUTE Germination over Days of Season ------------------------------
#  tau_spr is a (days x spp) matrix that spreads  germination over season 
#     one version where gmax[yr] is all on day tau_g50
#     one version where germination is Hill fxn w 50% germ on tau_g50
g_daily <- list()
g_cumulative <- list()
nH = 5  #Hill constant for germination by Hill Eqn
maxdailygerm <- 0 #use for plotting

for (yr in c(1:nyrs)){
  g1_byday <- data.frame(x=seq(0,days,by=dt),y = rep(0,days+1))
  g2_byday <- data.frame(x=seq(0,days,by=dt),y = rep(0,days+1))
  gc1 = rep(0,days+1)
  gc2 = rep(0,days+1) #cumulative germination on day d
  
  ##Version with all germination on one day (tau_g50)
  # ts1[ts1$x==tau_g50[yr,1],2] = gmax[1]
  # ts2[ts2$x==tau_g50[yr,2],2] = gmax[2]
  #Germination follows a Hill function (Dose-Response) with min = 0
  for (d in seq(1,days+1,by=1)){
    gc1[d]  <- ifelse(d==1,0,gmin[1] + (gmax[yr,1]-gmin[1]) * (d^nH / (tau_g50[yr,1]^nH + d^nH)))
    gc2[d]  <- ifelse(d==1,0,gmin[2] + (gmax[yr,2]-gmin[2]) * (d^nH / (tau_g50[yr,2]^nH + d^nH)))
    g1_byday$y[d] <- ifelse(d==1, 0, gc1[d] - gc1[d-1])
    g2_byday$y[d] <- ifelse(d==1, 0, gc2[d] - gc2[d-1])
  }
  maxdailygerm<- max(maxdailygerm,max(g1_byday$y),max(g2_byday$y))
  g_daily[[yr]] <- list(data.frame(g1_byday),data.frame(g2_byday))
  g_cumulative[[yr]] <- data.frame(gc1,gc2)
  #Test Plot for cumulative germination
 if(FALSE){
   if (yr==1) {
    plot(seq(0,days,1),rep(0,days+1),type="n",
         ylim=c(0,1.1),,xlim= c(0,80),
         xlab="days",ylab="cumulative germination by d",
         main="Cumulative Germination by sp, yrs",
         xaxs="i",yaxs="i")
  }
  lines(seq(0,days,1),gc1,type="l",col=1)
  lines(seq(0,days,1), gc2,type="l",col=2)
}
  }

#Test Plot for daily germination
if(FALSE){
for (yr in seq(1,nyrs,5)){
  if (yr==1) {
  plot(seq(0,days,1),rep(0,days+1),type="n",
       ylim=c(0,maxdailygerm+.01),xlim= c(0,80),
       xlab="days",ylab="germination on d",main="Daily Germination by sp, yrs",
       xaxs="i",yaxs="i")
  }
  lines(g_daily[[yr]][[1]]$x,g_daily[[yr]][[1]]$y,type="l",col="black")
  lines(g_daily[[yr]][[2]]$x,g_daily[[yr]][[2]]$y,type="l",col="blue")
  #print(paste("summed daily g sp1 = ", sum(g_daily[[yr]][[1]]$y),"gmax = ",gmax[yr,1]))
  #print(paste("summed daily g sp2 = ", sum(g_daily[[yr]][[2]]$y),"gmax = ",gmax[yr,2]))
}
}
