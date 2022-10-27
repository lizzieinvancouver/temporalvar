#Define the species characteristics for this run

nsp <- 2  #this is a 2-species model

#create arrays for within and between year dynamics
N <- matrix(rep(0), nyrs, nsp)   # number of seeds prior to winter
N0 <- rep(100,nsp)             # initial density of seeds
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


#germination timing tau_g (describes days of delay as a function of weeks of chilling)
#   tau_g is days of germ delay; avg min delay = 2; max delay ~15, depends on xi 
tau_start <- 2              # average start day, poisson
xi_tau <- runif(nsp,0,2)    # delay sensitivity to chill (up to 3 days per week of chill)
tau_delay <- t(xi_tau%*%t(xi)) # delay depends on chill
tau_g50 <-rpois(nsp*nyrs,tau_start+tau_delay)  #Day of 50% germination
dim(tau_g50) <- dim(tau_delay)

#test plot to show relationship between day of 50% germ and chilling (xi)
plot(xi,tau_g50[,1], xlim=c(0,ceiling(max(xi))), ylim=c(0,max(tau_g50)),
     xaxs="i",yaxs="i",
     ylab="Day of 50% germination",xlab="chilling (xi)",
     main="Day of 50% Germ ~ Chilling, 2 spp, nyrs")
points(tau_g50[,2]~xi, col="blue")

#germination fraction g (describes germination rate as a function of chilling)
#   g increases at rate gamma_g from gmin to gmax, where gmin, gmax, and gamma_g are species-specific
#2022/10/24 - perhaps we should think of 
#   gmin as a species char, defined as germ w 0 chilling
#   max possible germination for a species is 100%, but annual gmax will vary with chill
#   g_gamma as the sensitivity to chilling calculated from slope when
#   xi_100 is the chilling hours where species reaches 100% germination

# gmin.zero is the percent of species that will not germinate if zero chill
#    [QUESTION: do we care about this?  do we need a separate category here?]
# gmin is the minimum germination rate of species; 
# beta(1,8) has mean of 11% germination and 90%-ile of 25% germination
# gmax is the total germination fraction in year yr with chilling xi[yr]
gmin.zero <- 0
gmin <- rbeta(nsp, 1,10) * ifelse(runif(nsp, 0,1)<gmin.zero,0,1)
xi_100 <- rgamma(nsp,shape=3,scale=5)  #chilling reqd for 100% germ
gamma_g <- (1-gmin)/xi_100
gmax <- xi %*% t(gamma_g) + matrix(rep(gmin,nyrs),nrow=nyrs,byrow=TRUE)
gmax <- pmin(gmax,1)  

#test plots to show relationship between germ and chill
plot(x=seq(0,max(xi),length=10), y= seq(1,10), type="n",
     xlim=c(0,max(xi)),ylim=c(0,1), xaxs="i",yaxs="i",
     xlab="chiling",ylab="% germination with xi chilling",
     main = "Annual %Germ ~ Chilling, by species + rug of xi" )
points(x=xi,y=rep(0,length(xi)),col="cyan",pch=3)

for (i in seq(1,nyrs)){
  # gmin <- ifelse(runif(nsp, 0,1) < gmin.zero,0,1)*rbeta(nsp, 1,10)
  # xi_100 <- rgamma(nsp,shape=3,scale=5)
  # gamma_g <- (1-gmin)/xi_100
  abline(a=gmin[1],b=gamma_g[1])
  abline(a=gmin[2],b=gamma_g[2],col="blue")
}

#  tau_spr is a (days x spp) matrix that spreads  germination over season 
#     one version where gmax[yr] is all on day tau_g50
#     one version where germination is Hill fxn w 50% germ on tau_g50
g_daily <- list()
g_cumulative <- list()
nH = 5  #Hill constant for germination by Hill Eqn
maxdailygerm <-0 #use for plotting

for (yr in c(1:nyrs)){
  g1_byday <- data.frame(x=seq(0,days,by=dt),y = rep(0,days+1))
  g2_byday <- data.frame(x=seq(0,days,by=dt),y = rep(0,days+1))
  gc1 = rep(0,days+1)
  gc2 = rep(0,days+1) #cumulative germination on day d
  
  ##Version with all germination on one day (tau_g50)
  # ts1[ts1$x==tau_g50[yr,1],2] = gmax[1]
  # ts2[ts2$x==tau_g50[yr,2],2] = gmax[2]
  #Germination follows a Hill function (Dose-Response)
  for (d in seq(0,days,by=1)){
    gc1[d]  <- (gmin[1] + (gmax[yr,1]-gmin[1]) * (d^nH / (tau_g50[yr,1]^nH + d^nH)))*(d!=0)
    gc2[d]  <- (gmin[2] + (gmax[yr,2]-gmin[2]) * (d^nH / (tau_g50[yr,2]^nH + d^nH)))*(d!=0)
    g1_byday$y[d] <- ifelse(d==0, 0, gc1[d] - gc1[d-1])
    g2_byday$y[d] <- ifelse(d==0, 0, gc2[d] - gc2[d-1])
  }
  maxdailygerm<- max(maxdailygerm,max(g1_byday$y),max(g2_byday$y))
  g_daily[[yr]] <- list(data.frame(g1_byday),data.frame(g2_byday))
  g_cumulative[[yr]] <- data.frame(gc1,gc2)
  #Test Plot for cumulative germination
  if (yr==1) {
    plot(seq(0,days,1),rep(0,days+1),type="n",ylim=c(0,1.1),
         xlab="days",ylab="cumulative germination by d",
         xaxs="i",yaxs="i")
  }
  lines(seq(0,days,1),gc1,type="l",col="black")
  lines(seq(0,days,1), gc2,type="l",col="blue")
}

# #Test Plot for daily germination
# for (yr in seq(1,nyrs,5)){
#   if (yr==1) {
#   plot(seq(0,days,1),rep(0,days+1),type="n",ylim=c(0,maxdailygerm+.01),
#        xlab="days",ylab="germination",main="Daily Germination by sp, yrs",
#        xaxs="i",yaxs="i")
#   }
#   lines(g_daily[[yr]][[1]]$x,g_daily[[yr]][[1]]$y,type="l",col="black")
#   lines(g_daily[[yr]][[2]]$x,g_daily[[yr]][[2]]$y,type="l",col="blue")
#   #print(paste("summed daily g sp1 = ", sum(g_daily[[yr]][[1]]$y),"gmax = ",gmax[yr,1]))
#   #print(paste("summed daily g sp2 = ", sum(g_daily[[yr]][[2]]$y),"gmax = ",gmax[yr,2]))
# }
