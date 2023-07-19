#Define Environment for this run
#NOTE 1: the timescale of within year dynamics is determined by the time that R>R*
#       and depends on the relative values of R0.mu, eps, a,m,u.  
#       Megan adjusted species parms so that the season length is reasonable 
#       for the range of R0. Adjust eps, R0 to shorten/extend the season

#Growing season
days <- 120            # length of growing season in days
dt <- 1                # timestep for within year solver (4x per day)

#Resources in year y growing season
R0.mu <- log(100)                       # mean of start-of-season resource
R0.sigma <- log(2)                       # sd of start-of-season resources
R0 <- rlnorm(nyrs, R0.mu, R0.sigma)
eps <- 0.001                              # resource decay rate aside from uptake by 2 spp

#Weeks of chilling prior to growing season in year y

c_warm<- 0.5

if (runif(1,0,1)<c_warm) {


xi.mu <- log(8)                       # mean of chilling distribution ## 
xi.sigma <- .4                      # sd of chilling distribution
xi <- rlnorm(nyrs, xi.mu, xi.sigma)
}else {
  
  xi.mu <- log(4)                       # mean of chilling distribution ## 
  xi.sigma <- .4                      # sd of chilling distribution
  xi <- rlnorm(nyrs, xi.mu, xi.sigma)  
}# chilling accumulated before each season
# ADD CODE for copula to allow covariance between R0 and xi