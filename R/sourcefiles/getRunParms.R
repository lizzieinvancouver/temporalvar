#Runtime Parameters

#between year
nyrs <- 100  # number of yrs to run if nonsta=0
if (nonsta > 0) {
  nyrs <- nyrs + nonsta
}
y <- c(1:nyrs)
ext <- 1/10000  #Extinction Threshold:  1 seed/ha (assuming that initial density is 10 seeds per meter)

#within year
ndays <- 5  # number of days in a growing season
dt <- 0.005 # within yr timestep
tsteps <- ndays/dt
