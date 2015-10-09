#Runtime Parameters

#between year
nyrs <- 100  # number of yrs to run if nonsta=0
if (nonsta > 0) {
  nyrs <- 0 + nonsta # replacing nyrs with 0 on 6 Oct 2015 to do just a nonstationary run,
                      # otherwise could not get it to run
}
y <- c(1:nyrs)
ext <- 1/10000  #Extinction Threshold:  1 seed/ha (assuming that initial density is 10 seeds per meter)

#within year
ndays <- 5  # number of days in a growing season
dt <- 0.005 # within yr timestep
tsteps <- ndays/dt
