#Runtime Parameters

#between year
nyrs <- sum(nonsta)  # number of yrs to run if nonsta=0 or for initial period if nonsta>0

y <- c(1:nyrs)
ext <- 1/10000  #Extinction Threshold:  1 seed/ha (assuming that initial density is 10 seeds per meter)

#within year
ndays <- 5  # number of days in a growing season
dt <- 0.005 # within yr timestep
tsteps <- ndays/dt
