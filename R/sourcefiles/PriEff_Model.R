#Runs model for nyrs
#  year ends by counting the seed set in the fall = N[y-1]
#  year y starts after seed count in the fall N[y-1] 
#  --> overwinter survival in year y, s, & chilling in year y
#  --> germination, g[y,] --> conversion to seedling biomass, b
#  --> biomass entering growing season,B0[y,] --> growth
#  --> biomass at end of growing season, Bfin[y], converted to seeds by phi
#  --> count seeds at N[y,]

for (y in seq(1,nyrs)){
  B0[y,] <- ifelse(y==1, N0*s*g[y,]*b, N[y-1,]*s*g[y,]*b)
  R <- R0[y]
  B1 <- 0
  B2 <- 0
  State <- c(R=R,B1=B1,B2=B2)
  Time <- seq(0,days,by=1)
  tau.index <- c(which.min(tau_g[y,]),which.max(tau_g[y,]))
  germevents <- data.frame(var = Bnames[tau.index],
                           time = tau_g[y,tau.index],
                           value = B0[y,tau.index], 
                           method = c("add", "add"))
  
  Bout[[y]] <- as.data.frame(ode(func = RstarComp,
                           y = State,
                           parms = Pars,
                           times = Time, 
                           events=list(data=germevents),
                           rootfun=rootfun))
  source(here("R","sourcefiles","PriEff_PlotSeason.R"))

  #iterate N: seeds that survived but did not germinate + new seeds
  Bfin[y,] <- c(tail(Bout[[y]]$B1,1),tail(Bout[[y]]$B2,1))
  N[y,] <- Bfin[y,]*phi + ifelse(y==1,N0*s*(1-g[y,]), N[y-1,]*s*(1-g[y,]))  
  N[y,] <- N[y,]*(N[y,] > ext)            # call very low densities true zero
  if(isFALSE(sum(N)>0)) break             # if all species have gone extinct, stop
}
