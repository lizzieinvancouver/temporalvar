#Runs model for nyrs
#  year ends by counting the seed set in the fall = N[y-1]
#  year y starts after seed count in the fall N[y-1] 
#  --> overwinter survival in year y, s, & chilling in year y
#  --> germination, g[y,] --> conversion to seedling biomass, b
#  --> biomass entering growing season,B0[y,] --> growth
#  --> biomass at end of growing season, Bfin[y], converted to seeds by phi
#  --> count seeds at N[y,]

for (y in seq(1,nyrs)){
  if(y>1) N0<-N[y-1,]
  B0[y,] <- N0*s*g[y,]*b
  R <- R0[y]
  B1 <- 0
  B2 <- 0
  State <- c(R=R,B1=B1,B2=B2)
  Time <- seq(0,days,by=dt)
  Pars <- c(c=c,a=a,u=u,m=m,theta=theta,eps=eps,s=s, g=g[y,],b=b,N0=N0)
  Germinate1 <-approxfun(tau_spr[[y]][[1]], method = "constant",rule=2)
  Germinate2 <-approxfun(tau_spr[[y]][[2]], method = "constant",rule=2)
  Bout[[y]] <- ode(func = RstarComp,
                           y = State,
                           parms = Pars,
                           times = Time, 
                           Germinate1 = Germinate1,
                           Germinate2 = Germinate2)
  #step when R went below R* for each species
  ind.Rstar<- c(min(which(Bout[[y]][,"R"]<Rstar[1]),120),min(which(Bout[[y]][,"R"]<Rstar[2]),120))
  print(ind.Rstar)
  source(here("R","sourcefiles","PriEff_PlotSeason.R"))

  #iterate N: seeds that survived but did not germinate + new seeds
  #get Bfin the first step where R<R*min
  Bfin[y,] <- c(Bout[[y]][,"B1"][ind.Rstar[1]],Bout[[y]][,"B2"][ind.Rstar[2]])
  N[y,] <- Bfin[y,]*phi + ifelse(y==1,N0*s*(1-g[y,]), N[y-1,]*s*(1-g[y,]))  
  N[y,] <- N[y,]*(N[y,] > ext)            # call very low densities true zero
  if(isFALSE(sum(N)>0)) break             # if all species have gone extinct, stop
}
