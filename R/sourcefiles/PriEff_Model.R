#Runs model for nyrs
#  year starts with seed set in the fall

for (y in yrs){
  B0[y,] <- N[y,]*s*g[y,]*b        #fall popln * overwinter surv * spring germ rate * convert to seedling biomass
  Bini <- c(0,0)
  State<- c(R=R0[y],B=Bini)
  Time <- seq(0,days,by=dt) 
  Pars <- c(c=c,a=a,u=u,m=m,theta=theta,eps=eps)
  germevents <- data.frame(var = c("B1","B2"),
                          time = tau_g,
                          value = B0[y,], 
                          method = c("add", "add"))
  out <- as.data.frame(ode(y = State,
                           func = RstarComp,
                           parms = Pars,
                           times = Time, 
                           events=list(data=germevents),
                           rootfun=rootfun))
 plot(out) 

  N[y+1] <- N[y]+ Bfin[y]*phi  
  #CONSIDER: impose a threshold on min density
     #N[y,] <- N[y,]*(N[y,]>ext)  #if density does not exceed ext, set to zero
  if (!(sum(N[y,]>0))) break      #if all species have gone extinct, stop
    }

  Bfin[y,] <-  apply(Bout[[y]][3:(2+nsp)],2,FUN=max)  #final biomass
}

