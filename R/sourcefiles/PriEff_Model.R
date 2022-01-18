#Runs model for nyrs
#  year starts with seed set in the fall

for (y in yrs){
  B0[y,] <- N[y,]*s*g[y,]*b        #fall popln * overwinter surv * spring germ rate * convert to seedling biomass
  R=R0[y]
  B=c(0,0)
  State<- c(R=R0[y],B=B0[y,])#c(R=R,B=B)#
  Time <- seq(0,days,by=dt) 
  Pars <- c(c=c,a=a,u=u,m=m,theta=theta,eps=eps)
  tau.index <- c(which(tau_g==min(tau_g)),which(tau_g==max(tau_g)))
  germdates <- data.frame(var = paste0(rep("B",2),tau.index),
                          time = tau_g[tau.index],
                          value = B0[y,tau.index], 
                          method = c("add", "add"))
  out <- ode(y = State,func = RstarComp,parms = Pars,times = Time)
#  ,events=list(data=germdates),
#             rootfun=rootfun)
 plot(out) 

  N[y+1] <- N[y]+ Bfin[y]*phi  
  #CONSIDER: impose a threshold on min density
     #N[y,] <- N[y,]*(N[y,]>ext)  #if density does not exceed ext, set to zero
  if (!(sum(N[y,]>0))) break      #if all species have gone extinct, stop
    }

  Bfin[y,] <-  apply(Bout[[y]][3:(2+nsp)],2,FUN=max)  #final biomass
}

