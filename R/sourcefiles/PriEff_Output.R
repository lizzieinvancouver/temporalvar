###Dan tries
###run dynamics
mydat<-data.frame(sp1_Rstar=Rstar[1],sp2_Rstar=Rstar[2],
             
                   sp1_xi_tau=xi_tau[1],sp2_xi_tau=xi_tau[2],
                 run=nruns[j])
                   outputy<-rbind(mydat,outputy)
                 
###across year dynamics                     
                mydat2<-data.frame(sp1_gmax=gmax[,1],sp2_gmax=gmax[,2],
                                      sp1_g50=tau_g50[,1], sp2_g50=tau_g50[,2],
                                      sp1_B0=B0[,1],sp2_B0=B0[,2],
                                      sp1_Bfin=Bfin[,1],sp2_Bfin=Bfin[,2],year=1:nyrs,run=nruns[j])
                  
                 outputy2<-rbind(mydat2,outputy2)
                 
###within year
                 ###not sure how to set this only for certain yeats
Bout.trip <- mapply(cbind, Bout, "year"=1:nyrs,SIMPLIFY=F) #  assign that number
Bout.df<-do.call(rbind.data.frame, Bout.trip) ### make the list a data frame                 
Bout.df$run<-nruns[j]
outputy3<-rbind(Bout.df,outputy3)                 