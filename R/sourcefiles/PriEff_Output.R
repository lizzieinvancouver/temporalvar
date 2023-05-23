###Dan tries
###run dynamics
#mean max, var, gmax, xi_100, mean tau g50
mydat<-data.frame(sp1_Rstar=Rstar[1],sp2_Rstar=Rstar[2],
             sp1_meangmax=mean(gmax[,1]),sp1_SDgmax=sd(gmax[,1]),
             sp2_meangmax=mean(gmax[,2]),sp2_SDgmax=sd(gmax[,2]),
             sp1_xi100=xi_100[1],sp2_xi100=xi_100[2],
               sp1_xi_tau=xi_tau[1],sp2_xi_tau=xi_tau[2],
             sp1_mean_tau_g50=mean(tau_g50[,1]),sp1_SD_tau_g50=sd(tau_g50[,1]),
             sp2_mean_tau_g50=mean(tau_g50[,2]),sp2_SD_tau_g50=sd(tau_g50[,2]),
             sp1_ex=N[nyrs,1],sp2_ex=N[nyrs,2],
                 
             run=nruns[j])
                   outputy<-rbind(mydat,outputy)
       ###add all paramenters (env too)
                   ##write my dat _run number
###across year dynamics   
####add: amount of chilling, R_intial                   
                   
              #  mydat2<-data.frame(xi=xi,sp1_gmax=gmax[,1],sp2_gmax=gmax[,2],
  #                                    sp1_g50=tau_g50[,1], sp2_g50=tau_g50[,2],
   #                                   sp1_B0=B0[,1],sp2_B0=B0[,2],
    #                                  sp1_Bfin=Bfin[,1],sp2_Bfin=Bfin[,2],sp1_exist=N[,1],sp2_exist=N[,2],year=1:nyrs,run=nruns[j])
                  
               #  outputy2<-mydat2
                 
###within year
                 ###not sure how to set this only for certain yeats
#Bout.trip <- mapply(cbind, Bout, "year"=1:nyrs,SIMPLIFY=F) #  assign that number
#Bout.df<-do.call(rbind.data.frame, Bout.trip) ### make the list a data frame                 
#Bout.df$run<-nruns[j]
#outputy3<-rbind(Bout.df,outputy3)                 
###write it out and close
###
#outputy4<-filter(outputy3, year %in% seq(5,nyrs,by=10))

#path_out<- "R/output/within_year_dynams/"
#fileName <-paste(path_out, "iter",nruns[j],sep = '')
#write.csv(outputy4,fileName)

#rm(outputy4)
#head(outputy3)
write.csv(outputy,"R/output/prieff_params.csv",row.names=FALSE)
#write.csv(outputy2,"R/output/prieff_acrossyr_params.csv",row.names=FALSE)
print(paste("done",run=nruns[j]))
