##September 22, 2014
##The purpose of this file is to create a "map" of the ratio of BiomassIn/BiomassOut for each species
##as a function of BiomassIn_Sp1, BiomassIn_Sp2,R0
##This should give us a better grasp of the within season dynamics, so that we can ignore them and look at the 
##between season dynamics where the action is

library(deSolve)
nonsta=0
nsp = 2  
set.seed(1)

source("sourcefiles/getRunParms.R") #define runtime parameters
source("sourcefiles/getGraphParms.R")  #define graphics parameters
source("sourcefiles/getEnvt.R")  #get constant and time-varying envt parms
source("sourcefiles/getSpecies.R")  #get species characteristics and Rstar

Bin<- matrix(rep(seq(0,1,.1),nsp),ncol=2)
R0 <- seq(1,4,.1)
B1out <- array(NA,c(length(Bin)/nsp,length(Bin)/nsp,length(R0)))
B2out <- array(NA,c(length(Bin)/nsp,length(Bin)/nsp,length(R0)))

## Within-season dynamics set-up- ODE
source("sourcefiles/ResCompN.R") # define within-season ode solver AND defines Pars
source("sourcefiles/NoCompN.R")  # define within-season ode solver for no competition
Time <- seq(0,1,by=dt)

for (r in seq(1,length(R0))) {
  for (b1 in seq(1,length(Bin)/nsp)) {
    for (b2 in seq(1,length(Bin)/nsp)) {
      R<-R0[r]
      B<-c(Bin[b1,1],Bin[b2,2])
      State <- c(R=R,B=B)
      Brun <- as.data.frame(ode(func = ResCompN, y = State, parms = Pars, times = Time))
      B1out[b1,b2,r]<-max(Brun$B1)
      B2out[b1,b2,r]<-max(Brun$B2)      
    }
  }
}

#Confusing Contour plots of resulting biomass of sp1 and sp2 as a function of initial biomass of sp1, sp2, and R0
png("../figures/WithinYear_BiomassOutvBiomassInvR.png",width=6, height=8,units="in",res=400)
par(mfrow=c(2,2))
for (r in c(1,11,21,31)) {
  contour(Bin[,1],Bin[,2],B1out[,,r],xlab="BiomassIn_Sp1",ylab="BiomassIn_Sp2",main = paste("BiomassOut for R0 = ",R0[r]),col="blue",lwd=2)
  contour(Bin[,1],Bin[,2],B2out[,,r],add=TRUE,col="red",lwd=2)
  legend("bottomright",c("Bfin_Sp1","Bfin_Sp2"),bty="o", col=c("blue","red"),lwd=2,cex=0.7)
}
dev.off()




# #Plot dynamics of a single year
# par(mar = c(5,4,4,4+.3))
# plot(Brun$B1~Brun$time,type="l",col="red",xlim = c(0,4),ylab="Biomass",xlab="Time")
# points(Brun$B2~Brun$time,type="l", col="Blue")
# par(new=TRUE)
# plot(Brun$R~Brun$time,type="l",axes=FALSE,bty="n",xlab="",ylab="")
# axis(side=4,at=pretty(range(Brun$R)))
# mtext("Resource",side=4,line=3)

