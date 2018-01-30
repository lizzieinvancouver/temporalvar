#SIX output files are created each run. 
# 
#  RunParms_jobID.txt:  includes runparms is written in getRunParms.R
#  SpeciesParms_jobID.txt:  arrayID,runID,sppvars  is written in getSpecies.R
#  EnvtParms_jobID.txt:  arrayID,runID,yr,R0,tauP,eps,tauI_1,tauI_2,...,g_1,g_2...,
#                             is written in getSpecies.R
#TWO files are written at the end of each run year:
#  BfinN_jobID.txt:  arrayID,runID,yr,Bfin_1,Bfin_2...]
#  SummaryOut_jobID.txt:  
#BOUT file is written out every year, but there is a new file for each run
#  Bout_jobID.txt:  arrayID,runID,yr,t,B_1,B_2,...

#NOTE: outloc and suffix are defined in getRunParms.R

#write out headers only on first iteration


if (j == 1) {
  col.names.BfinN <- c("jobID","taskID","runID","yr","Ltstp",
                    paste0(rep("Bfin",nsp),c(1:nsp)),
                    paste0(rep("N",nsp),c(1:nsp)),
                    paste0(rep("rcrt",nsp),c(1:nsp)),
                    paste0(rep("coexist",nsp),c(1:nsp)))
  col.names.SummaryOut <- c("jobID","taskID","runID","yout","ncoexist",
                            paste0(rep("coexist",nsp),c(1:nsp)),
                            paste0(rep("alpha",nsp),c(1:nsp)),
                            paste0(rep("c",nsp),c(1:nsp)),
                            paste0(rep("Rstar",nsp),c(1:nsp)),
                            paste0(rep("g",nsp),c(1:nsp),rep("mean_pre",nsp)),
                            paste0(rep("tauIPini",nsp),c(1:nsp),rep("_pre",nsp)),
                            paste0(rep("tauIPns",nsp),c(1:nsp),rep("_pre",nsp)),
                            paste0(rep("tauIPfin",nsp),c(1:nsp),rep("_pre",nsp)))
} else {
  col.names.BfinN <- FALSE
  col.names.SummaryOut <- FALSE
}

#Calculate summary stats for output
L <- rep(0,yout)  #number of timesteps in Bout in year y
gmean_pre <- colMeans(gmax*exp(-h*(matrix(rep(tauP[1:yout],nsp),nrow=yout,ncol=nsp)-tauIhat[1:yout,])^2))  #germination fraction in year y

if (yout<=nonsta[1]) {
  tauIPini_pre <- colMeans(abs(tauP[1:yout] - tauIhat[1:yout,]))
  tauIPns_pre <- c(NA,NA)
  tauIPfin_pre <- c(NA,NA)
} else {
  if (yout <= nonsta[2]){
    tauIPini_pre <- tauIPini
    tauIPns_pre <- colMeans(abs(tauP[(nonsta[1]+1):yout] - tauIhat[(nonsta[1]+1):yout,]))
    tauIPfin_pre <- c(NA,NA)
  } else {
    tauIPini_pre <- tauIPini
    tauIPns_pre <- tauIPns
    tauIPfin_pre <- colMeans(abs(tauP[(nonsta[2]+1):yout] - tauIhat[(nonsta[2]+1):yout,]))
}
}
write.table(matrix(data=c(as.numeric(jobID[1]),as.numeric(jobID[2]),j,yout,sum(coexist[yout,]),coexist[yout,],
                          alpha,c,Rstar,
                          gmean_pre,tauIPini_pre,tauIPns_pre,tauIPfin_pre),nrow=1),
            file=paste0(outloc,"SummaryOut",suffix),
            col.names = col.names.SummaryOut, row.names = FALSE,
            append=TRUE,sep="\t", quote=FALSE)

#write Bout (by year) and Bfin
# note that y is the counter for years in PhenologyModel.r loop
# therefore, y is the number of years before loop breaks bc extinction
# m is the counter for year in this loop
# in Bfin, include number of timesteps in the ode run for each year (L[m])

#new Bout file gets written for each run
col.names.Bout <- FALSE

for (m in c(1:yout)) {
  L[m] <- dim(Bout[[m]])[1]
  ##only write headers in first year
  if (m==1) col.names.Bout <- c("jobID","taskID","runID","yr","time","R",paste0(rep("B",nsp),c(1:nsp)))
  write.table(matrix(data=c(rep(as.numeric(jobID[1]),L[m]),rep(as.numeric(jobID[2]),L[m]),
                            rep(j,L[m]),rep(m,L[m]),as.matrix(Bout[[m]])),
                     nrow=L[m],ncol=4+2+nsp),
              file = paste0(outloc,"Bout","_",jobID[1],"-",jobID[2],"-",j,".txt"),
              col.names = col.names.Bout,row.names = FALSE,
              append = TRUE, sep = "\t", quote=FALSE)
}
write.table(matrix(data=c(rep(as.numeric(jobID[1]),yout),rep(as.numeric(jobID[2]),yout),
                          rep(j,yout),c(1:yout),L,
                          Bfin[1:yout,],N[1:yout,],rcrt[1:yout,],coexist[1:yout,]),
                   nrow=yout,ncol=(nsp*4+5)),
            file=paste0(outloc,"BfinN",suffix),
            col.names = col.names.BfinN, row.names = FALSE,
            append=TRUE,sep="\t", quote=FALSE)
