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

#WRITE SummaryOut files at the end of each run; one line per period
col.names.SummaryOut<-c("jobID","taskID","runID","period","nperiods","yout","itertime","ncoexist",
                        paste0(rep("coexist",nsp),c(1:nsp)),
                        paste0(rep("alpha",nsp),c(1:nsp)),
                        paste0(rep("c",nsp),c(1:nsp)),
                        paste0(rep("Rstar",nsp),c(1:nsp)),
                        paste0(rep("tauI",nsp),c(1:nsp)),
                        paste0(rep("g",nsp),c(1:nsp),rep("mean",nsp)),
                        paste0(rep("tauIP",nsp),c(1:nsp),rep("_mean",nsp)),
                        paste0(c("wet","dry"),rep("ID",2)),"rho",
                        "R0_mean","R0_median","R0_autocor",
                        paste0(rep("Bfin",nsp),c(1:nsp)))

#Calculate summary stats for output at nyrs or yout, if species went extinct
#For runs with stationary and nonstationary periods, summary stats are written at the end of each period
nperiods <- sum(nonsta>0)
if (yout <= nonsta[1]) {  
  nst <- yout
} else if ((yout > nonsta[1]) && (yout <= (nonsta[1]+nonsta[2]))) {
    nst <- c(nonsta[1], yout) 
} else if ((yout > (nonsta[1] + nonsta[2])) && (yout <= sum(nonsta))) {
      if (nonsta[2]== 0) nst <- c(nonsta[1],yout)
      if (nonsta[2] > 0) nst <- c(nonsta[1],nonsta[1] + nonsta[2],yout)
}

#calculate gmean, tauIP, and R0 calcs for each time period
for (q in c(1:length(nst))) {
  ini <- ifelse(q==1, 1, nst[q-1] + 1)  # starting year for this period
  fin <- nst[q]                         # ending year for this period (or yout, if extinct before end of period)
  gmean <- colMeans(gmax*exp(-h*(matrix(rep(tauP[ini:fin],nsp),nrow=(fin-ini+1),ncol=nsp)
                                     -tauIhat[ini:fin,])^2))  #germination fraction in year y
  tauIP <- colMeans(abs(tauP[ini:fin] - tauIhat[ini:fin,]))
  R0mean <- mean(R0[ini:fin])
  R0median <- median(R0[ini:fin])
  R0autocor <- cor(R0[ini:(fin-1)],R0[(ini+1):fin])
  if (!((j==1)&&(q==1))) col.names.SummaryOut <- FALSE  #if run>1 or second period, col.names is FALSE
}
write.table(matrix(data=c(as.numeric(jobID[1]),as.numeric(jobID[2]),j,q,nperiods,fin,itertime,
                          sum(coexist[fin,]),coexist[fin,],alpha,c,Rstar,tauI,
                          gmean,tauIP,R0id,rho,R0mean,R0median,R0autocor, Bfin[fin,]),nrow=1),
            file=paste0(SummOut_loc,"/SummaryOut",suffix),
            col.names = col.names.SummaryOut, row.names = FALSE,
            append=TRUE,sep="\t", quote=FALSE)

#WRITE BOUT: this gives the within-year dynamics for every year; every run gets separate file
# note that y is the counter for years in PhenologyModel.r loop
# therefore, y is the number of years before loop breaks bc extinction
# m is the counter for year in this loop
# j is the run number;
# writeBout is the flag saying how frequently to write this file (0=never, n= every nth time)

L <- rep(0,yout)  #number of timesteps in Bout in year y; this goes in BfinN
col.names.Bout <- c("jobID","taskID","runID","yr","time","R",paste0(rep("B",nsp),c(1:nsp)))


for (m in c(1:yout)) {
  L[m] <- dim(Bout[[m]])[1]
  ##only write headers in first year
  if (writeBout > 0){
    if (j%%writeBout==0) {
      if (m>1) col.names.Bout <- FALSE
      write.table(matrix(data=c(rep(as.numeric(jobID[1]),L[m]),rep(as.numeric(jobID[2]),L[m]),
                            rep(j,L[m]),rep(m,L[m]),as.matrix(Bout[[m]])),
                     nrow=L[m],ncol=4+2+nsp),
              file = paste0(Bout_loc,"/Bout","_",jobID[1],"-",jobID[2],"-",j,".txt"),
              col.names = col.names.Bout,row.names = FALSE,
              append = TRUE, sep = "\t", quote=FALSE)
    }
  }
}

#WRITE BFIN at the end of each run
# in Bfin, include number of timesteps in the ode run for each year (L[m])
col.names.BfinN <- c("jobID","taskID","runID","yr","Ltstp",
                     paste0(rep("Bfin",nsp),c(1:nsp)),
                     paste0(rep("N",nsp),c(1:nsp)),
                     paste0(rep("rcrt",nsp),c(1:nsp)),
                     paste0(rep("coexist",nsp),c(1:nsp)))
write.table(matrix(data=c(rep(as.numeric(jobID[1]),yout),rep(as.numeric(jobID[2]),yout),
                          rep(j,yout),c(1:yout),L,
                          Bfin[1:yout,],N[1:yout,],rcrt[1:yout,],coexist[1:yout,]),
                   nrow=yout,ncol=(nsp*4+5)),
            file=paste0(OtherOut_loc,"/BfinN",suffix),
            col.names = col.names.BfinN, row.names = FALSE,
            append=TRUE,sep="\t", quote=FALSE)
