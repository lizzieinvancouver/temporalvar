#concatenate output files
#run in regal/temporalvar
#should save to n/wolkovich_lab/temporalvar/R, but currently saves to regal/temporalvar/R/output


print(getwd())
filelocIN <- "R/output/"
filelocOUT <-"ModelRuns/" #R/output/" #"/n/wolkovich_lab/temporalvar/R/modelruns"

nruns <- 100
narrays <- 100
jobID <- 78345799
prefix <- "Track_varR_2spp"
rem <- 0  #flag to indicate that small files should be deleted after concatenating
Boutflag <- 25  #flag indicating which Bout files to concatenate (0=never, n = every n runs)  

modelruns <- list()
for (a in c(1:narrays)){
  for (r in c(1:nruns)){
    load(paste(filelocIN,prefix,"_",jobID,"-",a,"-run",r,".Rdata",sep=""))
    #print(paste(filelocIN,prefix,"_",jobID,"-",a,"-run",r,".Rdata",sep=""))
    modelruns[[(a-1)*nruns + r]] <- list(jobID=jobID, arrayNum=a, runNum=r,sppvars=sppvars,
                                         envtvars=envtvars, tauIhat=tauIhat, Bfin=Bfin,g=g)
  }
}
save(modelruns,file=paste(filelocOUT,prefix,"_",jobID,".Rdata",sep=""))
rm(modelruns)

modelruns_Bout <- list()  
for (a in c(1:narrays)){
  for (r in c(1:nruns)){
    if (Boutflag >0 && r%%Boutflag==0) {
      load(paste(filelocIN,prefix,"_Bout_",jobID,"-",a,"-run",r,".Rdata",sep=""))  #new syntax
      #print(paste(filelocIN,prefix,"_Bout_",jobID,"-",a,"-run",r,".Rdata",sep=""))  #new syntax
      modelruns_Bout[[(a-1)*nruns + r]] <- list(jobID=jobID, arrayNum=a, runNum=r,Bout=Bout)
      rm(Bout)
    }
  }
}
save(modelruns_Bout,file=paste(filelocOUT,prefix,"_Bout_",jobID,".Rdata",sep=""))
rm(modelruns)


if (rem==1){
for (a in c(1:narrays)){
  for (r in c(1:nruns)){
    file.remove(paste(filelocIN,prefix,"_",jobID,"-",a,"-run",r,".Rdata",sep=""))
    #file.remove(paste(filelocIN,prefix,"_Bout_",jobID,"-",a,r,".Rdata",sep=""))  #old synatx
    file.remove(paste(filelocIN,prefix,"_Bout_",jobID,"-",a,"-run",r,".Rdata",sep=""))  #new syntax
  }
}
}
