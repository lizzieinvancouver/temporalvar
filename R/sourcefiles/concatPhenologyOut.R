#concatenate output files
#run in regal/temporalvar
#should save to n/wolkovich_lab/temporalvar/R, but currently saves to regal/temporalvar/R/output


print(getwd())
filelocIN <- "R/output/"
filelocOUT <-"R/output/" #"/n/wolkovich_lab/temporalvar/R/modelruns"

nruns <- 100
narrays <- 10
jobID <- 78345799
prefix <- "Track_varR_2spp"
rem <- 0  #flag to indicate that small files should be deleted after concatenating
Boutflag <- 1  #flag to indicate whether to include Bout files in modelruns list  

modelruns <- list()
for (a in c(1:narrays)){
  for (r in c(1:nruns)){
    load(paste(filelocIN,prefix,"_",jobID,"-",a,"-run",r,".Rdata",sep=""))
    print(paste(filelocIN,prefix,"_",jobID,"-",a,"-run",r,".Rdata",sep=""))
    if (Boutflag >0) {
      load(paste(filelocIN,prefix,"_Bout_",jobID,"-",a,"-run",r,".Rdata",sep=""))  #new syntax
      modelruns[[(a-1)*nruns + r]] <- list(jobID=jobID, arrayNum=a, runNum=r,sppvars=sppvars,
                                            tauI=tauI, envtvars=envtvars, Bfin=Bfin,Bout=Bout)
    } else {
      modelruns[[(a-1)*nruns + r]] <- list(jobID=jobID, arrayNum=a, runNum=r,sppvars=sppvars,
                                           tauI=tauI, envtvars=envtvars, Bfin=Bfin)
    }
  }
}
save(modelruns,file=paste(filelocOUT,prefix,"_",jobID,".Rdata",sep=""))
rm(modelruns)

if (rem==1){
for (a in c(1:narrays)){
  for (r in c(1:nruns)){
    file.remove(paste(filelocIN,prefix,"_",jobID,"-",a,"-run",r,".Rdata",sep=""))
    file.remove(paste(filelocIN,prefix,"_Bout_",jobID,"-",a,r,".Rdata",sep=""))  #old synatx
    #file.remove(paste(filelocIN,prefix,"_Bout_",jobID,"-",a,"-run",r,".Rdata",sep=""))  #new syntax
  }
}
}
