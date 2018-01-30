#concatenate output files
#run in regal/temporalvar
#should save to n/wolkovich_lab/temporalvar/R, but currently saves to regal/temporalvar/R/output

#create FIVE output files for each set of runs:
#  RunChars_jobID.txt:  includes runparms
#  Bfin_jobID.txt:  arrayID,runID,yr,Bfin_1,Bfin_2...]
#  StaticParms_jobID.txt:  arrayID,runID,sppvars
#  VaryParms_jobID.txt:  arrayID,runID,yr,R0,tauP,eps,tauI_1,tauI_2,...,g_1,g_2...,
#  Bout_jobID.txt:  arrayID,runID,yr,t,B_1,B_2,...

print(getwd())

nruns <- 100
narrays <- 100
jobID <- 78345799
prefix <- "Track_varR_2spp"

filelocIN <- "R/output/"
filelocOUT <-paste("ModelRuns/",jobID,"/",sep="") #R/output/" #"/n/wolkovich_lab/temporalvar/R/modelruns"
outfile=paste(filelocOUT,prefix,"_",jobID,".csv",sep="")
RunChars_outfile = paste(filelocOUT,"RunChars_",prefix,"_",jobID,".csv",sep=",")
Bfin_outfile = paste(filelocOUT,"RunChars_",prefix,"_",jobID,".csv",sep=",")
StaticParms_outfile = paste(filelocOUT,"StaticParms_",prefix,"_",jobID,".csv",sep=",")
VaryParms_outfile = paste(filelocOUT,"VaryParms_",prefix,"_",jobID,".csv",sep=",")
Bout_outfile = paste(filelocOUT,"Bout_",prefix,"_",jobID,".csv",sep=",")

rem <- 0  #flag to indicate that small files should be deleted after concatenating
Boutflag <-1  #flag indicating which Bout files to concatenate (0=never, n = every n runs)  

modelruns <- list()
for (a in c(1:narrays)){
  for (r in c(1:nruns)){
    filename = paste(filelocIN,prefix,"_",jobID,"-",a,"-run",r,".Rdata",sep="")
    if (file.exists(filename)) {
      load(filename)
      #modelruns[[(a-1)*nruns + r]] <- list(jobID=jobID, arrayNum=a, runNum=r,sppvars=sppvars,
      #                                   envtvars=envtvars, tauIhat=tauIhat, Bfin=Bfin,g=g)
      modelruns <- list(jobID=jobID, arrayNum=a, runNum=r,sppvars=sppvars,
                        envtvars=envtvars, tauIhat=tauIhat, Bfin=Bfin,g=g)
    }
  }
}
save(modelruns,file=paste(filelocOUT,prefix,"_",jobID,".Rdata",sep=""))
rm(modelruns)

modelruns_Bout <- list()  
for (a in c(1:narrays)){
  for (r in c(1:nruns)){
    if (Boutflag >0 && r%%Boutflag==0) {
      filename=paste(filelocIN,prefix,"_Bout_",jobID,"-",a,"-run",r,".Rdata",sep="")
      if (file.exists(filename)) {
        load(filename)  #new syntax
        modelruns_Bout[[(a-1)*nruns + r]] <- list(jobID=jobID, arrayNum=a, runNum=r,Bout=Bout)
        rm(Bout)
      }
    }
  }
}
save(modelruns_Bout,file=paste(filelocOUT,prefix,"_Bout_",jobID,".Rdata",sep=""))
rm(modelruns_Bout)


if (rem==1){
for (a in c(1:narrays)){
  for (r in c(1:nruns)){
    filename1=paste(filelocIN,prefix,"_Bout_",jobID,"-",a,"-run",r,".Rdata",sep="")
    filename2=paste(filelocIN,prefix,"_",jobID,"-",a,"-run",r,".Rdata",sep="")
    if (file.exists(filename1)) { file.remove(filename1) }
    if (file.exists(filename2)) { file.remove(filename2) }  #new syntax
  }
}
}
