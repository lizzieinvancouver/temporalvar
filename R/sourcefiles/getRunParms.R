#Runtime Parameters (if running locally, then only use first line of getInputParms)
print(paste("Sys.getenv(PHEN_RUNNUM) = ",Sys.getenv("PHEN_RUNNUM")))
#print(paste("JOB_ID = ", Sys.getenv("SLURM_ARRAY_JOB_ID")))
runflag <- ifelse(Sys.getenv("PHEN_RUNNUM")=="",1,as.numeric(Sys.getenv("PHEN_RUNNUM")))

datesuffix <- (paste0(format(Sys.time(),"%Y-%m-%d")))

#jobID: jobID & taskID if slurm; randomly gen 7digit nubmer starting with 999 if local
if(Sys.getenv("SLURM_ARRAY_JOB_ID")=="") {
  jobID <- c(paste0("999",trunc(runif(1,1000,9999))),"1")
} else {
  jobID <- Sys.getenv(c("SLURM_ARRAY_JOB_ID","SLURM_ARRAY_TASK_ID")) 
}
print(paste("jobID is ",jobID))

#output parms & folder locations
writeBout <- 1  #default=1; flag indicating how often Bout should be written (0=never, n = every n runs)

Bout_loc <- paste0(locOUT,"/Bout/",jobID[1],"/")
if(!dir.exists(file.path(Bout_loc))) dir.create(file.path(Bout_loc),recursive=TRUE)

SummOut_loc <- paste0(locOUT,"/SummaryFiles/",jobID[1],"/")
if(!dir.exists(file.path(SummOut_loc))) dir.create(file.path(SummOut_loc),recursive=TRUE)

OtherOut_loc <- paste0(locOUT,"/OtherOut/",jobID[1],"/")
if(!dir.exists(file.path(OtherOut_loc))) dir.create(file.path(OtherOut_loc),recursive=TRUE)

suffix <- paste0("_",jobID[1],"-",jobID[2],".txt") #unique for each array in batchfile

  inputs <- as.data.frame(read.table(file=paste0(locIN,"getInputParms.txt"),
                                     header=TRUE,stringsAsFactors=FALSE,sep="\t"))
  nruns <- inputs$nruns[runflag]
  nonsta <- as.numeric(unlist(strsplit(inputs$nonsta[runflag],",")))
  tracking <- inputs$tracking[runflag]
  varRstar <- ifelse(is.character(inputs$varRstar),
                     as.numeric(unlist(strsplit(inputs$varRstar[runflag],","))),
                     as.numeric(inputs$varRstar))
  if(length(varRstar)==1) varRstar <- c(varRstar,NA)
  vartauI <-inputs$vartauI[runflag]
  nsp <- inputs$nsp[runflag]

#between year
nyrs <- sum(nonsta)  # number of yrs to run if nonsta=0 or for initial period if nonsta>0

yrs <- c(1:nyrs)
ext <- 1/10000  #Extinction Threshold:  1 seed/ha (assuming that initial density is 10 seeds per meter)

#within year
ndays <- 5  # number of days in a growing season
dt <- 0.005 # within yr timestep
tsteps <- ndays/dt

runparms <- matrix(data= c(jobID[1],jobID[2],nruns,nsp,nyrs,nonsta,tracking,varRstar,vartauI,
                           writeBout,ext,ndays,dt,tsteps),nrow=1)

#write run conditions to RunParms
  col.names.runparms <- c("arrayID","taskID","nruns","nsp","nyrs",
                          paste0(rep("nonsta",3),c(1:3)),
                          "tracking",
                          paste0(rep("varRstar",2),c(1:2)),
                          "vartauI","writeBout","ext","ndays","dt","tsteps")
  fileparms <- paste0(OtherOut_loc,"RunParms_",jobID[1],".txt")
  write.table(runparms,file=fileparms,
              col.names = col.names.runparms, row.names = FALSE,
              sep= "\t",quote=FALSE)
  # if(!file.exists(file.path(fileparms))) {
  #   file.create(file.path(fileparms))
  #   col.names.Table_of_RunParms = col.names.runparms
  # } else {
  #   col.names.Table_of_RunParms = FALSE
  # }
  # write.table(runparms,file=fileparms,
  #             col.names = col.names.Table_of_RunParms,row.names = FALSE,
  #             append = TRUE, sep = "\t", quote=FALSE)
