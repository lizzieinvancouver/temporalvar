#Get run parameters; 
#runflag <- ifelse(localflag==1,1,as.numeric(Sys.getenv("PHEN_RUNNUM")))
runflag <- ifelse(localflag==1,101,as.numeric(Sys.getenv("PHEN_RUNNUM"))) #use this for running megaD locally
print(paste0("runflag is ", runflag))
megaD <- ifelse(runflag>100,TRUE,FALSE)  #use PHEN_RUNNUM as a flag for megaD runs
inputline <- ifelse(megaD==1,runflag - 100, runflag)

#jobID: if batch, then jobID & taskID from slurm; if local, randomly gen 999XXXX
if(localflag==1) {
   jobID <- c(paste0("999",trunc(runif(1,1000,9999))),"1")
} else {
   jobID <- c(Sys.getenv("SLURM_ARRAY_JOB_ID"),Sys.getenv("SLURM_ARRAY_TASK_ID"))
}

writeBout <- 0  #default=1; flag indicating how often Bout should be written (0=never, n = every n runs)

print(paste0(c(Sys.getenv("SLURM_ARRAY_JOB_ID"),Sys.getenv("SLURM_ARRAY_TASK_ID"))))
print(paste0("jobID = ", jobID))

#GET INPUT PARMS FOR THIS RUN
inputfile <- ifelse(megaD==1,paste0(locIN,"/getInputParms_megaD.txt"),paste0(locIN,"/getInputParms.txt"))
inputs <- as.data.frame(read.table(file=inputfile,
                                   header=TRUE,stringsAsFactors=FALSE,sep="\t"))
nruns <- inputs$nruns[inputline]
nonsta <- as.numeric(unlist(strsplit(inputs$nonsta[inputline],",")))
tracking <- inputs$tracking[inputline]
R0ns_flag <- inputs$nsR0[inputline]
varRstar <- ifelse(is.character(inputs$varRstar),
                   as.numeric(unlist(strsplit(inputs$varRstar[inputline],","))),
                   as.numeric(inputs$varRstar))
if(length(varRstar)==1) varRstar <- c(varRstar,NA)
if (is.character(inputs$vartauI[inputline])) {
  vartauI <- as.numeric(unlist(strsplit(inputs$vartauI[inputline],",")))
} else{
  vartauI <- as.numeric(inputs$vartauI[inputline])
}

nsp <- inputs$nsp[inputline]
rho <- inputs$megaD[inputline]  #rho is the value of the megaD in getInputParms and is corr(s,phi)
if (is.null(rho)) rho <- 0
if (rho==1) rho <- -0.5 #if megaD is treated as a flag (0/1) in the input file, then give rho standard value of -0.5
if (localflag==1 && !(rho==0)) megaD <-1   #local is special case: megaD is indicated only by getInputParms.

nyrs <- sum(nonsta)  # number of yrs to run if nonsta=0 or for initial period if nonsta>0
yrs <- c(1:nyrs)
ext <- 1/10000  #Extinction Threshold:  1 seed/ha (assuming that initial density is 10 seeds per meter)
ndays <- 5  # number of days in a growing season
dt <- 0.005 # within yr timestep
tsteps <- ndays/dt

runparms <- matrix(data= c(jobID[1],jobID[2],nruns,nsp,nyrs,nonsta,tracking,varRstar,inputs$vartauI[inputline],
                           megaD,rho,writeBout,ext,ndays,dt,tsteps),nrow=1)

#OUTPUT PARMS & FOLDER LOCATIONS
suffix <- paste0("_",jobID[1],"-",jobID[2],".txt") #unique for each array in batchfile
datesuffix <- (paste0(format(Sys.time(),"%Y-%m-%d")))

if (localflag==0){
  #locOUT is the fast-writing scratch space on the run node (mem limited)
  ##***I think jobID[1] = Sys.getenv("SLURM_ARRAY_JOB_ID) = Sys.getenv("SLURM_JOB_ID")
  locOUT <- paste0("/scratch/wolkovich_lab/temporalvar/",jobID[1],"-",jobID[2])
  #locSAVE is the permanent location to where runs are saved on the storage node
  locSAVE <- ifelse(megaD==1, 
                    "/n/wolkovich_lab/temporalvar/megadrought/output",
                    "/n/wolkovich_lab/temporalvar/R/output")
  #locMegaD in is the location of megadrought envt files from Ben  
  locMegaD <- "/n/wolkovich_lab/temporalvar/megadrought/fromBen"
} else {
  locOUT <- paste0("C:/Users/Megan/Documents/scratch/",jobID[1])
  locSAVE <- ifelse(megaD==1, 
                    "C:/Users/Megan/Documents/GitHub/temporalvar/megadrought/output",
                    "C:/Users/Megan/Documents/GitHub/temporalvar/R/output")
  locMegaD <- "C:/Users/Megan/Documents/GitHub/temporalvar/megadrought/fromBen"
}
if(!dir.exists(file.path(locOUT))) dir.create(file.path(locOUT),recursive=TRUE)

if (writeBout>0) {
  Bout_loc <- paste0(locOUT,"/Bout/",jobID[1])
  if(!dir.exists(file.path(Bout_loc))) dir.create(file.path(Bout_loc),recursive=TRUE)
}

SummOut_loc <- paste0(locOUT,"/SummaryFiles/",jobID[1])
if(!dir.exists(file.path(SummOut_loc))) dir.create(file.path(SummOut_loc),recursive=TRUE)

OtherOut_loc <- paste0(locOUT,"/OtherOut/",jobID[1])
if(!dir.exists(file.path(OtherOut_loc))) dir.create(file.path(OtherOut_loc),recursive=TRUE)

#WRITE RUN CONDITIONS TO RUNPARMS
#  runparms is identical for every run in a job, so every run in the job has the same filename
#  nonetheless, write this out for every run; it will overwrite, but just in case runs fail
  col.names.runparms <- c("arrayID","taskID","nruns","nsp","nyrs",
                          paste0(rep("nonsta",3),c(1:3)),
                          "tracking",
                          paste0(rep("varRstar",2),c(1:2)),
                          "vartauI","megaDflag","rho","writeBout","ext","ndays","dt","tsteps")
  fileparms <- paste0(OtherOut_loc,"/RunParms_",jobID[1],".txt")
  write.table(runparms,file=fileparms,
              col.names = col.names.runparms, row.names = FALSE,
              sep= "\t",quote=FALSE)
