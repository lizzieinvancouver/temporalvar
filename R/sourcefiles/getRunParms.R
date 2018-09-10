#Get run parameters; 
runflag <- ifelse(localflag==1,1,as.numeric(Sys.getenv("PHEN_RUNNUM")))

#jobID: if batch, then jobID & taskID from slurm; if local, randomly gen 999XXXX
jobID <- ifelse(localflag==1,
                c(paste0("999",trunc(runif(1,1000,9999))),"1"),
                Sys.getenv(c("SLURM_ARRAY_JOB_ID","SLURM_ARRAY_TASK_ID")))
print(paste("jobID is ",jobID))

#GET INPUT PARMS FOR THIS RUN
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
megaD <- ifelse(inputs$megaD[runflag]==0,0,1) #megaD is flag for megadrought run
rho <- inputs$megaD[runflag]  #rho is the value of the megaD in getInputParms and is corr(s,phi)

nyrs <- sum(nonsta)  # number of yrs to run if nonsta=0 or for initial period if nonsta>0
yrs <- c(1:nyrs)
ext <- 1/10000  #Extinction Threshold:  1 seed/ha (assuming that initial density is 10 seeds per meter)
ndays <- 5  # number of days in a growing season
dt <- 0.005 # within yr timestep
tsteps <- ndays/dt

runparms <- matrix(data= c(jobID[1],jobID[2],nruns,nsp,nyrs,nonsta,tracking,varRstar,vartauI,
                           megaD,writeBout,ext,ndays,dt,tsteps),nrow=1)

#OUTPUT PARMS & FOLDER LOCATIONS
writeBout <- 1  #default=1; flag indicating how often Bout should be written (0=never, n = every n runs)
suffix <- paste0("_",jobID[1],"-",jobID[2],".txt") #unique for each array in batchfile
datesuffix <- (paste0(format(Sys.time(),"%Y-%m-%d")))

if (localflag==0){
  #locOUT is the fast-writing scratch space on the run node (mem limited)
  ##***I think jobID[1] = Sys.getenv("SLURM_ARRAY_JOB_ID) = Sys.getenv("SLURM_JOB_ID")
  locOUT <- paste0("/scratch/wolkovich_lab/temporalvar/",jobID[1],"/")
  #locSAVE is the permanent location to where runs are saved on the storage node
  locSAVE <- ifelse(megaD==1, 
                    "/n/wolkovich_lab/temporalvar/megadrought/output/",
                    "/n/wolkovich_lab/temporalvar/R/output/")
  #locMegaD in is the location of megadrought envt files from Ben  
  locMegaD <- "/n/wolkovich_lab/temporalvar/megadrought/data/"
} else {
  locOUT <- paste0("C:/Users/Megan/Documents/scratch/",jobID[1],"/")
  locSAVE <- ifelse(megaD==1, 
                    "C:/Users/Megan/Documents/GitHub/temporalvar/megadrought/output/",
                    "C:/Users/Megan/Documents/GitHub/temporalvar/R/output/")
  locMegaD <- "C:/Users/Megan/Documents/GitHub/temporalvar/megadrought/fromBen/"
}
if(!dir.exists(file.path(locOUT))) dir.create(file.path(locOUT),recursive=TRUE)

Bout_loc <- paste0(locOUT,"/Bout/",jobID[1],"/")
if(!dir.exists(file.path(Bout_loc))) dir.create(file.path(Bout_loc),recursive=TRUE)

SummOut_loc <- paste0(locOUT,"/SummaryFiles/",jobID[1],"/")
if(!dir.exists(file.path(SummOut_loc))) dir.create(file.path(SummOut_loc),recursive=TRUE)

OtherOut_loc <- paste0(locOUT,"/OtherOut/",jobID[1],"/")
if(!dir.exists(file.path(OtherOut_loc))) dir.create(file.path(OtherOut_loc),recursive=TRUE)

#WRITE RUN CONDITIONS TO RUNPARMS
#  runparms is identical for every run in a job, so every run in the job has the same filename
#  nonetheless, write this out for every run; it will overwrite, but just in case runs fail
  col.names.runparms <- c("arrayID","taskID","nruns","nsp","nyrs",
                          paste0(rep("nonsta",3),c(1:3)),
                          "tracking",
                          paste0(rep("varRstar",2),c(1:2)),
                          "vartauI","megaD","writeBout","ext","ndays","dt","tsteps")
  fileparms <- paste0(OtherOut_loc,"RunParms_",jobID[1],".txt")
  write.table(runparms,file=fileparms,
              col.names = col.names.runparms, row.names = FALSE,
              sep= "\t",quote=FALSE)