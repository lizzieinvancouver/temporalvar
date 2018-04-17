#Runtime Parameters
paste0("Sys.getenv(PHEN_RUNNUM) = ",Sys.getenv("PHEN_RUNNUM"))
runflag <- ifelse(Sys.getenv("PHEN_RUNNUM")=="",0,as.numeric(Sys.getenv("PHEN_RUNNUM")))
#batch <- 1 #(runflag>0)*1  #flag if reading inputs from getInputParms.txt

#jobID: jobID & taskID if slurm; randomly gen 7digit nubmer starting with 999 if local
if(Sys.getenv("SLURM_ARRAY_JOB_ID")=="") {
  jobID <- c(paste0("999",trunc(runif(1,1000,9999))),"1")
} else {
  jobID <- Sys.getenv(c("SLURM_ARRAY_JOB_ID","SLURM_ARRAY_TASK_ID")) 
}
#output parms
writeBout <- 1  #default=1; flag indicating how often Bout should be written (0=never, n = every n runs)
#define/create output folder locations
Bout_loc <- paste0(loc,"output/Bout/",jobID[1],"/")
if(!dir.exists(file.path(Bout_loc))) dir.create(file.path(Bout_loc))

SummOut_loc <- paste0(loc,"output/SummaryFiles/",jobID[1],"/")
if(!dir.exists(file.path(SummOut_loc))) dir.create(file.path(SummOut_loc))

OtherOut_loc <- paste0(loc,"output/OtherOut/",jobID[1],"/")
if(!dir.exists(file.path(OtherOut_loc))) dir.create(file.path(OtherOut_loc))

suffix <- paste0("_",jobID[1],"-",jobID[2],".txt") #unique for each array in batchfile

sink(paste0(OtherOut_loc,"sink_",jobID[1],"-",jobID[2],".Rout"))

# if (batch==0){  #default run, not a batch process
#   nruns <- 2
#   nonsta <- c(200,0,0)  #number of [1] initial stationary,[2]nonstationary,[3]final nonstationary years
#   tracking <- 1         #tracking in these runs?
#   varRstar <- c(1,NA)   #flag for variation in Rstar; 
#                         #c(1,NA) draw c randomly for each spp so R* varies by species; 
#                         #c(0,NA) draw c randomly then assign to all species, so R* varies by run but equal bt spp; 
#                         #c(-1,NA) make c constant and equal to default (c=12 for all spp)
#                         #c(x,x) if varRstar is vector then it gives the specific values of c for each spp (only available for 2 spp)
#   vartauI <-1           #flag that indicates that tauI should vary (1) or be the same (0) between species; if vartauI is vector, then it is giving the tauI values for each species.
#   nsp <- 2
# 
# } else {
  inputs <- as.data.frame(read.table(file=paste0(loc,"getInputParms.txt"),
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
# }

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
col.names.runparms <- c("arrayID","taskID","nruns","nsp","nyrs",
                        paste0(rep("nonsta",3),c(1:3)),
                        "tracking",
                        paste0(rep("varRstar",2),c(1:2)),
                        "vartauI","writeBout","ext","ndays","dt","tsteps")
# write.table(runparms,
#             file=paste0(OtherOut_loc,"/RunParms",jobID[1],"/",suffix),
#             col.names = col.names.runparms, row.names = FALSE,
#             append = FALSE, sep= "\t",quote=FALSE)

#write run conditions to Table of RunParms
if(!file.exists(file.path(paste0(loc,"output/Table_of_RunParms.txt")))) {
  file.create(file.path(paste0(loc,"output/Table_of_RunParms.txt")))
  col.names.Table_of_RunParms = col.names.runparms
} else {
  col.names.Table_of_RunParms = FALSE
}
write.table(runparms,file=paste0(loc,"output/Table_of_RunParms.txt"),
            col.names = col.names.Table_of_RunParms,row.names = FALSE,
            append = TRUE, sep = "\t", quote=FALSE)
  
