#Test script

evar <- Sys.getenv(c("SLURM_ARRAY_JOB_ID","SLURM_ARRAY_TASK_ID"))
#evar <-Sys.getenv(c("HOST","R_PAPERSIZE"))

print(evar)

for (i in c(1:5)){
  #save(evar,file=paste("output/testbatch_host-",evar[1],"_psize-",evar[2],"_iter", i,".Rdata",sep=""))
  save(evar,file=paste("output/testbatch_jobid-",evar[1],"_taskid-",evar[2],"_iter", i,".Rdata",sep="")) 
}