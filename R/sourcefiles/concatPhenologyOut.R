#concatenate output files

nruns=10

modelruns <- list()
for (j in c(1:nruns)){
  load(paste("out_",j,".Rdata",sep=""))
  modelruns[[j]]<-list("sppvars"=sppvars, "tauI"=tauI, "tauP"=tauP, "Bfin"=Bfin)
}