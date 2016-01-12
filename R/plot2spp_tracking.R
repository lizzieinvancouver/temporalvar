#plot tauIini vs Rstar with points coded for coexistence


## Helpful reminders on indexing within lists from Lizzie!
# how to call something (given that you have at least 2 runs) ....
modelruns[[2]][1][[1]]$gmax
# how to call internal list by name
modelruns[[2]]["sppvars"]
# when slicing out bits of vector within a list you have to get beyond the $ 
# this is often something I forget so here are two examples
modelruns[[2]]["sppvars"][[1]]$gmax
modelruns[[2]]["tauP"][[1]][1:10]

gtauIPini <- matrix(0,nruns,nsp)
gRstar <- matrix(0,nruns,nsp)

for (i in seq(1,nruns)){
  gtauIPini[i,] <- modelruns[[i]][[1]]$tauIPini
  gRstar[i,] <- modelruns[[i]][[1]]$Rstar
}
plot(gRstar[,1]-gRstar[,2], gtauIPini[,1]-gtauIPini[,2],col=numcoexist+1)
