# At the end of the job, move output from scratch space on runnode to storage area

#Bout files
if (writeBout>0){
  from = file.path(Bout_loc,list.files(Bout_loc))
  todir = file.path(paste0(locSAVE,"/Bout/",jobID[1]))
  if(!dir.exists(todir)) dir.create(todir,recursive=TRUE)
  to = file.path(paste0(todir),list.files(Bout_loc))
  file.copy(from,to)
}

#Summary Files
from = file.path(SummOut_loc,list.files(SummOut_loc))
todir = file.path(paste0(locSAVE,"/SummaryFiles/",jobID[1]))
if(!dir.exists(todir)) dir.create(todir,recursive=TRUE)
to = file.path(paste0(todir),list.files(SummOut_loc))
file.copy(from,to)

#OtherOut Files
from = file.path(OtherOut_loc,list.files(OtherOut_loc))
todir = file.path(paste0(locSAVE,"/OtherOut/",jobID[1]))
if(!dir.exists(todir)) dir.create(todir,recursive=TRUE)
to = file.path(paste0(todir),list.files(OtherOut_loc))
file.copy(from,to)

#CoexistPlots Files
from = file.path(PlotOut_loc,list.files(PlotOut_loc))
todir = file.path(paste0(locSAVE,"/CoexistPlots/",jobID[1]))
if(!dir.exists(todir)) dir.create(todir,recursive=TRUE)
to = file.path(paste0(todir),list.files(PlotOut_loc))
file.copy(from,to)
