## Started 30 December 2018 ##
## f(x)s to read in and plot the temporalvar runs for the MEGAdrought runs ##
## based off runanalysistxs.R ##

########################################
## f(x)s for reading and manipulating ##
########################################

getfiles <- function(folderID, file.names, colnameshere){
    # numhere <- as.numeric(gsub("^[^-]*-([^.]+).*", "\\1", file.names))
    filepack <- lapply(file.names, function(file.names) {
    filename <- paste("output/megadroughts/SummaryFiles/", folderID, "/", file.names, sep="")
    dat <- read.table(filename, skip=1)
    names(dat) <- colnameshere
    return(data.frame(dat))
    })
    datahere <- do.call("rbind", filepack)
}


makediffs <- function(df){
    dathere <- df
    dathere$diff.c <-  dathere$c1-dathere$c2
    dathere$diff.rstar <-  dathere$Rstar1-dathere$Rstar2
    dathere$diff.tauI <-  dathere$tauI1-dathere$tauI2
    dathere$diff.tauIP <- dathere$tauIP1_mean-dathere$tauIP2_mean
    dathere$diff.alpha <-dathere$alpha1-dathere$alpha2
    dathere$ratio.c <-  dathere$c1/dathere$c2
    dathere$ratio.rstar <-  dathere$Rstar1/dathere$Rstar2
    dathere$ratio.tauI <-  dathere$tauI1/dathere$tauI2
    dathere$ratio.tauIP <- dathere$tauIP1_mean/dathere$tauIP2_mean
    dathere$ratio.alpha <-dathere$alpha1/dathere$alpha2
    return(dathere)
    }
