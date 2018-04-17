## Started 19 March 2018 ##
## f(x)s to read in and plot the temporalvar runs ##

getfiles <- function(folderID, filenamestart, numhere, colnameshere){
    filepack <- lapply(numhere, function(numhere) {
    filename <- paste("output/SummaryFiles/", folderID, "/", filenamestart, numhere, ".txt", sep="")
    dat <- read.table(filename, skip=1)
    names(dat) <- colnameshere
    return(data.frame(dat))
    })
    datahere <- do.call("rbind", filepack)
}

getBoutfiles <- function(folderID, filenamestart, numhere, colnameshere){
    filepack <- lapply(numhere, function(numhere) {
    filename <- paste("output/", folderID, "/", filenamestart, numhere, ".txt", sep="")
    dat <- read.table(filename, skip=1)
    names(dat) <- colnameshere
    return(data.frame(dat))
    })
    datahere <- do.call("rbind", filepack)
}

# potentially superior to above as you don't need to tell it...
# ... the numbers to expect, but returns a list
# aka different format than everything else returns
getBoutfiles.list <- function(folderID, boutfilenamestart){
    boutfilenamestart <- 
    filename <- paste("output/", folderID, "/", "Bout/", "Bout*", ".txt", sep="")
    dat <- lapply(Sys.glob(filename), function(i) read.table(i, header=TRUE))
    datahere <- do.call("rbind", dat)
    return(dat)
}

# add tauI! 
makediffs <- function(df){
    dathere <- df
    dathere$c.rstar <-  as.numeric(dathere$c1)/as.numeric(dathere$c2)
    dathere$ratio.rstar <-  as.numeric(dathere$Rstar1)/as.numeric(dathere$Rstar2)
    dathere$ratio.tauIP <- as.numeric(dathere$tauIP1_mean)/as.numeric(dathere$tauIP2_mean)
    dathere$ratio.alpha <- as.numeric(dathere$alpha1)/as.numeric(dathere$alpha2)
    return(dathere)
    }

makediffs.ns <- function(df){
    dathere <- df
    dathere$ratio.rstar <-  dathere$Rstar1/dathere$Rstar2
    dathere$ratio.tauIini_pre <- dathere$tauIPini1_pre/dathere$tauIPini2_pre
    dathere$ratio.tauIPns_pre <-dathere$tauIPns1_pre/dathere$tauIPns2_pre
    dathere$ratio.alpha <- dathere$alpha1/dathere$alpha2
    return(dathere)
    }

####################
## plotting f(x)s ##
####################

plot.params <- function(df, df.coexist, colname, colname.sp1, colname.sp2){
    pdf(paste("graphs/modelruns/params/runs_", folderID, colname, "_compare.pdf", sep=""),
        width=9, height=4)
    plot.allruns <- ggplot(data=df, aes(x=df[colname.sp1], y=df[colname.sp2],
        colour=ncoexist)) +
        geom_point() +
        labs(x = colname.sp1, y=colname.sp2)
    plot.coexistruns <- ggplot(data=df.coexist, aes(x=df.coexist[colname.sp1], 
        y=df.coexist[colname.sp2], colour=ncoexist)) +
        geom_point() +
        labs(x = colname.sp1, y=colname.sp2)
    multiplot(plot.allruns, plot.coexistruns, cols=2)
    dev.off()
}


plot.paramdiffs <- function(df, df.coexist, figname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/runs_", folderID, figname, ".pdf", sep=""),
        width=9, height=4)
    plot.allruns <- ggplot(data=df, aes(x=df[colname.x], y=df[colname.y],
        colour=ncoexist)) +
        geom_point() +
        labs(x = colname.x, y=colname.y)
    plot.coexistruns <- ggplot(data=df.coexist, aes(x=df.coexist[colname.x], 
        y=df.coexist[colname.y], colour=ncoexist)) +
        geom_point() +
        labs(x = colname.x, y=colname.y)
    multiplot(plot.allruns, plot.coexistruns, cols=2)
    dev.off()
}


