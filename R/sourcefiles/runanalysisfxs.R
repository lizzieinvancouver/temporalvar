## Started 19 March 2018 ##
## f(x)s to read in and plot the temporalvar runs ##

getfiles <- function(folderID, filenamestart, numhere, colnameshere){
    filepack <- lapply(numhere, function(numhere) {
    filename <- paste("output/", folderID, "/", filenamestart, numhere, ".txt", sep="")
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

makediffs <- function(df){
    dathere <- df
    dathere$diff.alpha <- dathere$alpha1-dathere$alpha2
    dathere$diff.c <- dathere$c1-dathere$c2
    dathere$diff.rstar <- dathere$Rstar1-dathere$Rstar2
    dathere$diff.gmean <- dathere$g1mean_pre-dathere$g2mean_pre
    dathere$diff.tauIPini <- dathere$tauIPini1_pre-dathere$tauIPini2_pre
    dathere$ratio.rstar <-  dathere$Rstar1/dathere$Rstar2
    dathere$ratio.tauIini <- dathere$tauIPini1_pre/dathere$tauIPini2_pre
    dathere$ratio.alpha <- dathere$alpha1/dathere$alpha2
    return(dathere)
    }

makediffs.ns <- function(df){
    dathere <- df
    dathere$diff.alpha <- dathere$alpha1-dathere$alpha2
    dathere$diff.c <- dathere$c1-dathere$c2
    dathere$diff.rstar <- dathere$Rstar1-dathere$Rstar2
    dathere$diff.gmean <- dathere$g1mean_pre-dathere$g2mean_pre
    dathere$diff.tauIPini_pre <- dathere$tauIPini1_pre-dathere$tauIPini2_pre
    dathere$diff.tauIPns_pre <- dathere$tauIPns1_pre-dathere$tauIPns2_pre
    dathere$diff.tauIPfin_pre <- dathere$tauIPfin1_pre-dathere$tauIPfin2_pre
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


