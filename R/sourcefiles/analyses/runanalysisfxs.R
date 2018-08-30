## Started 19 March 2018 ##
## f(x)s to read in and plot the temporalvar runs ##

getfiles <- function(folderID, file.names, colnameshere){
    # numhere <- as.numeric(gsub("^[^-]*-([^.]+).*", "\\1", file.names))
    filepack <- lapply(file.names, function(file.names) {
    filename <- paste("output/SummaryFiles/", folderID, "/", file.names, sep="")
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


####################
## plotting f(x)s ##
####################

add.alpha <- function(col, alpha=1){ # Stolen from Mage's blog
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

plot.paramdiffs <- function(df, figname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/run_", folderID, figname, ".pdf", sep=""),
        width=5, height=4)
        df0 <- subset(df, ncoexist.t2==0)
        df1 <- subset(df, ncoexist.t2==1)
        df2 <- subset(df, ncoexist.t2==2)
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="")
        points(df0[[colname.x]], unlist(df0[colname.y]),
           col=coexist3col[1], pch=16)
        points(df1[[colname.x]], unlist(df1[colname.y]),
           col=coexist3col[2], pch=16)
        points(df2[[colname.x]], unlist(df2[colname.y]),
           col=coexist3col[3], pch=16)
        legend("topright", leg.txt, pch=16, col=coexist3col, bty="n")
    dev.off()
}

plot.paramdiffs.yourownrunanme <- function(df, figname, runname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/", runname, folderID, figname, ".pdf", sep=""),
        width=5, height=4)
        df0 <- subset(df, ncoexist.t2==0)
        df1 <- subset(df, ncoexist.t2==1)
        df2 <- subset(df, ncoexist.t2==2)
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="")
        points(df0[[colname.x]], unlist(df0[colname.y]),
           col=coexist3col[1], pch=16)
        points(df1[[colname.x]], unlist(df1[colname.y]),
           col=coexist3col[2], pch=16)
        points(df2[[colname.x]], unlist(df2[colname.y]),
           col=coexist3col[3], pch=16)
        legend("topright", leg.txt, pch=16, col=coexist3col, bty="n")
    dev.off()
}


plot.histograms <- function(df1, df2, figname, colname.x1, colname.x2,
    collist, breaknum, xlim, ylim){
    pdf(paste("graphs/modelruns/histograms/run_", folderID, figname, ".pdf", sep=""),
        width=5, height=4)
        hist(c(unlist(df1[colname.x1]), unlist(df1[colname.x2])), xlim=xlim, ylim=ylim,
            breaks=breaknum, col=collist[3], main="", xlab=colname.x1)
        par(new=TRUE)
        hist(c(unlist(df2[colname.x1]), unlist(df2[colname.x2])), xlim=xlim, ylim=ylim,
            breaks=breaknum, col=collist[1], main="", xlab="", ylab="") 
    dev.off()
}

#######################################
## old code (pre April trip to Oahu) ##
#######################################

if(FALSE){
    
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

    
makediffs.ns <- function(df){
    dathere <- df
    dathere$ratio.rstar <-  dathere$Rstar1/dathere$Rstar2
    dathere$ratio.tauIini_pre <- dathere$tauIPini1_pre/dathere$tauIPini2_pre
    dathere$ratio.tauIPns_pre <-dathere$tauIPns1_pre/dathere$tauIPns2_pre
    dathere$ratio.alpha <- dathere$alpha1/dathere$alpha2
    return(dathere)
    }

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
}


