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

calcbesttauI <- function(df){
    dathere <- df
    dathere.sm <- dathere[c("tauIP1_mean", "tauIP2_mean")]
    bettertauI <- colnames(dathere.sm)[max.col(dathere.sm, ties.method="first")]
    dathere.tauI <- dathere[c("tauI1", "tauI2")]
    dathere.tauI$besttauI <- NA
    dathere.tauI$besttauI[which(bettertauI=="tauIP1_mean")] <-
        dathere.tauI$tauI1[which(bettertauI=="tauIP1_mean")]
    dathere.tauI$besttauI[which(bettertauI=="tauIP2_mean")] <-
        dathere.tauI$tauI2[which(bettertauI=="tauIP2_mean")]
    dathere$besttauI <- dathere.tauI$besttauI
    return(dathere)
    }


calcsp.besttauI <- function(df){
    dathere <- df
    dathere.sm <- dathere[c("tauIP1_mean", "tauIP2_mean")]
    bettertauI <- colnames(dathere.sm)[max.col(dathere.sm, ties.method="first")]
    dathere.tauI <- dathere[c("tauI1", "tauI2")]
    dathere.tauI$sp.besttauI <- NA
    dathere.tauI$sp.besttauI[which(bettertauI=="tauIP1_mean")] <-1
    dathere.tauI$sp.besttauI[which(bettertauI=="tauIP2_mean")] <-2
    dathere$sp.besttauI <- dathere.tauI$sp.besttauI 
    return(dathere)
    }

calcsp.bestrstar <- function(df){
    dathere <- df
    dathere.sm <- dathere[c("Rstar1", "Rstar2")]
    worserstar <- colnames(dathere.sm)[max.col(dathere.sm, ties.method="first")]
    dathere.tauI <- dathere[c("Rstar1", "Rstar2")]
    dathere.tauI$sp.bestrstar <- NA
    dathere.tauI$sp.bestrstar[which(worserstar=="Rstar1")] <-2
    dathere.tauI$sp.bestrstar[which(worserstar=="Rstar2")] <-1
    dathere$sp.bestrstar <- dathere.tauI$sp.bestrstar 
    return(dathere)
    }

calcsp.bestalpha <- function(df){
    dathere <- df
    dathere.sm <- dathere[c("alpha1", "alpha2")]
    bestalpha <- colnames(dathere.sm)[max.col(dathere.sm, ties.method="first")]
    dathere.tauI <- dathere[c("alpha1", "alpha2")]
    dathere.tauI$sp.bestalpha <- NA
    dathere.tauI$sp.bestalpha[which(bestalpha=="alpha1")] <-1
    dathere.tauI$sp.bestalpha[which(bestalpha=="alpha2")] <-2
    dathere$sp.bestalpha <- dathere.tauI$sp.bestalpha 
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

plot.paramdiffs.onepanel <- function(df, figname, runname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/", runname, figname, "1p.pdf", sep=""),
        width=5, height=4)
        par(mfrow=c(1,1))
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


plot.paramdiffs.twopanel <- function(df, figname, runname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/", figname, runname, "2p.pdf", sep=""),
        width=5, height=8)
        par(mfrow=c(2,1))
        df0 <- subset(df, ncoexist.t2==0)
        df1 <- subset(df, ncoexist.t2==1)
        df2 <- subset(df, ncoexist.t2==2)
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="survived after stat: colored by non-stat")
        points(df0[[colname.x]], unlist(df0[colname.y]),
           col=coexist3col[1], pch=16)
        points(df1[[colname.x]], unlist(df1[colname.y]),
           col=coexist3col[2], pch=16)
        points(df2[[colname.x]], unlist(df2[colname.y]),
           col=coexist3col[3], pch=16)
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="survived after nonstat")
        points(df2[[colname.x]], unlist(df2[colname.y]),
           col=coexist3col[3], pch=16)
        legend("topright", leg.txt, pch=16, col=coexist3col, bty="n")
    dev.off()
}

plot.paramdiffs.fixedxy <- function(dfmultivar, figname, runname, colname.x, colname.y,
    dfother){
    pdf(paste("graphs/modelruns/paramdiffs/", figname, runname, "4pXY.pdf", sep=""),
        width=10, height=8)
        par(mfrow=c(2,2))
        df0 <- subset(dfmultivar, ncoexist.t2==0)
        df1 <- subset(dfmultivar, ncoexist.t2==1)
        df2 <- subset(dfmultivar, ncoexist.t2==2)
    xlimhere <- c(min(dfmultivar[[colname.x]]), max(dfmultivar[[colname.x]]))
    ylimhere <- c(min(dfmultivar[[colname.y]]), max(dfmultivar[[colname.y]]))
    # first type of run
        plot(unlist(df1[colname.x]), unlist(df1[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="3 traits vary: survived after stat - color by ns", xlim=xlimhere, ylim=ylimhere)
        points(df0[[colname.x]], unlist(df0[colname.y]),
           col=coexist3col[1], pch=16, xlim=xlimhere, ylim=ylimhere)
        points(df1[[colname.x]], unlist(df1[colname.y]),
           col=coexist3col[2], pch=16, xlim=xlimhere, ylim=ylimhere)
        points(df2[[colname.x]], unlist(df2[colname.y]),
           col=coexist3col[3], pch=16, xlim=xlimhere, ylim=ylimhere)
        legend("topright", leg.txt, pch=16, col=coexist3col, bty="n")
        plot(unlist(dfmultivar[colname.x]), unlist(dfmultivar[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="3 traits vary: survived after nonstat", xlim=xlimhere, ylim=ylimhere)
        points(df2[[colname.x]], unlist(df2[colname.y]),
           col=coexist3col[3], pch=16, xlim=xlimhere, ylim=ylimhere)
    # second type of run
        df20 <- subset(dfother, ncoexist.t2==0)
        df21 <- subset(dfother, ncoexist.t2==1)
        df22 <- subset(dfother, ncoexist.t2==2)
        plot(unlist(dfother[colname.x]), unlist(dfother[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="2 traits vary: survived after stat - color by ns", xlim=xlimhere, ylim=ylimhere)
        points(df20[[colname.x]], unlist(df20[colname.y]),
           col=coexist3col[1], pch=16, xlim=xlimhere, ylim=ylimhere)
        points(df21[[colname.x]], unlist(df21[colname.y]),
           col=coexist3col[2], pch=16, xlim=xlimhere, ylim=ylimhere)
        points(df22[[colname.x]], unlist(df22[colname.y]),
           col=coexist3col[3], pch=16, xlim=xlimhere, ylim=ylimhere)
        legend("topright", leg.txt, pch=16, col=coexist3col, bty="n")
        plot(unlist(dfother[colname.x]), unlist(dfother[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="2 traits vary: survived after nonstat", xlim=xlimhere, ylim=ylimhere)
        points(df22[[colname.x]], unlist(df22[colname.y]),
           col=coexist3col[3], pch=16, xlim=xlimhere, ylim=ylimhere)
    dev.off()
}


plot.paramdiffs.colorbyzvar <- function(df, figname, runname, colname.x, colname.y,
    colname.z, figtitle, midpt){
    plothere <- ggplot(df, aes(df[[colname.x]], df[[colname.y]])) +
        geom_point(aes(color=df[[colname.z]])) +
        scale_colour_gradient2(midpoint = midpt) + 
        labs(colour = colname.z, x = colname.x, y=colname.y,
             title=figtitle)
    # scale_colour_gradient2(low = "white", mid ="white", high = "darkred")
    ggsave(paste("graphs/modelruns/paramdiffs/", figname, runname, ".zvar.pdf", sep=""),
        width=10, height=8)
}

plot.histograms.bothspp <- function(df, figname, colname.x1, colname.x2,
    collist, varcollist, ylim){
    df2here <- subset(df, ncoexist.t2==2)
    pdf(paste("graphs/modelruns/histograms/", figname, "bothspcolorbyzvarp.pdf", sep=""),
        width=4.5, height=4)
    plot(nhere, tauPfin, type="l", ylim=ylim, xlab="tauP and trait", ylab="density", main="trait for both spp")
    polygon(nhere, tauPfin, col=collist[2])
    polygon(nhere, tauP, col=collist[1])
    lines(density(c(df[[colname.x1]], df[[colname.x2]])), col=varcollist[1], lwd=2)
    lines(density(c(df2here[[colname.x1]], df2here[[colname.x2]])), col=varcollist[2], lwd=2)
    dev.off()
}

plot.histograms.max <- function(df, figname, colname.x1, colname.x2,
    collist, varcollist, ylim){
    df2here <- subset(df, ncoexist.t2==2)
    pdf(paste("graphs/modelruns/histograms/", figname, "max.pdf", sep=""),
        width=4.5, height=4)
    plot(nhere, tauPfin, type="l", ylim=ylim, xlab="tauP and trait", ylab="density",  main="max trait")
    polygon(nhere, tauPfin, col=collist[2])
    polygon(nhere, tauP, col=collist[1])
    lines(density(pmax(df[[colname.x1]], df[[colname.x2]])), col=varcollist[1], lwd=2)
    lines(density(pmax(df2here[[colname.x1]], df2here[[colname.x2]])), col=varcollist[2], lwd=2)
    dev.off()
}

# Do we need this one?
plot.histograms.min <- function(df, figname, colname.x1, colname.x2,
    collist, varcollist, ylim){
    df2here <- subset(df, ncoexist.t2==2)
    pdf(paste("graphs/modelruns/histograms/", figname, "min.pdf", sep=""),
        width=4.5, height=4)
    plot(nhere, tauPfin, type="l", ylim=ylim, xlab="tauP and trait", ylab="density",  main="min trait")
    polygon(nhere, tauPfin, col=collist[2])
    polygon(nhere, tauP, col=collist[1])
    lines(density(pmin(df[[colname.x1]], df[[colname.x2]])), col=varcollist[1], lwd=2)
    lines(density(pmin(df2here[[colname.x1]], df2here[[colname.x2]])), col=varcollist[2], lwd=2)
    dev.off()
}

plot.histograms.onespp.skipnonstat <- function(df, figname, colname.x,
    collist, varcollist, ylim=ylim){
    df2here <- subset(df, ncoexist.t2==2)
    pdf(paste("graphs/modelruns/histograms/", figname, "best.pdf", sep=""),
        width=4.5, height=4)
    plot(nhere, tauPfin, type="l", ylim=ylim, xlab="tauP and trait", ylab="density",  main="best tauI")
    polygon(nhere, tauPfin, col=collist[2])
    polygon(nhere, tauP, col=collist[1])
    lines(density(df[[colname.x]]), col=varcollist[1], lwd=2)
    # lines(density(df2here[[colname.x]]), col=varcollist[2], lwd=2)
    dev.off()
}


plot.histograms.onespp <- function(df, figname, colname.x,
    collist, varcollist, ylim=ylim){
    df2here <- subset(df, ncoexist.t2==2)
    pdf(paste("graphs/modelruns/histograms/", figname, "best.pdf", sep=""),
        width=4.5, height=4)
    plot(nhere, tauPfin, type="l", ylim=ylim, xlab="tauP and trait", ylab="density",  main="best tauI")
    polygon(nhere, tauPfin, col=collist[2])
    polygon(nhere, tauP, col=collist[1])
    lines(density(df[[colname.x]]), col=varcollist[1], lwd=2)
    lines(density(df2here[[colname.x]]), col=varcollist[2], lwd=2)
    dev.off()
}


plot.histograms.bars.onespp.skipnonstat <- function(df, figname, colname.x,
    collist, varcollist, breaks){
    df2here <- subset(df, ncoexist.t2==2)
    pdf(paste("graphs/modelruns/histograms/", figname, "barsbest.pdf", sep=""),
        width=4.5, height=4)
    plot(nhere, tauPfin, type="l", xlab="tauP and trait", ylab="density",  main="best tauI")
    polygon(nhere, tauPfin, col=collist[2])
    polygon(nhere, tauP, col=collist[1])
    hist((df[[colname.x]]), col=varcollist[1], breaks=breaks)
    dev.off()
}


## Plots for just the stationary period ###

plot.paramdiffs.stat.onepanel <- function(df, figname, runname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/", runname, figname, ".stat.1p.pdf", sep=""),
        width=5, height=4)
        par(mfrow=c(1,1))
        df0 <- subset(df, ncoexist==0)
        df1 <- subset(df, ncoexist==1)
        df2 <- subset(df, ncoexist==2)
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="Stationary period only")
        points(df0[[colname.x]], unlist(df0[colname.y]),
           col=coexist3col[1], pch=16)
        points(df1[[colname.x]], unlist(df1[colname.y]),
           col=coexist3col[2], pch=16)
        points(df2[[colname.x]], unlist(df2[colname.y]),
           col=coexist3col[3], pch=16)
        legend("topright", leg.txt, pch=16, col=coexist3col, bty="n")
    dev.off()
}

plot.rstar.winnersp.stat <- function(df, figname, runname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/", runname, figname, ".Rstar.winnerofstat.1p.pdf", sep=""),
        width=5, height=4)
        par(mfrow=c(1,1))
        df1 <- subset(df, ncoexist==1)
        df1.sp1wins <- subset(df1, coexist1==1 & sp.bestrstar==1)
        df1.sp2wins <- subset(df1, coexist2==1 & sp.bestrstar==2)
        other1 <- subset(df1, coexist1==1 & sp.bestrstar==2)
        other2 <- subset(df1, coexist2==1 & sp.bestrstar==1)
        otherdf <- rbind(other1, other2)
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="Stationary period only")
        points(df1.sp1wins[[colname.x]], unlist(df1.sp1wins[colname.y]),
           col=coexist3col[1], pch=16)
        points(df1.sp2wins[[colname.x]], unlist(df1.sp2wins[colname.y]),
           col=coexist3col[2], pch=16)
         points(otherdf[[colname.x]], unlist(otherdf[colname.y]),
           col=coexist3col[3], pch=16)
        legend("topright", c("sp1 wins and has better Rstar", "sp2 wins and has better Rstar",
            "other sp wins"), pch=16, col=coexist3col, bty="n")
    dev.off()
}

plot.tauI.winnersp.stat <- function(df, figname, runname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/", runname, figname, ".tauI.winnerofstat.1p.pdf", sep=""),
        width=5, height=4)
        par(mfrow=c(1,1))
        df1 <- subset(df, ncoexist==1)
        df1.sp1wins <- subset(df1, coexist1==1 & sp.besttauI==1)
        df1.sp2wins <- subset(df1, coexist2==1 & sp.besttauI==2)
        other1 <- subset(df1, coexist1==1 & sp.besttauI==2)
        other2 <- subset(df1, coexist2==1 & sp.besttauI==1)
        otherdf <- rbind(other1, other2)
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="Stationary period only")
        points(otherdf[[colname.x]], unlist(otherdf[colname.y]),
           col=coexist3col[3], pch=16)
        points(df1.sp1wins[[colname.x]], unlist(df1.sp1wins[colname.y]),
           col=coexist3col[1], pch=16)
        points(df1.sp2wins[[colname.x]], unlist(df1.sp2wins[colname.y]),
           col=coexist3col[2], pch=16)
        legend("topright", c("sp1 wins and has better tauI", "sp2 wins and has better tauI",
            "other sp wins"), pch=16, col=coexist3col, bty="n")
    dev.off()
}

plot.tauI.winnersp.stat.alt <- function(df, figname, runname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/", runname, figname, ".tauI.winnerofstat.1p.alt.pdf", sep=""),
        width=5, height=4)
        par(mfrow=c(1,1))
        df1 <- subset(df, ncoexist==1)
        df1.sp1wins <- subset(df1, coexist1==1 & sp.besttauI==1)
        df1.sp2wins <- subset(df1, coexist2==1 & sp.besttauI==2)
        other1 <- subset(df1, coexist1==1 & sp.besttauI==2)
        other2 <- subset(df1, coexist2==1 & sp.besttauI==1)
        otherdf <- rbind(other1, other2)
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="Stationary period only")
        points(df1.sp1wins[[colname.x]], unlist(df1.sp1wins[colname.y]),
           col=coexist3col[1], pch=16)
        points(df1.sp2wins[[colname.x]], unlist(df1.sp2wins[colname.y]),
           col=coexist3col[2], pch=16)
        points(otherdf[[colname.x]], unlist(otherdf[colname.y]),
           col=coexist3col[3], pch=16)
        legend("topright", c("sp1 wins and has better tauI", "sp2 wins and has better tauI",
            "other sp wins"), pch=16, col=coexist3col, bty="n")
    dev.off()
}

plot.alpha.winnersp.stat <- function(df, figname, runname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/", runname, figname, ".alpha.winnerofstat.1p.pdf", sep=""),
        width=5, height=4)
        par(mfrow=c(1,1))
        df1 <- subset(df, ncoexist==1)
        df1.sp1wins <- subset(df1, coexist1==1 & sp.bestalpha==1)
        df1.sp2wins <- subset(df1, coexist2==1 & sp.bestalpha==2)
        other1 <- subset(df1, coexist1==1 & sp.bestalpha==2)
        other2 <- subset(df1, coexist2==1 & sp.bestalpha==1)
        otherdf <- rbind(other1, other2)
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="Stationary period only")
        points(df1.sp1wins[[colname.x]], unlist(df1.sp1wins[colname.y]),
           col=coexist3col[1], pch=16)
        points(df1.sp2wins[[colname.x]], unlist(df1.sp2wins[colname.y]),
           col=coexist3col[2], pch=16)
        points(otherdf[[colname.x]], unlist(otherdf[colname.y]),
           col=coexist3col[3], pch=16)
        legend("topright", c("sp1 wins and has better alpha", "sp2 wins and has better alpha",
            "other sp wins"), pch=16, col=coexist3col, bty="n")
    dev.off()
}

plot.alpha.winnersp.stat.alt <- function(df, figname, runname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/", runname, figname, ".alpha.winnerofstat.1p.alt.pdf", sep=""),
        width=5, height=4)
        par(mfrow=c(1,1))
        df1 <- subset(df, ncoexist==1)
        df1.sp1wins <- subset(df1, coexist1==1 & sp.bestalpha==1)
        df1.sp2wins <- subset(df1, coexist2==1 & sp.bestalpha==2)
        other1 <- subset(df1, coexist1==1 & sp.bestalpha==2)
        other2 <- subset(df1, coexist2==1 & sp.bestalpha==1)
        otherdf <- rbind(other1, other2)
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="Stationary period only")
        points(otherdf[[colname.x]], unlist(otherdf[colname.y]),
           col=coexist3col[3], pch=16)
        points(df1.sp1wins[[colname.x]], unlist(df1.sp1wins[colname.y]),
           col=coexist3col[1], pch=16)
        points(df1.sp2wins[[colname.x]], unlist(df1.sp2wins[colname.y]),
           col=coexist3col[2], pch=16)
        legend("topright", c("sp1 wins and has better alpha", "sp2 wins and has better alpha",
            "other sp wins"), pch=16, col=coexist3col, bty="n")
    dev.off()
}

plot.paramdiffs.stat.onepanel <- function(df, figname, runname, colname.x, colname.y){
    pdf(paste("graphs/modelruns/paramdiffs/", runname, figname, ".stat.1p.pdf", sep=""),
        width=5, height=4)
        par(mfrow=c(1,1))
        df0 <- subset(df, ncoexist==0)
        df1 <- subset(df, ncoexist==1)
        df2 <- subset(df, ncoexist==2)
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="Stationary period only")
        points(df0[[colname.x]], unlist(df0[colname.y]),
           col=coexist3col[1], pch=16)
        points(df1[[colname.x]], unlist(df1[colname.y]),
           col=coexist3col[2], pch=16)
        points(df2[[colname.x]], unlist(df2[colname.y]),
           col=coexist3col[3], pch=16)
        legend("topright", leg.txt, pch=16, col=coexist3col, bty="n")
    dev.off()
}



########################################
## old code (pre August 2018 meeting) ##
########################################

if(FALSE){
# This is similar to other plotting code (above)
    # but was designed to run on each run, now we're generally merging runs more
plot.histograms.perrun <- function(df1, df2, figname, colname.x1, colname.x2,
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
    
plot.paramdiffs.perrun <- function(df, figname, colname.x, colname.y){
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

}

############################################
## old code (pre April 2018 trip to Oahu) ##
############################################

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


