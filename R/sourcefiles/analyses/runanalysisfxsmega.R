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
    dathere$diff.s <-  dathere$s1-dathere$s2
    dathere$diff.phi <- dathere$phi1-dathere$phi2
    dathere$ratio.s <-  dathere$s1/dathere$s2
    dathere$ratio.phi <- dathere$phi1/dathere$phi2
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

# below taken from:
# https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
 
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
 
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
 
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
 # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
 
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
 
  x1 <- switch(pos,
    topleft     =x[1] + sw, 
    left        =x[1] + sw,
    bottomleft  =x[1] + sw,
    top         =(x[1] + x[2])/2,
    center      =(x[1] + x[2])/2,
    bottom      =(x[1] + x[2])/2,
    topright    =x[2] - sw,
    right       =x[2] - sw,
    bottomright =x[2] - sw)
 
  y1 <- switch(pos,
    topleft     =y[2] - sh,
    top         =y[2] - sh,
    topright    =y[2] - sh,
    left        =(y[1] + y[2])/2,
    center      =(y[1] + y[2])/2,
    right       =(y[1] + y[2])/2,
    bottomleft  =y[1] + sh,
    bottom      =y[1] + sh,
    bottomright =y[1] + sh)
    
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
 
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}



plot.paramdiffs.twopanel <- function(df, runname, figname, colname.x, colname.y, cex, pch,
    corner1.text, corner1.pos, corner2.text, corner2.pos){
    pdf(paste("graphs/megadroughts/paramdiffs/", runname, figname, "2p.pdf", sep=""),
        width=5, height=8)
        par(mfrow=c(2,1))
        df0 <- subset(df, ncoexist.t2==0)
        df1 <- subset(df, ncoexist.t2==1)
        df2 <- subset(df, ncoexist.t2==2)
        df1.sp1 <- subset(df1, coexist1.t2==1)
        df1.sp2 <- subset(df1, coexist2.t2==1)
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="survived after stat: colored by non-stat")
        abline(v=1)
        abline(h=1)
        fig_label(text=corner1.text, region="plot", pos=corner1.pos)
        fig_label(text=corner2.text, region="plot", pos=corner2.pos)
        points(df0[[colname.x]], unlist(df0[colname.y]),
           col=coexist3col[1],pch=pch, cex=cex)
        points(df1[[colname.x]], unlist(df1[colname.y]),
           col=coexist3col[2],pch=pch, cex=cex)
        points(df2[[colname.x]], unlist(df2[colname.y]),
           col=coexist3col[3],pch=pch, cex=cex)
        legend("topright", leg.txt, pch=pch, col=coexist3col, bty="n")
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="survived after nonstat")
        abline(v=1)
        abline(h=1)
        fig_label(text=corner1.text, region="plot", pos=corner1.pos)
        fig_label(text=corner2.text, region="plot", pos=corner2.pos)
        points(df2[[colname.x]], unlist(df2[colname.y]),
           col=coexistmocol[3], pch=pch, cex=cex)
        points(df1.sp1[[colname.x]], unlist(df1.sp1[colname.y]),
           col=coexistmocol[4], pch=pch, cex=cex)
        points(df1.sp2[[colname.x]], unlist(df1.sp2[colname.y]),
           col=coexistmocol[5], pch=pch, cex=cex)
        legend("topright", c("both survived", "sp1 left", "sp2 left"), pch=pch, col=coexistmocol[3:5], bty="n")
    dev.off()
}

plot.paramdiffs.twopanel.fixedxy <- function(df, runname, figname, colname.x, colname.y, cex, pch, xlim, ylim,
    corner1.text, corner1.pos, corner2.text, corner2.pos){
    pdf(paste("graphs/megadroughts/paramdiffs/", runname, figname, "2pXY.pdf", sep=""),
        width=5, height=8)
        par(mfrow=c(2,1))
        df0 <- subset(df, ncoexist.t2==0)
        df1 <- subset(df, ncoexist.t2==1)
        df2 <- subset(df, ncoexist.t2==2)
        df1.sp1 <- subset(df1, coexist1.t2==1)
        df1.sp2 <- subset(df1, coexist2.t2==1)
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="survived after stat: colored by non-stat", xlim=xlim, ylim=ylim)
        abline(v=1)
        abline(h=1)
        fig_label(text=corner1.text, region="plot", pos=corner1.pos)
        fig_label(text=corner2.text, region="plot", pos=corner2.pos)
        points(df0[[colname.x]], unlist(df0[colname.y]),
           col=coexist3col[1],pch=pch, cex=cex, xlim=xlim, ylim=ylim)
        points(df1[[colname.x]], unlist(df1[colname.y]),
           col=coexist3col[2],pch=pch, cex=cex, xlim=xlim, ylim=ylim)
        points(df2[[colname.x]], unlist(df2[colname.y]),
           col=coexist3col[3],pch=pch, cex=cex, xlim=xlim, ylim=ylim)
        legend("topright", leg.txt, pch=pch, col=coexist3col, bty="n")
        plot(unlist(df[colname.x]), unlist(df[colname.y]), type="n", xlab=colname.x,
           ylab=colname.y, main="survived after nonstat", xlim=xlim, ylim=ylim)
        abline(v=1)
        abline(h=1)
        fig_label(text=corner1.text, region="plot", pos=corner1.pos)
        fig_label(text=corner2.text, region="plot", pos=corner2.pos)
        points(df2[[colname.x]], unlist(df2[colname.y]),
           col=coexistmocol[3], pch=pch, cex=cex, xlim=xlim, ylim=ylim)
        points(df1.sp1[[colname.x]], unlist(df1.sp1[colname.y]),
           col=coexistmocol[4], pch=pch, cex=cex, xlim=xlim, ylim=ylim)
        points(df1.sp2[[colname.x]], unlist(df1.sp2[colname.y]),
           col=coexistmocol[5], pch=pch, cex=cex, xlim=xlim, ylim=ylim)
        legend("topright", c("both survived", "sp1 left", "sp2 left"), pch=pch, col=coexistmocol[3:5], bty="n")
    dev.off()
}
