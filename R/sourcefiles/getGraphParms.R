# set up graphics parameters
colerz <- topo.colors(nsp)
rcol <- "darkslateblue"
linez <- rep(1:6, 100) # enough for 600 species for now
lspbyrs <- 1
lresbyrs <- 2
lwd=2

#years to plot for within-year dynamics
plotyrs<- seq(1, nyrs, by=floor(nyrs/8))
