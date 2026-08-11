rm(list = ls())
gc()
set.seed(7)
library(deSolve)
library(doFuture)
library(foreach)


wd <- '/home/victor/projects/temporalvar/R'

# Simulation settings
nsp    <- 2
nonsta <- c(100, 100, 100)
nyrs   <- sum(nonsta)
yrs    <- 1:nyrs
ext    <- 1e-3
ndays  <- 5
dt     <- 0.005

# Fixed environment
p <- 10; q <- 10
tauP <- rbeta(nyrs, p, q)
# R0   <- c(rnorm(100, 5, 1.5), rnorm(100, 2, 2), rnorm(100, 5, 1.5))
# R0[R0 < 0] <- 0

xDrought <- 0.5
locMegaD <- "/home/victor/projects/temporalvar/megadrought/fromBen"
R0wet1 <- scan(file=paste0(locMegaD,"/data_wetresamp.csv"),skip=(1+trunc(runif(1,0,1)*10000)),nlines=1,sep=",",quiet=TRUE)
R0wet2 <- scan(file=paste0(locMegaD,"/data_wetresamp.csv"),skip=(1+trunc(runif(1,0,1)*10000)),nlines=1,sep=",",quiet=TRUE)
R0dry <- scan(file=paste0(locMegaD,"/data_dryresamp.csv"),skip=(1+trunc(runif(1,0,1)*10000)),nlines=1,sep=",",quiet=TRUE)
# R0 <- c(rep(R0wet[2:length(R0wet)],length.out=nonsta[1]),c(rep(xDrought*R0dry[2:length(R0dry)],length.out=nonsta[2])))
R0 <- c(rep(R0wet1[2:length(R0wet1)],length.out=nonsta[1]),
        rep(xDrought*R0dry[2:length(R0dry)],length.out=nonsta[2]), 
        rep(R0wet2[2:length(R0wet2)],length.out=nonsta[3]))


eps  <- 1

# Fixed pecies parameter
s <- c(0.8, 0.8)
b <- rep(1, nsp)
h <- rep(100, nsp)
phi <- rep(0.25, nsp)
tauI <- c(0.5, 0.5) # could change this
a <- rep(20, nsp)
u <- rep(1, nsp)
theta<- rep(1, nsp)
m <- rep(0.05, nsp)
N0 <- rep(20, nsp)

source(file.path(wd, "sourcefiles/ResCompN.R"), local = TRUE)

# Function to run simulation (sp1 is exotic, sp2 is native)
run_sim <- function(cvec, gvec, ae) {
  gmax  <- gvec[order(-gvec)] # assign higher to exotic
  alpha <- c(ae, 0)
  c <- cvec[order(cvec)] # assign higher to native              
  Rstar <- (m / (a * (c - m * u)))^(1 / theta)
  tauIhat <- matrix(alpha, nyrs, nsp, byrow = TRUE) * tauP + matrix((1 - alpha) * tauI, nyrs, nsp, byrow = TRUE)
  g <- gmax * exp(-h * (matrix(rep(tauP, nsp), nyrs, nsp) - tauIhat)^2)
  
  N <- matrix(0, nyrs, nsp)
  Bfin <- matrix(0, nyrs, nsp)
  Pars <- list(c = c, a = a, u = u, m = m, theta = theta, eps = eps)
  
  for (y in yrs) {
    if (y == 1) {
      N[y, ] <- N0
    } else {
      N[y, ] <- N[y-1, ] * (1 - g[y-1, ]) * s + phi * Bfin[y-1, ] * s
      N[y, ] <- N[y, ] * (N[y, ] > ext)
      if (sum(N[y, ] > 0) == 0) break
    }
    B0y <- b * g[y, ] * N[y, ]
    if (all(B0y <= 0)) { Bfin[y, ] <- 0; next }
    State <- c(R = R0[y], B = B0y)
    Rstarmin <- min(Rstar[N[y, ] != 0])
    rootfun  <- function(Time, State, Pars) State[1] - Rstarmin
    Bout <- as.data.frame(ode(func = ResCompN, y = State, parms = Pars,
                              times = seq(0, ndays, by = dt), rootfun = rootfun))
    Bfin[y, ] <- apply(Bout[3:(2 + nsp)], 2, max)
    Bfin[y, ] <- Bfin[y, ] * (Rstar < R0[y])
  }
  
  return(
    c(N1dry = N[200, 1], N2dry = N[200, 2], # end of dry period
      N1wet = N[100, 1], N2wet = N[100, 2])) # end of first wet period
}

# Just a wraper to classify the output
classify <- function(e, n) {
  if(e >  0 & n <= 0){
    return(1) # exotic wins
  }else if(e <= 0 & n > 0){
    return(2) # native wins 
  }else if(e > 0 & n > 0){
    return(3) # coexistence
  }else{
    return(0) # everything die, sad
  }
}

# Sample some parameters
n <- 50000
alphas <- runif(n, 0, 1)
gmaxs <- lapply(1:n, function(i) runif(2, 0.5, 0.75))
cs <- lapply(1:n, function(i) runif(2, 15, 30))

# Parallelize everything! Henrik Bengtsson you're great
plan(multisession, workers = 20)
cat <- foreach(k = 1:n) %dofuture%{
  r <- run_sim(cs[[k]], gmaxs[[k]], alphas[k])
  list(dry = classify(r["N1dry"], r["N2dry"]), wet = classify(r["N1wet"], r["N2wet"]))
}
plan(sequential)

# Make plots
catwet <- sapply(cat, function(c) c$wet)
catdry <- sapply(cat, function(c) c$dry)

cratio <- sapply(cs, function(i){
  c <- i[order(i)]
  c[1]/c[2]
})
gratio <- sapply(gmaxs, function(i){
  g <- i[order(-i)]
  g[1]/g[2]
})

par(mfrow = c(1,2), mar = c(4,4,1,1), cex.main = 1)
cols <- c('0'='#444444', '1'='#B97C7C', '2'='#487575', '3'='#d9c67a') 
plot(x = NULL, y = NULL, xlim = c(0.5,1), ylim = c(0,1),
     xlab = 'c ratio (exo/nat)', ylab = 'alpha (exo)', main = 'End of wet period',
     xaxs = "i", yaxs = "i")
for(i in 1:n){
  points(x = cratio[i], y = alphas[i], col = cols[as.character(catwet[i])], pch = 20, cex = 0.1)
}

bins <- cut(alphas, seq(0, 1, 0.01))
ys <- seq(0.005, 0.995, 0.01)

edge <- tapply(1:n, bins, function(i)
  quantile( cratio[i][catwet[i] == 3], 0.01))
lines(predict(loess(edge ~ ys)), ys, lwd = 6, col = 'white')
lines(predict(loess(edge ~ ys)), ys, lwd = 1, col = '#d9c67a')

edge <- tapply(1:n, bins, function(i)
  quantile( cratio[i][catwet[i] == 3], 0.99))
lines(predict(loess(edge ~ ys)), ys, lwd = 6, col = 'white')
lines(predict(loess(edge ~ ys)), ys, lwd = 1, col = '#d9c67a')

plot(x = NULL, y = NULL, xlim = c(0.5,1), ylim = c(0,1),
     xlab = 'c ratio (exo/nat)', ylab = 'alpha (exo)', main = 'End of dry period',
     xaxs = "i", yaxs = "i")
for(i in 1:n){
  points(x = cratio[i], y = alphas[i], col = cols[as.character(catdry[i])], pch = 20, cex = 0.1)
}


edge <- tapply(1:n, bins, function(i)
  quantile( cratio[i][catdry[i] == 3], 0.01))
lines(predict(loess(edge ~ ys)), ys, lwd = 6, col = 'white')
lines(predict(loess(edge ~ ys)), ys, lwd = 2, col = '#d9c67a')


edge <- tapply(1:n, bins, function(i)
  quantile( cratio[i][catdry[i] == 3], 0.99))
lines(predict(loess(edge ~ ys)), ys, lwd = 6, col = 'white')
lines(predict(loess(edge ~ ys)), ys, lwd = 2, col = '#d9c67a')

layout(mat = matrix(c(1,2,3), ncol = 3), widths = c(0.5,1,0.5))
par(mar = c(4,4,1,1), cex.main = 1)
plot.new()
plot(x = NULL, y = NULL, xlim = c(0.5,1), ylim = c(0,1),
     xlab = 'c ratio (exo/nat)', ylab = 'alpha (exo)', main = '',
     xaxs = "i", yaxs = "i")

edge <- tapply(1:n, bins, function(i)
  quantile( cratio[i][catwet[i] == 3], 0.01))
line1 <- predict(loess(edge ~ ys))
lines(line1, ys, lwd = 2, col = '#d9c67a')

edge <- tapply(1:n, bins, function(i)
  quantile( cratio[i][catdry[i] == 3], 0.01))
line2 <- predict(loess(edge ~ ys))
lines(line2, ys, lwd = 2, col = '#d9c67a')


polygon(c(line2, rev(line1)),
        c(ys,      rev(ys)),
        density = 4,      # nb de traits par pouce
        angle   = 45,
        col     = "#d9c67a",
        border  = NA,
        lwd = 2)

# abline(v = 0.66, lty = 2)
# abline(v = 0.8, lty = 2)
