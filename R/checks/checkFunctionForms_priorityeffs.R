#testing some functional forms for priority effect model

#forcing is growing degree days until germination
  force <- seq(1,100,10)
  
#chilling is days between 4-10 C
  chill <- seq(1,100,10)
  
#environment is a function of both forcing and chilling
#  E = F^e_f * C^e_c creates a nonlinear tradeoff between F and C
  
  e_f1 <- 1.
  e_c1 <- .9
  E1 <- outer(force^e_f1, chill^e_c1,"*")
  contour(force,chill,E1)
  e_f2 <- 0.95
  e_c2 <- 1.0
  E2 <-outer(force^e_f2, chill^e_c2,"*")
  contour(force, chill, E2, add = TRUE, col = "red")
  
  
#create functional forms for reltaionship betweeen:
# tauI = day of highest germination for species i
# h = variability of germination timing for species i
# gmax = maximum germination rate for species i
# gtot = total germination for species i  (integral of g from day 0 to day T)
# d is day in germination season d = 1...D (e.g., D=30)
#for illustration, create 5 "species"  
  D <- 25
  d <- seq(1,D,1)
  nsp=5
  g <- matrix(NA,nrow=nsp, ncol=D)
  tauD = seq(1,20,length=nsp)
  h = c(1,5,2,3,.5)
  gmax = c(1,.2,.8,1.2,1)
  plot(d,rep(max(gmax),D),type="n", ylim = c(0,1.25), ylab = "germination fraction",xlab = "day")
  for (i in seq(1,5,1)) {
    g[i,] <- gmax[i]*exp(-((d-tauI[i])^2)/h[i])  #germination fraction at time T
    lines(d,g[i,],type="l", col=i)  
  }
  #Calculate G = total germination.  This is an important constraint
 G <-  rowSums(g)
  
  #next, need to make tauI, h, and/or gmax a function of envt
  