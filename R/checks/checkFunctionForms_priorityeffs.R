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
# tauD = day of highest germination for species i
# h = variability of germination timing for species i
# gmax = maximum germination rate for species i
# gtot = total germination for species i  (integral of g from day 0 to day T)
# d is day in germination season d = 1...D (e.g., D=30)
  
  tauD = 5
  D = 30
  h = 10
  d = seq(1,D,1)
  gmax = 1
  g <- gmax*exp(-((d-tauD)^2)/h)  #germination fraction at time T
  plot(d,g,type="l")  
  
  
  #next, need to make tauD, h, and/or gmax a function of envt
  