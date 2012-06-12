reulermultinom <- function (n = 1, size, rate, dt)
  .Call(R_Euler_Multinom,n,size,rate,dt)

deulermultinom <- function (x, size, rate, dt, log = FALSE)
  .Call(D_Euler_Multinom,as.matrix(x),size,rate,dt,log)

rgammawn <- function (n = 1, sigma, dt)
  .Call(R_GammaWN,n,sigma,dt)
