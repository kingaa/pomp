library(pomp)

set.seed(2130639172L)

reulermultinom(n=5,size=100,rate=c(1,2,3),dt=0.1)
reulermultinom(n=5,size=-3,rate=c(1,2,3),dt=0.1)
reulermultinom(n=5,size=100,rate=c(1,-2,3),dt=0.1)
reulermultinom(n=5,size=100,rate=c(1,NA,3),dt=0.1)
reulermultinom(n=5,size=100.3,rate=c(1,2,3),dt=0.1)
reulermultinom(n=0,size=100,rate=c(1,2,3),dt=0.1)
try(reulermultinom(n=-2,size=100,rate=c(1,2,3),dt=0.1))
reulermultinom(n=1,size=100,rate=c(1,2e400,3),dt=0.1)

