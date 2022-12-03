options(digits=3)
library(pomp)
set.seed(1147838039)

x <- c(1,1,1,2,2,3,3,3,3,4,5,5,6,6,6)
quantile(x)
wquant(x,weights=rep(1,length(x)))
wquant(c(1,2,3,4,5,6),weights=c(3,2,4,1,2,3))
wquant(c(1,2,3,4,5),c(1,0,0,1,1))
wquant(c(1,2,3,4,5),c(0,1,0,1,1))
wquant(c(1,2,3,4,5),c(0,1,0,1,0))
wquant(c(1,2,3,4,5),c(0,1,0,0,1))
wquant(c(1,1,2,2),c(1,1,1,1))
wquant(c(1,2),c(2,2))

try(wquant(c(1,NA),c(1,2)))
try(wquant(c(1,2),c(NA,1)))
try(wquant(c(1,2,3),c(1,2)))
try(wquant(c(1,2,3),c(1,2,-1)))
try(wquant(c(1,2,3),c(1,1,1),probs=c(0.1,NA)))
try(wquant(c(1,2,3),c(1,2,3),probs=c(0.1,2)))

x <- rnorm(n=10000)
stopifnot(
  all.equal(
    wquant(x,probs=seq(0.1,0.9,by=0.1)),
    quantile(x,probs=seq(0.1,0.9,by=0.1),names=FALSE),
    tolerance=0.01
  )
)
y <- seq(-4,4,by=0.01)
p <- diff(pnorm(y))
y <- 0.5*(head(y,-1)+tail(y,-1))
stopifnot(
  all.equal(
    wquant(y,weights=p,probs=c(0.1,0.5,0.9)),
    qnorm(p=c(0.1,0.5,0.9)),
    tolerance=0.01
  )
)
