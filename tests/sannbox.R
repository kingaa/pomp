options(digits=3)

library(pomp2)
library(magrittr)

set.seed(240311270)

ripple <- function (x) {
  r <- sqrt(sum(x^2))
  1-exp(-r^2)*cos(10*r)^2
}

sannbox(par=c(1),fn=ripple,control=list(lower=1,upper=5)) -> out1
sannbox(par=c(1),fn=ripple,control=list(lower=NULL,upper=NULL)) -> out2
stopifnot(out1$convergence==0, out2$convergence==0)
try(sannbox(par=c(1),fn=ripple,control=list(
  lower=NULL,upper=NULL,
  sched=seq(10,1,length=10))))
sannbox(par=c(1),fn=ripple,control=list(
  lower=1,upper=5,maxit=100,sched=seq(10,1,length=100))) -> out3
stopifnot(out3$convergence == 0)
capture.output(
  sannbox(par=c(1),fn=ripple,control=list(lower=1,upper=5,maxit=10,trace=4)) -> out5
) -> tout
stopifnot(
  sum(grepl("current params",tout))==10,
  sum(grepl("proposed params",tout))==10,
  sum(grepl("accept",tout))==10,
  sum(grepl("best",tout))==1
)

f <- function (x) {
  if (sqrt(sum(x^2))<1) ripple(x)
  else Inf
}
invisible(sannbox(par=c(2,2),fn=f))

try(sannbox(par=c(2,2),fn=f,control=list(candidate.dist="bob")))
try(sannbox(par=c(2,2),fn=f,control=list(
  candidate.dist=function(x)x)))
