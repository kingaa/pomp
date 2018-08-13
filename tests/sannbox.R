library(pomp)
set.seed(240311270)

ripple <- function (x) {
  r <- sqrt(sum(x^2))
  1-exp(-r^2)*cos(10*r)^2
}

invisible(sannbox(par=c(1),fn=ripple))

f <- function (x) {
  if (sqrt(sum(x^2))<1) ripple(x)
  else Inf
}
invisible(sannbox(par=c(2,2),fn=f))

try(sannbox(par=c(2,2),fn=f,control=list(candidate.dist="bob")))
try(sannbox(par=c(2,2),fn=f,control=list(
  candidate.dist=function(x)x)))
