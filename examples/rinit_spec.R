## We set up a trivial process model:

trivial <- function (X, Y, ...) {
  c(X = X+1, Y = Y-1)
}

## We specify \code{rinit} with a function that
## sets state variables X and Y to the values in
## parameters X0, Y0:

f <- function (X0, Y0, ...) {
  c(X = X0, Y = Y0)
}

plot(simulate(times=1:5,t0=0,params=c(X0=3,Y0=-7),
  rinit=f,rprocess=onestep(trivial)))

## A function that depends on covariate P and
## time t0, as well as parameter X0:

g <- function (t0, X0, P, ...) {
  c(X = X0, Y = P + sin(2*pi*t0))
}

plot(simulate(times=1:5,t0=0,params=c(X0=3,Y0=-7),
  covar=covariate_table(t=0:10,P=3:13,times="t"),
  rinit=g,rprocess=onestep(trivial)))
