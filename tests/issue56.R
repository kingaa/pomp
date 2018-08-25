library(pomp)

create_example <- function(times, t0 = 0, mu = 0.001, N_0 = 1) {

  rmeasure <- function (ct, ...) {
    c(y=rpois(n=1,ct))
  }

  rate.fun <- function(j, x, t, params, covars, ...) {
    switch(j, params["mu"]*x["N"], stop("unrecognized event ",j))
  }

  rprocess <- gillespie.sim(rate.fun = rate.fun, v=rbind(N=-1, ct=1))

  rinit <- function(params, t0, ...) c(N=N_0,ct=12)

  simulate(times = times, t0 = t0, params = c(mu=mu),
    rprocess = rprocess, rinit = rinit, rmeasure = rmeasure, zeronames = "ct",
    paramnames = "mu", statenames = c("N","ct"), obsnames = "y",
    covar = covariate_table(x=c(0,1),times=c(0,52)), format = "data.frame")
}

create_example(times = 1) -> x1
create_example(times = c(1,2)) -> x2
create_example(times = 0)-> x3
create_example(times = c(0,1)) -> x4
stopifnot(names(x1)==names(x2),x3[1,]==x4[1,])
