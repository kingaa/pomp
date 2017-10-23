library(pomp)
set.seed(34996688L)

create_example <- function(times = c(1,2), t0 = 0, mu = 0.001, N_0 = 1,
                           covar = data.frame(x = c(0, 1), time = c(0, 52))) {
  data <- data.frame(time = times, reports = NA)
  mmodel <- reports ~ pois(ct)
  rate.fun <- function(j, x, t, params, covars, ...) {
    switch(j, params["mu"]*x["N"], stop("unrecognized event ",j))
  }
  rprocess <- gillespie.sim(rate.fun = rate.fun, v=rbind(N=-1, ct=1))
  initializer <- function(params, t0, ...) {
    c(N=N_0,ct=12)
  }
  pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
       rprocess = rprocess, initializer = initializer,
       zeronames = "ct", paramnames = "mu", statenames = c("N","ct"),
       measurement.model = mmodel, covar = covar, tcovar = "time")
}

names(simulate(create_example(times = 1), as.data.frame=TRUE)) -> nm1
names(simulate(create_example(times = c(1,2)), as.data.frame=TRUE)) -> nm2
stopifnot(nm1==nm2)
