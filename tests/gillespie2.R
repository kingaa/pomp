
library(pomp)

deprecated_args <- function(times = 1, t0 = 0, mu = 0.001, N_0 = 1) {
    data <- data.frame(time = times, reports = NA)
    v <- cbind(death = c(-1,1))
    d <- cbind(c(1,0))
    rprocess <- gillespie.sim(Csnippet("rate = mu * N;"), v = v, d=d)
    initializer <- function(params, t0, ...) {
        c(N=N_0,ct=12)
    }
    pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
         rprocess = rprocess, initializer = initializer, zeronames="ct",
         paramnames=c("mu"), statenames=c("N","ct"))
}
invisible(deprecated_args())

dup_vnames <- function(times = 1, t0 = 0, mu = 0.001, N_0 = 1) {
    data <- data.frame(time = times, reports = NA)
    v <- cbind(death = c(N=-1,N=1))
    rprocess <- gillespie.sim(Csnippet("rate = mu * N;"), v = v)
    initializer <- function(params, t0, ...) {
        c(N=N_0,ct=12)
    }
    pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
         rprocess = rprocess, initializer = initializer, zeronames="ct",
         paramnames=c("mu"), statenames=c("N","ct"))
}
try(simulate(dup_vnames(), states=TRUE))

dup_statenames <- function(times = 1, t0 = 0, mu = 0.001, N_0 = 1) {
    data <- data.frame(time = times, reports = NA)
    v <- cbind(death = c(N=-1,ct=1))
    rprocess <- gillespie.sim(Csnippet("rate = mu * N;"), v = v)
    initializer <- function(params, t0, ...) {
        c(N=N_0,ct=12)
    }
    pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
         rprocess = rprocess, initializer = initializer, zeronames="ct",
         paramnames=c("mu"), statenames=c("N","ct", "N"))
}
try(simulate(dup_statenames(), states=TRUE))

dup_paramnames <- function(times = 1, t0 = 0, mu = 0.001, N_0 = 1) {
    data <- data.frame(time = times, reports = NA)
    v <- cbind(death = c(N=-1,ct=1))
    rprocess <- gillespie.sim(Csnippet("rate = mu * N;"), v = v)
    initializer <- function(params, t0, ...) {
        c(N=N_0,ct=12)
    }
    pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
         rprocess = rprocess, initializer = initializer, zeronames="ct",
         paramnames=c("mu", "mu"), statenames=c("N","ct"))
}
try(simulate(dup_paramnames(), states=TRUE))

bad_code <- function(times = 1, t0 = 0, mu = 0.001, N_0 = 1) {
    data <- data.frame(time = times, reports = NA)
    rprocess <- gillespie.hl.sim(list(3L, c(N=-1, ct=1)))
    initializer <- function(params, t0, ...) {
        c(N=N_0,ct=12)
    }
    pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
         rprocess = rprocess, initializer = initializer, zeronames="ct",
         paramnames=c("mu"), statenames=c("N","ct"))
}
try(simulate(bad_code(), states=TRUE))

bad_code <- function(times = 1, t0 = 0, mu = 0.001, N_0 = 1) {
    data <- data.frame(time = times, reports = NA)
    rprocess <- gillespie.hl.sim(list(c("rate = mu * N;", "mistake"), c(N=-1, ct=1)))
    initializer <- function(params, t0, ...) {
        c(N=N_0,ct=12)
    }
    pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
         rprocess = rprocess, initializer = initializer, zeronames="ct",
         paramnames=c("mu"), statenames=c("N","ct"))
}
try(simulate(bad_code(), states=TRUE))

bad_changes <- function(times = 1, t0 = 0, mu = 0.001, N_0 = 1) {
    data <- data.frame(time = times, reports = NA)
    rprocess <- gillespie.hl.sim(list(Csnippet("rate = mu * N;"), c(N="-1", ct="1")))
    initializer <- function(params, t0, ...) {
        c(N=N_0,ct=12)
    }
    pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
         rprocess = rprocess, initializer = initializer, zeronames="ct",
         paramnames=c("mu"), statenames=c("N","ct"))
}
try(simulate(bad_changes(), states=TRUE))

mismatch_nvar <- function(times = 1, t0 = 0, mu = 0.001, N_0 = 1) {
    data <- data.frame(time = times, reports = NA)
    rprocess <- gillespie.hl.sim(list(Csnippet("rate = mu * N;"), c(N=-1)))
    initializer <- function(params, t0, ...) {
        c(N=N_0,ct=12)
    }
    pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
         rprocess = rprocess, initializer = initializer, zeronames="ct",
         paramnames=c("mu"), statenames=c("N","ct"))
}
try(simulate(mismatch_nvar(), states=TRUE))

wrong_vnames <- function(times = 1, t0 = 0, mu = 0.001, N_0 = 1) {
    data <- data.frame(time = times, reports = NA)
    rprocess <- gillespie.hl.sim(list(Csnippet("rate = mu * N;"), c(N=-1, Ct=1)))
    initializer <- function(params, t0, ...) {
        c(N=N_0,ct=12)
    }
    pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
         rprocess = rprocess, initializer = initializer, zeronames="ct",
         paramnames=c("mu"), statenames=c("N","ct"))
}
try(simulate(wrong_vnames(), states=TRUE))

no_vnames <- function(times = 1, t0 = 0, mu = 0.001, N_0 = 1) {
    data <- data.frame(time = times, reports = NA)
    rprocess <- gillespie.hl.sim(list(Csnippet("rate = mu * N;"), c(-1, 1)))
    initializer <- function(params, t0, ...) {
        c(N=N_0,ct=12)
    }
    pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
         rprocess = rprocess, initializer = initializer, zeronames="ct",
         paramnames=c("mu"), statenames=c("N","ct"))
}
try(simulate(no_vnames(), states=TRUE))

## vnames in different order than statevars
## order of names for each elementary event different
## Cnsippets in hl.sim, .sim, and R function in .sim all are consistent
