library(magrittr)
library(pomp)

v <- cbind(death = c(-1,1))

works <- function(times = c(1), t0 = 0, mu = 0.001, N_0 = 1) {
    data <- data.frame(time = times, reports = NA)
    rmeasure <- Csnippet("reports = ct;")
    rprocess <- gillespie.sim(Csnippet("rate = mu * N;"), v = v)
    initializer <- function(params, t0, ...) {
        c(N=N_0,ct=12)
    }
    pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
         rprocess = rprocess, rmeasure = rmeasure, initializer = initializer,
         zeronames="ct", paramnames=c("mu"), statenames=c("N","ct"))
}
mwe <- simulate(works(), seed=1L)

d <- cbind(c(1,0))
try(mwe %>% pomp(rprocess=gillespie.sim(Csnippet("rate = mu * N;"),v=v, d=d),
                 zeronames="ct", paramnames=c("mu"), statenames=c("N","ct")) %>%
    simulate() %>% invisible())

vdup <- cbind(death = c(N=-1,N=1))
try(mwe %>% pomp(rprocess=gillespie.sim(Csnippet("rate = mu * N;"),v=vdup),
                 zeronames="ct", paramnames=c("mu"), statenames=c("N","ct")) %>%
    simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.sim(Csnippet("rate = mu * N;"),v=vdup),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct", "N")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.sim(Csnippet("rate = mu * N;"),v=vdup),
                 zeronames="ct", paramnames=c("mu", "mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.sim(Csnippet("rate = mu * N;"),v=vdup),
                 zeronames=c("ct", "ct"), paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

mwe %>% pomp(rprocess=gillespie.hl.sim(list(Csnippet("rate = mu * N;"), c(N=-1, ct=1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible()

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(3L, c(N=-1, ct=1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(c("rate = mu * N;", "mistake"),
                                                c(N=-1, ct=1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(c("rate = mu * N;"),
                                                c(N="-1", ct=1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(c("rate = mu * N;"), c(N=-1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(c("rate = mu * N;"),
                                                c(N=-1, Ct=1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(c("rate = mu * N;"),
                                                c(-1, 1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

## vnames in different order than statevars
## order of names for each elementary event different
## Cnsippets in hl.sim, .sim, and R function in .sim all are consistent
