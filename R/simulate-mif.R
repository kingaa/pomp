## simulate the partially-observed Markov process contained in a mif object, using the fitted parameters
simulate.mif <- function (object, nsim = 1, seed = NULL, ...) {
  p <- particles(object,Np=1,center=coef(object),sd=0)
  simulate(
           as(object,'pomp'),
           nsim=nsim,
           seed=seed,
           params=p,
           ...
           )
}

setMethod('simulate','mif',simulate.mif)
