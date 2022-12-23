\donttest{
  ## Starting with an existing pomp object

  verhulst() -> po
  
  ## we add or change the initial-state simulator,
  ## rinit, using the 'rinit' argument in any 'pomp'
  ## elementary or estimation function (or in the
  ## 'pomp' constructor itself).
  ## Here, we pass the rinit specification to 'simulate'
  ## as an R function.

  po |>
    simulate(
      rinit=function (n_0, ...) {
        c(n=rpois(n=1,lambda=n_0))
      }
    ) -> sim

  ## We can also pass it as a C snippet:

  po |>
    simulate(
      rinit=Csnippet("n = rpois(n_0);"),
      paramnames="n_0",
      statenames="n"
    ) -> sim

}
