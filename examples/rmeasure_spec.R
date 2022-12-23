\donttest{
  ## We start with the pre-built Ricker example:
  
  ricker() -> po

  ## To change the measurement model simulator, rmeasure,
  ## we use the 'rmeasure' argument in any 'pomp'
  ## elementary or estimation function.
  ## Here, we pass the rmeasure specification to 'simulate'
  ## as an R function.

  po |>
    simulate(
      rmeasure=function (N, phi, ...) {
        c(y=rpois(n=1,lambda=phi*N))
      }
    ) -> sim

  ## We can also pass it as a C snippet:

  po |>
    simulate(
      rmeasure=Csnippet("y = rpois(phi*N);"),
      paramnames="phi",
      statenames="N"
    ) -> sim

}
