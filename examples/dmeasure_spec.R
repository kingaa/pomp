\donttest{
  ## We start with the pre-built Ricker example:
  ricker() -> po

  ## To change the measurement model density, dmeasure,
  ## we use the 'dmeasure' argument in any 'pomp'
  ## elementary or estimation function.
  ## Here, we pass the dmeasure specification to 'pfilter'
  ## as an R function.

  po %>%
    pfilter(
      dmeasure=function (y, N, phi, ..., log) {
        dpois(y,lambda=phi*N,log=log)
      },
      Np=100
    ) -> pf

  ## We can also pass it as a C snippet:

  po %>%
    pfilter(
      dmeasure=Csnippet("lik = dpois(y,phi*N,give_log);"),
      paramnames="phi",
      statenames="N",
      Np=100
    ) -> pf

}
