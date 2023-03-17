\donttest{ # takes too long for R CMD check
  ## an estimate of the log likelihood:
  ricker() |>
    pfilter(Np=1000) |>
    logLik() |>
    replicate(n=5) -> ll
  logmeanexp(ll)
  ## with standard error:
  logmeanexp(ll,se=TRUE)
  ## with effective sample size
  logmeanexp(ll,ess=TRUE)
}
