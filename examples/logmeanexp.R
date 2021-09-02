\dontrun{
  ## an estimate of the log likelihood:
  po <- ricker()
  ll <- replicate(n=5,logLik(pfilter(po,Np=1000)))
  logmeanexp(ll)
  ## with standard error:
  logmeanexp(ll,se=TRUE)
}
