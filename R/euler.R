euler.simulate <- function (xstart, times, params,
                            euler.step.fun, delta.t,
                            statenames = character(0),
                            paramnames = character(0),
                            zeronames = character(0),
                            tcovar, covar, PACKAGE)
{
  if (missing(tcovar))
    tcovar <- range(times)
  if (missing(covar))
    covar <- matrix(nrow=2,ncol=0)
  if (missing(PACKAGE))
    PACKAGE <- ""
  efun <- getNativeSymbolInfo(euler.step.fun,PACKAGE)$address
  .Call(
        euler_model_simulator,
        efun,
        xstart,
        times,
        params,
        delta.t,
        statenames,
        paramnames,
        zeronames,
        tcovar,
        covar
        )
}
