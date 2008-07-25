euler.simulate <- function (xstart, times, params,
                            step.fun, delta.t, ...,
                            statenames = character(0),
                            paramnames = character(0),
                            covarnames = character(0),
                            zeronames = character(0),
                            tcovar, covar, PACKAGE)
{
  if (is.character(step.fun)) {
    efun <- try(
                getNativeSymbolInfo(step.fun,PACKAGE)$address,
                silent=FALSE
                )
    if (inherits(efun,'try-error')) {
      stop("no symbol named ",step.fun," in package ",PACKAGE)
    }
  } else if (is.function(step.fun)) {
    if (!all(c('x','t','params','delta.t','...')%in%names(formals(step.fun))))
      stop(sQuote("step.fun")," must be a function of prototype ",sQuote("step.fun(x,t,params,delta.t,...)"))
    efun <- step.fun
  } else {
    stop(sQuote("step.fun")," must be either a function or the name of a compiled routine")
  }

  .Call(
        euler_model_simulator,
        efun,
        xstart,
        times,
        params,
        delta.t,
        statenames,
        paramnames,
        covarnames,
        zeronames,
        tcovar,
        covar,
        args=pairlist(...)
        )
}

euler.density <- function (x, times, params,
                           dens.fun, ...,
                           statenames = character(0),
                           paramnames = character(0),
                           covarnames = character(0),
                           tcovar, covar, log = FALSE,
                           PACKAGE)
{
  if (is.character(dens.fun)) {
    efun <- try(
                getNativeSymbolInfo(dens.fun,PACKAGE)$address,
                silent=FALSE
                )
    if (inherits(efun,'try-error')) {
      stop("no symbol named ",dens.fun," in package ",PACKAGE)
    }
  } else if (is.function(dens.fun)) {
    if (!all(c('x1','x2','t1','t2','params','...')%in%names(formals(dens.fun))))
      stop(sQuote("dens.fun")," must be a function of prototype ",sQuote("dens.fun(x1,x2,t1,t2,params,...)"))
    efun <- dens.fun
  } else {
    stop(sQuote("dens.fun")," must be either a function or the name of a compiled routine")
  }

  .Call(
        euler_model_density,
        efun,
        x,
        times,
        params,
        statenames,
        paramnames,
        covarnames,
        tcovar,
        covar,
        log,
        args=pairlist(...)
        )
}
