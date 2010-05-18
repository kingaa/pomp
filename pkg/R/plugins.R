onestep.sim <- function (step.fun, PACKAGE) {
  efun <- pomp.fun(
                   f=step.fun,
                   PACKAGE=PACKAGE,
                   proto="step.fun(x,t,params,delta.t,...)"
                   )
  function (xstart, times, params, ...,
            statenames = character(0),
            paramnames = character(0),
            covarnames = character(0),
            zeronames = character(0),
            tcovar, covar) {
    .Call(
          euler_model_simulator,
          func=efun,
          xstart=xstart,
          times=times,
          params=params,
          dt=0,
          method=1L,
          statenames=statenames,
          paramnames=paramnames,
          covarnames=covarnames,
          zeronames=zeronames,
          tcovar=tcovar,
          covar=covar,
          args=pairlist(...)
          )
  }
}

discrete.time.sim <- function (step.fun, delta.t = 1, PACKAGE) {
  efun <- pomp.fun(
                   f=step.fun,
                   PACKAGE=PACKAGE,
                   proto="step.fun(x,t,params,delta.t,...)"
                   )
  function (xstart, times, params, ...,
            statenames = character(0),
            paramnames = character(0),
            covarnames = character(0),
            zeronames = character(0),
            tcovar, covar) {
    .Call(
          euler_model_simulator,
          func=efun,
          xstart=xstart,
          times=times,
          params=params,
          dt=delta.t,
          method=0L,
          statenames=statenames,
          paramnames=paramnames,
          covarnames=covarnames,
          zeronames=zeronames,
          tcovar=tcovar,
          covar=covar,
          args=pairlist(...)
          )
  }
}

euler.sim <- function (step.fun, delta.t, PACKAGE) {
  efun <- pomp.fun(
                   f=step.fun,
                   PACKAGE=PACKAGE,
                   proto="step.fun(x,t,params,delta.t,...)"
                   )
  function (xstart, times, params, ...,
            statenames = character(0),
            paramnames = character(0),
            covarnames = character(0),
            zeronames = character(0),
            tcovar, covar) {
    .Call(
          euler_model_simulator,
          func=efun,
          xstart=xstart,
          times=times,
          params=params,
          dt=delta.t,
          method=0L,
          statenames=statenames,
          paramnames=paramnames,
          covarnames=covarnames,
          zeronames=zeronames,
          tcovar=tcovar,
          covar=covar,
          args=pairlist(...)
          )
  }
}

onestep.dens <- function (dens.fun, PACKAGE) {
  efun <- pomp.fun(
                   f=dens.fun,
                   PACKAGE=PACKAGE,
                   proto="dens.fun(x1,x2,t1,t2,params,...)"
                   )
  function (x, times, params, ...,
            statenames = character(0),
            paramnames = character(0),
            covarnames = character(0),
            tcovar, covar, log = FALSE) {
    .Call(
          euler_model_density,
          func=efun,
          x=x,
          times=times,
          params=params,
          statenames=statenames,
          paramnames=paramnames,
          covarnames=covarnames,
          tcovar=tcovar,
          covar=covar,
          log=log,
          args=pairlist(...)
          )
  }
}

gillespie.sim <- function (rate.fun, v, d, PACKAGE) {
  if (!(is.matrix(d)&&is.matrix(v))) {
    stop(sQuote("v")," and ",sQuote("d")," must be matrices")
  }
  nvar <- nrow(v)
  nevent <- ncol(v)
  if ((nvar!=nrow(d))||(nevent!=ncol(d)))
    stop(sQuote("v")," and ",sQuote("d")," must agree in dimension")

  efun <- pomp.fun(
                   f=rate.fun,
                   PACKAGE=PACKAGE,
                   proto="rate.fun(j,x,t,params,...)"
                   )
  function (xstart, times, params,
            statenames = character(0),
            paramnames = character(0),
            covarnames = character(0),
            zeronames = character(0),
            tcovar, covar, ...) {
    .Call(
          SSA_simulator,
          func=efun,
          mflag=0L, ## Gillespie's algorithm
          xstart=xstart,
          times=times,
          params=params,
          e=rep(0,nvar),
          vmatrix=v,
          dmatrix=d,
          tcovar=tcovar,
          covar=covar,
          statenames=statenames,
          paramnames=paramnames,
          covarnames=covarnames,
          zeronames=zeronames,
          args=pairlist(...)
          )
  }
}
