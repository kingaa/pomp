onestep.simulate <- function (xstart, times, params,
                              step.fun, ...,
                              statenames = character(0),
                              paramnames = character(0),
                              covarnames = character(0),
                              zeronames = character(0),
                              tcovar, covar, PACKAGE)
{
  .Deprecated(new="onestep.sim",package="pomp")
  efun <- pomp.fun(
                   f=step.fun,
                   PACKAGE=PACKAGE,
                   proto="step.fun(x,t,params,delta.t,...)"
                   )
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

euler.simulate <- function (xstart, times, params,
                            step.fun, delta.t, ...,
                            statenames = character(0),
                            paramnames = character(0),
                            covarnames = character(0),
                            zeronames = character(0),
                            tcovar, covar, PACKAGE)
{
  .Deprecated(new="euler.sim",package="pomp")
  efun <- pomp.fun(
                   f=step.fun,
                   PACKAGE=PACKAGE,
                   proto="step.fun(x,t,params,delta.t,...)"
                   )
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

onestep.density <- function (x, times, params,
                             dens.fun, ...,
                             statenames = character(0),
                             paramnames = character(0),
                             covarnames = character(0),
                             tcovar, covar, log = FALSE,
                             PACKAGE)
{
  .Deprecated(new="onestep.dens",package="pomp")
  efun <- pomp.fun(
                   f=dens.fun,
                   PACKAGE=PACKAGE,
                   proto="dens.fun(x1,x2,t1,t2,params,...)"
                   )
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
