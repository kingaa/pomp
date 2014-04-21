setClass(
         "pompPlugin",
         slots=c(
           csnippet='logical',
           slotname='character',
           PACKAGE='character'
           ),
         prototype=prototype(
           csnippet=FALSE,
           slotname=character(0),
           PACKAGE=character(0)
           )
         )

setClass(
         "onestepRprocessPlugin",
         contains="pompPlugin",
         slots=c(
           step.fn="ANY"
           )
         )

setClass(
         "discreteRprocessPlugin",
         contains="pompPlugin",
         slots=c(
           step.fn="ANY",
           delta.t="numeric"
           )
         )
           
setClass(
         "eulerRprocessPlugin",
         contains="pompPlugin",
         slots=c(
           step.fn="ANY",
           delta.t="numeric"
           )
         )
           
setClass(
         "gillespieRprocessPlugin",
         contains="pompPlugin",
         slots=c(
           rate.fn="ANY",
           v="matrix",
           d="matrix"
           )
         )

setClass(
         "onestepDprocessPlugin",
         contains="pompPlugin",
         slots=c(
           dens.fn="ANY"
           )
         )

onestep.sim <- function (step.fun, PACKAGE) {
  if (missing(PACKAGE)) PACKAGE <- character(0)
  new("onestepRprocessPlugin",
      step.fn=step.fun,
      slotname="step.fn",
      PACKAGE=PACKAGE)
}

discrete.time.sim <- function (step.fun, delta.t = 1, PACKAGE) {
  if (missing(PACKAGE)) PACKAGE <- character(0)
  new("discreteRprocessPlugin",
      step.fn=step.fun,delta.t=delta.t,
      slotname="step.fn",
      csnippet=is(step.fun,"Csnippet"),
      PACKAGE=PACKAGE)
}

euler.sim <- function (step.fun, delta.t, PACKAGE) {
  if (missing(PACKAGE)) PACKAGE <- character(0)
  new("eulerRprocessPlugin",
      step.fn=step.fun,delta.t=delta.t,
      slotname="step.fn",
      csnippet=is(step.fun,"Csnippet"),
      PACKAGE=PACKAGE)
}

gillespie.sim <- function (rate.fun, v, d, PACKAGE) {
  if (missing(PACKAGE)) PACKAGE <- character(0)
  new("gillespieRprocessPlugin",
      rate.fn=rate.fun,v=v,d=d,
      slotname="rate.fn",
      csnippet=is(rate.fun,"Csnippet"),
      PACKAGE=PACKAGE)
}

onestep.dens <- function (dens.fun, PACKAGE) {
  if (missing(PACKAGE)) PACKAGE <- character(0)
  new("onestepDprocessPlugin",
      dens.fn=dens.fun,
      slotname="dens.fn",
      csnippet=is(dens.fun,"Csnippet"),
      PACKAGE=PACKAGE)
}

setMethod(
          "plugin.handler",
          signature=signature(object='function'),
          definition=function (object, ...) {
            object
          }
          )

setMethod(
          "plugin.handler",
          signature=signature(object='ANY'),
          definition=function (object, ...) {
            stop("plugin has an invalid form")
          }
          )

setMethod(
          "plugin.handler",
          signature=signature(object='onestepRprocessPlugin'),
          definition=function (object, ...) {
            efun <- pomp.fun(
                             f=object@step.fn,
                             PACKAGE=object@PACKAGE,
                             proto=quote(step.fun(x,t,params,delta.t,...)),
                             slotname=object@slotname,
                             ...
                             )
            function (xstart, times, params, ...,
                      zeronames = character(0),
                      tcovar, covar,
                      .getnativesymbolinfo = TRUE) {
              .Call(
                    euler_model_simulator,
                    func=efun,
                    xstart=xstart,
                    times=times,
                    params=params,
                    dt=0,
                    method=1L,
                    zeronames=zeronames,
                    tcovar=tcovar,
                    covar=covar,
                    args=pairlist(...),
                    gnsi=.getnativesymbolinfo
                    )
            }
          }
          )

setMethod(
          "plugin.handler",
          signature=signature(object='discreteRprocessPlugin'),
          definition=function (object, ...) {
            efun <- pomp.fun(
                             f=object@step.fn,
                             PACKAGE=object@PACKAGE,
                             proto=quote(step.fun(x,t,params,...)),
                             slotname=object@slotname,
                             ...
                             )
            function (xstart, times, params, ...,
                      zeronames = character(0),
                      tcovar, covar,
                      .getnativesymbolinfo = TRUE) {
              .Call(
                    euler_model_simulator,
                    func=efun,
                    xstart=xstart,
                    times=times,
                    params=params,
                    dt=object@delta.t,
                    method=2L,
                    zeronames=zeronames,
                    tcovar=tcovar,
                    covar=covar,
                    args=pairlist(...),
                    gnsi=.getnativesymbolinfo
                    )
            }
          }
          )

setMethod(
          "plugin.handler",
          signature=signature(object='eulerRprocessPlugin'),
          definition=function (object, ...) {
            efun <- pomp.fun(
                             f=object@step.fn,
                             PACKAGE=object@PACKAGE,
                             proto=quote(step.fun(x,t,params,delta.t,...)),
                             slotname=object@slotname,
                             ...
                             )
            function (xstart, times, params, ...,
                      zeronames = character(0),
                      tcovar, covar,
                      .getnativesymbolinfo = TRUE) {
              .Call(
                    euler_model_simulator,
                    func=efun,
                    xstart=xstart,
                    times=times,
                    params=params,
                    dt=object@delta.t,
                    method=0L,
                    zeronames=zeronames,
                    tcovar=tcovar,
                    covar=covar,
                    args=pairlist(...),
                    gnsi=.getnativesymbolinfo
                    )
            }
          }
          )

setMethod(
          "plugin.handler",
          signature=signature(object='gillespieRprocessPlugin'),
          definition=function (object, ...) {
            if (!(is.matrix(object@d)&&is.matrix(object@v))) {
              stop(sQuote("v")," and ",sQuote("d")," must be matrices")
            }
            nvar <- nrow(object@v)
            nevent <- ncol(object@v)
            if ((nvar!=nrow(object@d))||(nevent!=ncol(object@d)))
              stop(sQuote("v")," and ",sQuote("d")," must agree in dimension")
            efun <- pomp.fun(
                             f=object@rate.fn,
                             PACKAGE=object@PACKAGE,
                             proto=quote(rate.fun(j,x,t,params,...)),
                             slotname=object@slotname,
                             ...
                             )
            function (xstart, times, params,
                      zeronames = character(0),
                      tcovar, covar,
                      .getnativesymbolinfo = TRUE,
                      ...) {
              .Call(
                    SSA_simulator,
                    func=efun,
                    mflag=0L, ## Gillespie's algorithm
                    xstart=xstart,
                    times=times,
                    params=params,
                    e=rep(0,nvar),
                    vmatrix=object@v,
                    dmatrix=object@d,
                    tcovar=tcovar,
                    covar=covar,
                    zeronames=zeronames,
                    args=pairlist(...),
                    gnsi=.getnativesymbolinfo
                    )
            }
          }
          )

setMethod(
          "plugin.handler",
          signature=signature(object='onestepDprocessPlugin'),
          definition=function (object, ...) {
            efun <- pomp.fun(
                             f=object@dens.fn,
                             PACKAGE=object@PACKAGE,
                             proto=quote(dens.fun(x1,x2,t1,t2,params,...)),
                             slotname=object@slotname,
                             ...
                             )
            function (x, times, params, ...,
                      tcovar, covar, log = FALSE,
                      .getnativesymbolinfo = TRUE) {
              .Call(
                    euler_model_density,
                    func=efun,
                    x=x,
                    times=times,
                    params=params,
                    tcovar=tcovar,
                    covar=covar,
                    log=log,
                    args=pairlist(...),
                    gnsi=.getnativesymbolinfo
                    )
            }
          }
          )
