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
    hmax="numeric",
    v="matrix"
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
    csnippet=is(step.fun,"Csnippet"),
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

gillespie.sim <- function (rate.fun, v, d, hmax = Inf, PACKAGE) {
  ep <- paste0("in ",sQuote("gillespie.sim")," plugin: ")
  if (missing(PACKAGE)) PACKAGE <- character(0)
  if (!missing(d)){
    warning(ep, "argument", sQuote("d"), "is deprecated; updates to the simulation",
            "algorithm have made it unnecessary", call. = FALSE)
  }
  if (!is.matrix(v)) {
    stop(ep,sQuote("v")," must be a matrix.",
      call.=FALSE)
  }
  new("gillespieRprocessPlugin",
    rate.fn=rate.fun,v=v, hmax=hmax,
    slotname="rate.fn",
    csnippet=is(rate.fun,"Csnippet"),
    PACKAGE=PACKAGE)
}

gillespie.snip.sim <- function(..., _pre = "", _post = "", hmax = Inf){
    ep <- paste0("in ",sQuote("gillespie.snip.sim")," plugin: ")
    PACKAGE <- character(0) # TODO need this?
    args <- list(...)
    if (anyDuplicated(names(args))) { # TODO does this checker serve any purpose?
        stop(ep,"event arguments must have unique names",call.=FALSE)
    }
    code <- lapply(args, "[[", 1)
    codecheck <- function(x) {
    if(!inherits(x, what = c("Csnippet", "character"))) {
        stop(ep,"the first list element of each event argument should be a"
             " Csnippet or string", call.=FALSE)
    }
    if (length(x) != 1){
        stop(ep,"the length of the first list element of each event",
             "argument should be equal to 1", call.=FALSE)
    }
    lapply(code, codecheck)
    codecheck(_pre)
    codecheck(_post)
    stoich <- lapply(args, "[[", 2)
    stoichcheck <- function(x){
        if (! typeof(x) %in%  c("integer", "double")){
            stop(ep,"the second list element of each event argument should be",
                 "a numeric or integer vector", call.=FALSE)
        }
    }
    lapply(stoich, stoichcheck)
    ## By coercing the vectors to a data frame and then using rbind,
    ## we can ensure that all stochiometric coefficients for the same
    ## state variables are in the same column even if the vectors in
    ## stoich have differently ordered names. Also, rbind will fail if
    ## the set of variables in each data frame is not the same.
    stoichdf <- sapply(stoich, function(x) data.frame(as.list(x)))
    v <- do.call(rbind, stoichdf)
    new("gillespieRprocessPlugin",
      rate.fn=code,v=v, hmax=hmax, # TODO need to check names of v for match to statevars in plugin handler
      slotname="rate.fn",
      csnippet=is(rate.fun,"Csnippet"), # TODO need this?
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
  definition=function (object, purpose = "\b", ...) {
    stop(purpose," plugin has an invalid form.",call.=FALSE)
  }
)

setMethod(
  "plugin.handler",
  signature=signature(object='onestepRprocessPlugin'),
  definition=function (object, ...) {
    ep <- paste0("in ",sQuote("onestep.sim")," plugin: ")
    efun <- tryCatch(
      pomp.fun(
        f=object@step.fn,
        PACKAGE=object@PACKAGE,
        proto=quote(step.fun(x,t,params,delta.t,...)),
        slotname=object@slotname,
        ...
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
    function (xstart, times, params, ...,
      zeronames = character(0),
      tcovar, covar,
      .getnativesymbolinfo = TRUE) {
      tryCatch(
        .Call(
          euler_model_simulator,
          func=efun,
          xstart=xstart,
          times=times,
          params=params,
          deltat=1,
          method=1L,
          zeronames=zeronames,
          tcovar=tcovar,
          covar=covar,
          args=pairlist(...),
          gnsi=.getnativesymbolinfo
        ),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE)
        }
      )
    }
  }
)

setMethod(
  "plugin.handler",
  signature=signature(object='eulerRprocessPlugin'),
  definition=function (object, ...) {
    ep <- paste0("in ",sQuote("euler.sim")," plugin: ")
    efun <- tryCatch(
      pomp.fun(
        f=object@step.fn,
        PACKAGE=object@PACKAGE,
        proto=quote(step.fun(x,t,params,delta.t,...)),
        slotname=object@slotname,
        ...
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
    function (xstart, times, params, ...,
      zeronames = character(0),
      tcovar, covar,
      .getnativesymbolinfo = TRUE) {
      tryCatch(
        .Call(
          euler_model_simulator,
          func=efun,
          xstart=xstart,
          times=times,
          params=params,
          deltat=object@delta.t,
          method=0L,
          zeronames=zeronames,
          tcovar=tcovar,
          covar=covar,
          args=pairlist(...),
          gnsi=.getnativesymbolinfo
        ),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE)
        }
      )
    }
  }
)

setMethod(
  "plugin.handler",
  signature=signature(object='discreteRprocessPlugin'),
  definition=function (object, ...) {
    ep <- paste0("in ",sQuote("discrete.time.sim")," plugin: ")
    efun <- tryCatch(
      pomp.fun(
        f=object@step.fn,
        PACKAGE=object@PACKAGE,
        proto=quote(step.fun(x,t,params,...)),
        slotname=object@slotname,
        ...
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
    function (xstart, times, params, ...,
      zeronames = character(0),
      tcovar, covar,
      .getnativesymbolinfo = TRUE) {
      tryCatch(
        .Call(
          euler_model_simulator,
          func=efun,
          xstart=xstart,
          times=times,
          params=params,
          deltat=object@delta.t,
          method=2L,
          zeronames=zeronames,
          tcovar=tcovar,
          covar=covar,
          args=pairlist(...),
          gnsi=.getnativesymbolinfo
        ),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE)
        }
      )
    }
  }
)

setMethod(
  "plugin.handler",
  signature=signature(object='gillespieRprocessPlugin'),
  definition=function (object, ...) {
    ep <- paste0("in ",sQuote("gillespie.sim")," plugin: ")
    efun <- tryCatch(
      pomp.fun(
        f=object@rate.fn,
        PACKAGE=object@PACKAGE,
        proto=quote(rate.fun(j,x,t,params,...)),
        slotname=object@slotname,
        ...
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
    function (xstart, times, params,
      zeronames = character(0),
      tcovar, covar,
      .getnativesymbolinfo = TRUE,
      ...) {
      tryCatch(
        .Call(
          SSA_simulator,
          func=efun,
          xstart=xstart,
          times=times,
          params=params,
          e=numeric(0),
          vmatrix=object@v,
          dmatrix=object@d,
          deps=integer(0),
          tcovar=tcovar,
          covar=covar,
          zeronames=zeronames,
          hmax=object@hmax,
          args=pairlist(...),
          gnsi=.getnativesymbolinfo
        ),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE)
        }
      )
    }
  }
)


setMethod(
  "plugin.handler",
  signature=signature(object='onestepDprocessPlugin'),
  definition=function (object, ...) {
    ep <- paste0("in ",sQuote("onestep.dens")," plugin: ")
    efun <- tryCatch(
      pomp.fun(
        f=object@dens.fn,
        PACKAGE=object@PACKAGE,
        proto=quote(dens.fun(x1,x2,t1,t2,params,...)),
        slotname=object@slotname,
        ...
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
    function (x, times, params, ...,
      tcovar, covar, log = FALSE,
      .getnativesymbolinfo = TRUE) {
      tryCatch(
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
        ),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE)
        }
      )
    }
  }
)
