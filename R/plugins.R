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
    "kleapRprocessPlugin",
    contains="pompPlugin",
    slots=c(
        rate.fn="ANY",
        e="numeric",
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
    ep <- paste0("in ",sQuote("gillespie.sim")," plugin: ")
    if (missing(PACKAGE)) PACKAGE <- character(0)
    if (!is.matrix(v) || !is.matrix(d)) {
        stop(ep,sQuote("v")," and ",sQuote("d")," must be matrices.",
             call.=FALSE)
    }
    nvar <- nrow(v)
    nevent <- ncol(v)
    if ((nvar!=nrow(d))||(nevent!=ncol(d)))
        stop(ep,sQuote("v")," and ",sQuote("d")," must agree in dimension.",
             call.=FALSE)
    new("gillespieRprocessPlugin",
        rate.fn=rate.fun,v=v,d=d,
        slotname="rate.fn",
        csnippet=is(rate.fun,"Csnippet"),
        PACKAGE=PACKAGE)
}

kleap.sim <- function (rate.fun, e, v, d, PACKAGE) {
    ep <- paste0("in ",sQuote("kleap.sim")," plugin: ")
    if (missing(PACKAGE)) PACKAGE <- character(0)
    if (!is.matrix(v) || !is.matrix(d)) {
        stop(ep,sQuote("v")," and ",sQuote("d")," must be matrices.",
             call.=FALSE)
    }
    nvar <- nrow(v)
    nevent <- ncol(v)
    if ((nvar!=nrow(d))||(nevent!=ncol(d)))
        stop(ep,sQuote("v")," and ",sQuote("d")," must agree in dimension.",
             call.=FALSE)
    e <- as.numeric(e)
    if (nvar!=length(e))
        stop(ep,sQuote("e")," must have one entry for each state variable.",
             call.=FALSE)
    if (any((e>1)|(e<0)))
        stop(ep,"each element of ",sQuote("e")," must be in [0,1].",call.=FALSE)
    new("kleapRprocessPlugin",
        rate.fn=rate.fun,e=e,v=v,d=d,
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
                    mflag=0L, ## Gillespie's algorithm
                    xstart=xstart,
                    times=times,
                    params=params,
                    e=numeric(0),
                    vmatrix=object@v,
                    dmatrix=object@d,
                    tcovar=tcovar,
                    covar=covar,
                    zeronames=zeronames,
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
    signature=signature(object='kleapRprocessPlugin'),
    definition=function (object, ...) {
        ep <- paste0("in ",sQuote("kleap.sim")," plugin: ")
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
                    mflag=1L, ## K-leap algorithm
                    xstart=xstart,
                    times=times,
                    params=params,
                    e=object@e,
                    vmatrix=object@v,
                    dmatrix=object@d,
                    tcovar=tcovar,
                    covar=covar,
                    zeronames=zeronames,
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
