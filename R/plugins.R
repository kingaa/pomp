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
    warning("argument", sQuote("d"), "is deprecated; updates to the simulation",
            " algorithm have made it unnecessary", call. = FALSE)
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

gillespie.hl.sim <- function(..., .pre = "", .post = "", hmax = Inf){
    ep <- paste0("in ",sQuote("gillespie.hl.sim")," plugin: ")
    PACKAGE <- character(0)
    args <- list(...)

    codeChunks <- lapply(args, "[[", 1)
    checkCode <- function(x) {
      inh <- inherits(x, what = c("Csnippet", "character"), which = TRUE)
      if(!any(inh)) {
        stop(ep,"the first list element of each event argument should be a",
             " Csnippet or string", call.=FALSE)
      }
      if (length(x) != 1){
        stop(ep,"the length of the first list element of each event",
             " argument should be equal to 1", call.=FALSE)
      }
      if (inh[1] == 1) {
        x@text
      } else {
        x
      }
    }
    codeChunks <- lapply(codeChunks, checkCode)
    .pre <- checkCode(.pre)
    .post <- checkCode(.post)

    stoich <- lapply(args, "[[", 2)
    checkStoic <- function(x){
        if (! typeof(x) %in%  c("integer", "double")){
            stop(ep,"the second list element of each event argument should be",
                 " a numeric or integer vector", call.=FALSE)
        }
        if (is.null(names(x))) {
            stop(ep,"the second list element of each event argument should be",
                 " a named vector", call.=FALSE)
        }
    }
    lapply(stoich, checkStoic)

    ## Create C snippet of switch statement
    inds <- seq_along(args)
    makeCase <- function(i, chunk){
      label <- paste0("case ", i, ":{")
      endcase <- "break;\n}"
      paste(label, chunk, endcase, sep="\n")
    }
    body <- Map(makeCase, inds, codeChunks)
    body <- paste(body, collapse="\n")
    lastCase <- "default:\nerror(\"unrecognized event %d\",j);\nbreak;\n"
    header <- paste0(.pre, "\nswitch (j) {")
    footer <- paste0(lastCase, "}\n", .post)
    rate.fn <- Csnippet(paste(header, body, footer, sep="\n"))

    ### Create v matrix
    ## By coercing the vectors to a data frame and then using rbind,
    ## we can ensure that all stoichiometric coefficients for the same
    ## state variables are in the same column even if the vectors in
    ## stoich have differently ordered names. Also, rbind will fail if
    ## the set of variables in each data frame is not the same.
    stoichdf <- lapply(stoich, function(x) data.frame(as.list(x)))
    v <- t(data.matrix(do.call(rbind, stoichdf)))
    new("gillespieRprocessPlugin",
      rate.fn=rate.fn, v=v, hmax=hmax,
      slotname="rate.fn",
      csnippet=TRUE,
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
    args <- list(...)
    function (xstart, times, params,
      zeronames = character(0),
      tcovar, covar,
      .getnativesymbolinfo = TRUE,
              ...) {
      if (anyDuplicated(rownames(object@v))){
        stop(ep,"duplicates in rownames of ",sQuote("v"), call.=FALSE)
      }
      tryCatch(
        .Call(
          SSA_simulator,
          func=efun,
          xstart=xstart,
          times=times,
          params=params,
          vmatrix=object@v,
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
