## this file defines methods for the 'pmcmcd_pomp' and 'pmcmcList' classes

## pmcmcList class
setClass(
  "pmcmcList",
  contains="list",
  validity=function (object) {
    if (length(object) > 0) {
      if (!all(vapply(object,is,logical(1),"pmcmcd_pomp"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
      d <- sapply(object,function(x)dim(x@traces))
      if (!all(apply(d,1,diff)==0)) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": to be combined, ",sQuote("pmcmcd_pomp"),
          " objects must have chains of equal length"
        )
        return(retval)
      }
    }
    TRUE
  }
)

setClassUnion("Pmcmc",c("pmcmcd_pomp","pmcmcList"))

setMethod(
  "concat",
  signature=signature(...="Pmcmc"),
  definition=function (...) {
    y <- lapply(
      list(...),
      function (z) {
        if (is(z,"list"))
          setNames(as(z,"list"),names(z))
        else
          z
      }
    )
    new("pmcmcList",unlist(y))
  }
)

c.Pmcmc <- concat

setMethod(
  "[",
  signature=signature(x="pmcmcList"),
  definition=function(x, i, ...) {
    y <- as(x,"list")
    names(y) <- names(x)
    y <- unlist(y[i])
    if (is.null(y)) {
      list(NULL)
    } else {
      new("pmcmcList",y)
    }
  }
)

## extract the convergence record as a coda::mcmc object
setMethod(
  'traces',
  signature=signature(object='pmcmcd_pomp'),
  function (object, pars, ...) {
    if (missing(pars)) pars <- colnames(object@traces)
    coda::mcmc(object@traces[,pars,drop=FALSE])
  }
)

## extract the convergence records as a coda::mcmc.list object
setMethod(
  'traces',
  signature=signature(object='pmcmcList'),
  definition=function (object, ...) {
    coda::mcmc.list(lapply(object,traces,...))
  }
)

## extract the filtered trajectories from a pmcmcd_pomp
setMethod(
  'filter.traj',
  signature=signature(object='pmcmcd_pomp'),
  definition=function (object, ...) {
    filter.traj(as(object,"pfilterd_pomp"),...)
  }
)

## extract the filtered trajectories from a pmcmcList
setMethod(
  'filter.traj',
  signature=signature(object='pmcmcList'),
  definition=function (object, ...) {
    lapply(object,filter.traj,...)
  }
)

setMethod(
  "plot",
  signature=signature(x='Pmcmc'),
  function (x, pars, ...) {
    plot(traces(x,pars),...)
  }
)

setMethod(
  "show",
  signature=signature(object="pmcmcd_pomp"),
  definition=function (object) {
    cat("<object of class ",sQuote("abc"),">\n",sep="")
    invisible(NULL)
  }
)

setMethod(
  "show",
  signature=signature(object="pmcmcList"),
  definition=function (object) {
    y <- as(object,"list")
    names(y) <- names(object)
    show(y)
  }
)

## extract the estimated log likelihood
setMethod('logLik','pmcmcd_pomp',function(object,...)object@loglik)
setMethod('logLik','pmcmcList',function(object,...)sapply(object,slot,"loglik"))
