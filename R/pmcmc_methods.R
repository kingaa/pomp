## this file defines methods for the 'pmcmc' and 'pmcmcList' classes

## pmcmcList class
setClass(
  "pmcmcList",
  contains="list",
  validity=function (object) {
    if (length(object) > 0) {
      if (!all(vapply(object,is,logical(1),"pmcmc"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
      d <- sapply(object,function(x)dim(x@conv.rec))
      if (!all(apply(d,1,diff)==0)) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": to be combined, ",sQuote("pmcmc"),
          " objects must have chains of equal length"
        )
        return(retval)
      }
    }
    TRUE
  }
)

setClassUnion("Pmcmc",c("pmcmc","pmcmcList"))

setMethod(
  "cc",
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

c.Pmcmc <- cc

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
  'conv.rec',
  signature=signature(object='pmcmc'),
  function (object, pars, ...) {
    if (missing(pars)) pars <- colnames(object@conv.rec)
    coda::mcmc(object@conv.rec[,pars,drop=FALSE])
  }
)

## extract the convergence records as a coda::mcmc.list object
setMethod(
  'conv.rec',
  signature=signature(object='pmcmcList'),
  definition=function (object, ...) {
    coda::mcmc.list(lapply(object,conv.rec,...))
  }
)

## extract the filtered trajectories from a pmcmc
setMethod(
  'filter.traj',
  signature=signature(object='pmcmc'),
  definition=function (object, ...) {
    filter.traj(as(object,"pfilterd.pomp"),...)
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

## plot pmcmc object
setMethod(
  "plot",
  signature=signature(x='Pmcmc'),
  function (x, pars, ...) {
    plot(conv.rec(x,pars),...)
  }
)

setMethod(
  "show",
  signature=signature(object="pmcmc"),
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
setMethod('logLik','pmcmc',function(object,...)object@loglik)
setMethod('logLik','pmcmcList',function(object,...)sapply(object,slot,"loglik"))
