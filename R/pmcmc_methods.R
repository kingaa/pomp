## this file defines methods for the 'pmcmc' and 'pmcmcList' classes

## pmcmcList class
setClass(
  "pmcmcList",
  contains='list',
  validity=function (object) {
    if (!all(sapply(object,is,"pmcmc"))) {
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
    TRUE
  }
)

setClassUnion("Pmcmc",c("pmcmc","pmcmcList"))

setMethod(
  "c",
  signature=signature(x="Pmcmc"),
  definition=function (x, ...) {
    y <- lapply(
      list(x,...),
      function (z) {
        if (is(z,"list"))
          as(z,"list")
        else
          list(z)
      }
    )
    new("pmcmcList",do.call(c,y))
  }
)

setMethod(
  "[",
  signature=signature(x="pmcmcList"),
  definition=function(x, i, ...)
    new("pmcmcList",as(x,"list")[i])
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
  "print",
  signature=signature(x="pmcmc"),
  definition=function (x, ...) {
    cat("<object of class ",sQuote("pmcmc"),">\n",sep="")
    invisible(x)
  }
)

setMethod(
  "print",
  signature=signature(x="pmcmcList"),
  definition=function (x, ...) {
    lapply(as(x,"list"),print)
    invisible(x)
  }
)


## extract the estimated log likelihood
setMethod('logLik','pmcmc',function(object,...)object@loglik)
setMethod('logLik','pmcmcList',function(object,...)sapply(object,slot,"loglik"))
