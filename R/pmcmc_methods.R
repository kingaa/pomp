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

## extract the convergence record as a coda::mcmc object
setMethod(
  "traces",
  signature=signature(object="pmcmcd_pomp"),
  function (object, pars, ...) {
    if (missing(pars)) pars <- colnames(object@traces)
    coda::mcmc(object@traces[,pars,drop=FALSE])
  }
)

## extract the convergence records as a coda::mcmc.list object
setMethod(
  "traces",
  signature=signature(object="pmcmcList"),
  definition=function (object, ...) {
    coda::mcmc.list(lapply(object,traces,...))
  }
)

## extract the filtered trajectories from a pmcmcd_pomp
setMethod(
  "filter.traj",
  signature=signature(object="pmcmcd_pomp"),
  definition=function (object, ...) {
    filter.traj(as(object,"pfilterd_pomp"),...)
  }
)

## extract the filtered trajectories from a pmcmcList
setMethod(
  "filter.traj",
  signature=signature(object="pmcmcList"),
  definition=function (object, ...) {
    fts <- lapply(object,filter.traj,...)
    d <- dim(fts[[1]])
    nm <- dimnames(fts[[1]])
    x <- do.call(c,fts)
    dim(x) <- c(d,length(fts))
    dimnames(x) <- c(nm,list(chain=names(fts)))
    x
  }
)

setMethod(
  "plot",
  signature=signature(x="Pmcmc"),
  function (x, pars, ...) {
    plot(traces(x,pars),...)
  }
)

## extract the estimated log likelihood
setMethod(
  "logLik",
  signature=signature(object="pmcmcd_pomp"),
  definition=function (object, ...)
    object@loglik
)
