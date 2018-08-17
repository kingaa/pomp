##' Traces
##'
##' Retrieve the history of an iterative calculation.
##'
##' @name Traces
##' @rdname traces
##' @aliases traces traces,missing-method traces,ANY-method
##' @include pmcmc.R mif2.R abc.R
NULL

setGeneric(
  "traces",
  function (object, ...)
    standardGeneric("traces")
)

setMethod(
  "traces",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("traces"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "traces",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("traces")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' @name traces-mif2d_pomp
##' @aliases traces,mif2d_pomp-method
##' @rdname traces
##'
##' @param object an object of class extending \sQuote{pomp}, the result of the application of a parameter estimation algorithm
##' @param pars names of parameters
##' @param transform logical; should the traces be transformed back onto the natural scale?
##' @param \dots ignored or (in the case of the listies, passed to the more primitive function)
##'
##' @return
##' When \code{object} is the result of a \code{\link{mif2}} calculation,
##' \code{traces(object, pars, transform = FALSE)} returns the traces of the parameters named in \code{pars}.
##' By default, the traces of all parameters are returned.
##' Note that, if the computation was performed with transformed parameters, the traces are on the estimation scale.
##' If \code{transform=TRUE}, the parameters are transformed from the estimation scale onto the natural scale.
##'
setMethod(
  "traces",
  signature=signature(object="mif2d_pomp"),
  definition=function (object, pars, transform = FALSE, ...) {
    traces.internal(object=object,pars=pars,transform=transform,...)
  }
)

##' @name traces-mif2List
##' @aliases traces,mif2List-method traces-Mif2 traces,Mif2-method
##' @rdname traces
setMethod(
  "traces",
  signature=signature(object="mif2List"),
  definition=function (object, ...) {
    lapply(object,traces,...)
  }
)

##' @name traces-abcd_pomp
##' @aliases traces,abcd_pomp-method
##' @rdname traces
##' @return
##' When \code{object} is a \sQuote{abcd_pomp}, \code{traces(object)}
##' extracts the traces as a \code{coda::mcmc}.
setMethod(
  "traces",
  signature=signature(object="abcd_pomp"),
  definition=function (object, pars, ...) {
    if (missing(pars)) pars <- colnames(object@traces)
    coda::mcmc(object@traces[,pars,drop=FALSE])
  }
)

##' @name traces-abcList
##' @aliases traces,abcList-method traces-Abc traces,Abc-method
##' @rdname traces
##' @return
##' When \code{object} is a \sQuote{abcList}, \code{traces(object)}
##' extracts the traces as a \code{coda::mcmc.list}.
setMethod(
  "traces",
  signature=signature(object="abcList"),
  definition=function (object, ...) {
    coda::mcmc.list(lapply(object,traces,...))
  }
)

##' @name traces-pmcmcd_pomp
##' @aliases traces,pmcmcd_pomp-method
##' @rdname traces
##' @return
##' When \code{object} is a \sQuote{pmcmcd_pomp}, \code{traces(object)}
##' extracts the traces as a \code{coda::mcmc}.
##' @details
##' Note that \code{\link{pmcmc}} does not currently support parameter transformations.
setMethod(
  "traces",
  signature=signature(object="pmcmcd_pomp"),
  function (object, pars, ...) {
    if (missing(pars)) pars <- colnames(object@traces)
    coda::mcmc(object@traces[,pars,drop=FALSE])
  }
)

##' @name traces-pmcmcList
##' @aliases traces,pmcmcList-method traces-Pmcmc traces,Pmcmc-method
##' @rdname traces
##' @return
##' When \code{object} is a \sQuote{pmcmcList}, \code{traces(object)}
##' extracts the traces as a \code{coda::mcmc.list}.
setMethod(
  "traces",
  signature=signature(object="pmcmcList"),
  definition=function (object, ...) {
    coda::mcmc.list(lapply(object,traces,...))
  }
)

traces.internal <- function (object, pars, transform = FALSE, ...) {
  transform <- as.logical(transform)
  if (transform) {
    retval <- cbind(
      object@traces[,c(1,2)],
      t(
        partrans(
          object,
          params=t(object@traces)[-c(1,2),,drop=FALSE],
          dir="fromEst"
        )
      )
    )
    names(dimnames(retval)) <- names(dimnames(object@traces))
  } else {
    retval <- object@traces
  }
  if (missing(pars)) {
    retval
  } else {
    pars <- as.character(pars)
    bad.pars <- setdiff(pars,colnames(retval))
    if (length(bad.pars)>0)
      stop(
        "in ",sQuote("traces"),": name(s) ",
        paste(sQuote(bad.pars),collapse=","),
        " correspond to no parameter(s) in ",
        if (transform) sQuote("traces(object,transform=TRUE)")
        else sQuote("traces(object,transform=FALSE)"),
        call.=FALSE
      )
    retval[,pars,drop=FALSE]
  }
}
