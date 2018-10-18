##' Traces
##'
##' Retrieve the history of an iterative calculation.
##'
##' @name traces
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
    reqd_arg("traces","object")
  }
)

setMethod(
  "traces",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("traces",object)
  }
)

##' @name traces-mif2d_pomp
##' @aliases traces,mif2d_pomp-method
##' @rdname traces
##'
##' @importFrom coda mcmc mcmc.list
##'
##' @param object an object of class extending \sQuote{pomp}, the result of the application of a parameter estimation algorithm
##' @param pars names of parameters
##' @param transform logical; should the traces be transformed back onto the natural scale?
##' @param \dots ignored or (in the case of the listie, passed to the more primitive function)
##'
##' @return
##' When \code{object} is the result of a \code{\link{mif2}} calculation,
##' \code{traces(object, pars, transform = FALSE)} returns the traces of the parameters named in \code{pars}.
##' By default, the traces of all parameters are returned.
##' Note that, if the computation was performed with transformed parameters, the traces are on the estimation scale.
##' If \code{transform=TRUE}, the parameters are transformed from the estimation scale onto the natural scale.
##'
##' @export
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
##' @export
setMethod(
  "traces",
  signature=signature(object="mif2List"),
  definition=function (object, pars, ...) {
    lapply(object,traces,pars,...)
  }
)

##' @name traces-abcd_pomp
##' @aliases traces,abcd_pomp-method
##' @rdname traces
##' @return
##' When \code{object} is a \sQuote{abcd_pomp}, \code{traces(object)}
##' extracts the traces as a \code{coda::mcmc}.
##' @export
setMethod(
  "traces",
  signature=signature(object="abcd_pomp"),
  definition=function (object, pars, ...) {
    tr <- traces.internal(object,pars=pars)
    coda::mcmc(tr)
  }
)

##' @name traces-abcList
##' @aliases traces,abcList-method traces-Abc traces,Abc-method
##' @rdname traces
##' @return
##' When \code{object} is a \sQuote{abcList}, \code{traces(object)}
##' extracts the traces as a \code{coda::mcmc.list}.
##' @export
setMethod(
  "traces",
  signature=signature(object="abcList"),
  definition=function (object, pars, ...) {
    coda::mcmc.list(lapply(object,traces,pars=pars))
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
##' @export
setMethod(
  "traces",
  signature=signature(object="pmcmcd_pomp"),
  function (object, pars, ...) {
    tr <- traces.internal(object,pars=pars)
    coda::mcmc(tr)
  }
)

##' @name traces-pmcmcList
##' @aliases traces,pmcmcList-method traces-Pmcmc traces,Pmcmc-method
##' @rdname traces
##' @return
##' When \code{object} is a \sQuote{pmcmcList}, \code{traces(object)}
##' extracts the traces as a \code{coda::mcmc.list}.
##' @export
setMethod(
  "traces",
  signature=signature(object="pmcmcList"),
  definition=function (object, pars, ...) {
    coda::mcmc.list(lapply(object,traces,pars=pars))
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
      pStop("traces",
        ngettext(length(bad.pars),"parameter ","parameters "),
        paste(sQuote(bad.pars),collapse=",")," not found.")

    retval[,pars,drop=FALSE]

  }
}
