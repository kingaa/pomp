##' Estimate a covariance matrix from algorithm traces
##'
##' A helper function to extract a covariance matrix.
##'
##' @name covmat
##' @aliases covmat,missing-method covmat,ANY-method
##' @include pmcmc.R abc.R probe.R
##' @rdname covmat
##' @family extraction methods
##'
##' @seealso \link[=proposals]{MCMC proposals}.
NULL

setGeneric(
  "covmat",
  function (object, ...)
    standardGeneric("covmat")
)

##' @rdname covmat
##' @param object an object extending \sQuote{pomp}
##' @param start the first iteration number to be used in estimating the covariance matrix.
##' Setting \code{thin > 1} allows for a burn-in period.
##' @param thin factor by which the chains are to be thinned
##' @param expand the expansion factor
##' @param \dots ignored
##' @return
##' When \code{object} is the result of a \code{pmcmc} or \code{abc} computation,
##' \code{covmat(object)} gives the covariance matrix of the chains.
##' This can be useful, for example, in tuning the proposal distribution.
##' @export
setMethod(
  "covmat",
  signature=signature(object="pmcmcd_pomp"),
  definition=function (object, start = 1, thin = 1,
    expand = 2.38, ...) {
    covmat.internal(traces=as.matrix(traces(object,object@pars)),
      start=start,thin=thin,expand=expand)
  })

##' @rdname covmat
##' @export
setMethod(
  "covmat",
  signature=signature(object="pmcmcList"),
  definition=function (object, start = 1, thin = 1,
    expand = 2.38, ...) {
    pars <- unique(c(sapply(object,slot,"pars")))
    covmat.internal(traces=as.array(traces(object,pars)),
      start=start,thin=thin,expand=expand)
  })

##' @rdname covmat
##' @export
setMethod(
  "covmat",
  signature=signature(object="abcd_pomp"),
  definition=function (object, start = 1, thin = 1,
    expand = 2.38, ...) {
    covmat.internal(traces=as.matrix(traces(object,object@pars)),
      start=start,thin=thin,expand=expand)
  })

##' @rdname covmat
##' @export
setMethod(
  "covmat",
  signature=signature(object="abcList"),
  definition=function (object, start = 1, thin = 1,
    expand = 2.38, ...) {
    pars <- unique(c(sapply(object,slot,"pars")))
    covmat.internal(traces=as.array(traces(object,pars)),
      start=start,thin=thin,expand=expand)
  })

##' @rdname covmat
##' @importFrom stats var
##' @return
##' When \code{object} is a \sQuote{probed_pomp} object (i.e., the result
##' of a \code{probe} computation), \code{covmat(object)} returns the
##' covariance matrix of the probes, as applied to simulated data.
##' @export
setMethod(
  "covmat",
  signature=signature(object="probed_pomp"),
  definition=function (object, ...) {
    n <- nrow(object@simvals)
    crossprod(object@simvals)/(n-1)
  })

covmat.internal <- function (traces, start, thin, expand = 2.38, ...) {
  dd <- dim(traces)
  nms <- colnames(traces)
  keep <- seq.int(from=as.integer(start),to=dd[1],by=as.integer(thin))
  if (length(keep) < 100)
    pWarn("covmat","only ",length(keep),
      ngettext(length(keep)," point is"," points are"),
      " being used to estimate covariance matrix.")
  if (length(dd)==2L) {
    traces <- traces[keep,,drop=FALSE]
  } else if (length(dd)==3L) {
    traces <- aperm(traces[keep,,,drop=FALSE],c(1L,3L,2L))
    dd <- dim(traces)
    dim(traces) <- c(dd[1L]*dd[2L],dd[3L])
  }
  v <- var(traces)
  keep <- which(diag(v)>0)
  v <- expand^2*v[keep,keep]/length(keep)
  dimnames(v) <- list(nms[keep],nms[keep])
  v
}
