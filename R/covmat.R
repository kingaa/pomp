## helper function to extract a covariance matrix

setMethod(
  "covmat",
  signature=signature(object="pmcmcd_pomp"),
  definition=function (object, start = 1, thin = 1,
    expand = 2.38, ...) {
    covmat.internal(traces=as.matrix(traces(object,object@pars)),
      start=start,thin=thin,
      expand=expand)
  })

setMethod(
  "covmat",
  signature=signature(object="pmcmcList"),
  definition=function (object, start = 1, thin = 1,
    expand = 2.38, ...) {
    pars <- unique(c(sapply(object,slot,"pars")))
    covmat.internal(traces=as.array(traces(object,pars)),
      start=start,thin=thin,
      expand=expand)
  })

setMethod(
  "covmat",
  signature=signature(object="abcd_pomp"),
  definition=function (object, start = 1, thin = 1,
    expand = 2.38, ...) {
    covmat.internal(traces=as.matrix(traces(object,object@pars)),
      start=start,thin=thin,
      expand=expand)
  })

setMethod(
  "covmat",
  signature=signature(object="abcList"),
  definition=function (object, start = 1, thin = 1,
    expand = 2.38, ...) {
    pars <- unique(c(sapply(object,slot,"pars")))
    covmat.internal(traces=as.array(traces(object,pars)),
      start=start,thin=thin,
      expand=expand)
  })

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
    warning("in ",sQuote("covmat"),": only ",length(keep),
      " points being used to estimate covariance matrix",call.=FALSE)
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
