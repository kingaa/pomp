## history of an iterative calculation

setGeneric(
    "traces",
    function (object, ...)
        standardGeneric("traces")
)


setMethod(
  "traces",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("traces")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "traces",
  signature=signature(object="mif2d_pomp"),
  definition=function (object, pars, transform = FALSE, ...) {
    traces.internal(object=object,pars=pars,transform=transform,...)
  }
)

setMethod(
  "traces",
  signature=signature(object="mif2List"),
  definition=function (object, ...) {
    lapply(object,traces,...)
  }
)

## extract the traces as a coda::mcmc object
setMethod(
  "traces",
  signature=signature(object="abcd_pomp"),
  definition=function (object, pars, ...) {
    if (missing(pars)) pars <- colnames(object@traces)
    coda::mcmc(object@traces[,pars,drop=FALSE])
  }
)

## extract the traces as a coda::mcmc.list object
setMethod(
  "traces",
  signature=signature(object="abcList"),
  definition=function (object, ...) {
    coda::mcmc.list(lapply(object,traces,...))
  }
)

## extract the traces as a coda::mcmc object
setMethod(
  "traces",
  signature=signature(object="pmcmcd_pomp"),
  function (object, pars, ...) {
    if (missing(pars)) pars <- colnames(object@traces)
    coda::mcmc(object@traces[,pars,drop=FALSE])
  }
)

## extract the traces as a coda::mcmc.list object
setMethod(
  "traces",
  signature=signature(object="pmcmcList"),
  definition=function (object, ...) {
    coda::mcmc.list(lapply(object,traces,...))
  }
)

## extract the traces
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
