##' Filtering trajectories
##'
##' Trajectories drawn from the smoothing distribution
##'
##' @name Filtering trajectories
##' @aliases filter.traj filter.traj,ANY-method filter.traj,missing-method
##' @include pfilter.R pmcmc.R
##' @rdname filter_traj
NULL

setGeneric(
  "filter.traj",
  function (object,...) standardGeneric("filter.traj")
)

setMethod(
  "filter.traj",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("filter.traj"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "filter.traj",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("filter.traj")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' @name filter.traj-pfilterd_pomp
##' @aliases filter.traj,pfilterd_pomp-method
##' @rdname filter_traj
##'
##' @inheritParams filter.mean-kalmand_pomp
##'
setMethod(
  "filter.traj",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@filter.traj)
    object@filter.traj[vars,,,drop=FALSE]
  }
)

##' @name filter.traj-pfilterList
##' @aliases filter.traj,pfilterList-method
##' @rdname filter_traj
setMethod(
  "filter.traj",
  signature=signature(object="pfilterList"),
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

##' @name filter.traj-pmcmcd_pomp
##' @aliases filter.traj,pmcmcd_pomp-method
##' @rdname filter_traj
setMethod(
  "filter.traj",
  signature=signature(object="pmcmcd_pomp"),
  definition=function (object, vars, ...) {
    filter.traj(as(object,"pfilterd_pomp"),vars,...)
  }
)

##' @name filter.traj-pmcmcList
##' @aliases filter.traj,pmcmcList-method
##' @rdname filter_traj
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
