setGeneric(
    "filter.traj",
    function (object,...) standardGeneric("filter.traj")
)

setMethod(
  "filter.traj",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@filter.traj)
    object@filter.traj[vars,,,drop=FALSE]
  }
)

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

setMethod(
  "filter.traj",
  signature=signature(object="pmcmcd_pomp"),
  definition=function (object, ...) {
    filter.traj(as(object,"pfilterd_pomp"),...)
  }
)

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
  "filter.traj",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("filter.traj")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)
