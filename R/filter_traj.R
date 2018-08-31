##' Filtering trajectories
##'
##' Trajectories drawn from the smoothing distribution
##'
##' The smoothing distribution is the distribution of
##' \deqn{X_t | Y_1=y^*_1, \dots, Y_T=y^*_T,}{Xt | Y1=y1*, \dots, YT=yT*,}
##' where \eqn{X_t}{Xt} is the latent state process, \eqn{Y_t}{Yt} is the observable process, \eqn{t} is time, and \eqn{T} is the time of the final observation.
##'
##' In a particle filter, the trajectories of the individual particles are not independent of one another, since they share ancestry.
##' However, a randomly sampled particle trajectory \eqn{X_1,\dots,X_T} is a draw from the smoothing distribution.
##' Seting \code{filter.traj = TRUE} in \code{\link{pfilter}} causes one such trajectory to be sampled.
##' By running multiple independent \code{pfilter} operations, one can thus build up a picture of the smoothing distribution.
##'
##' In particle MCMC (\code{\link{pmcmc}}), this operation is performed at each MCMC iteration.
##' Assuming the MCMC chain has converged, and after proper measures are taken to assure approximate independence of samples, \code{filter.traj} allows one to extract a sample from the smoothing distribution.
##'
##' @name filter.traj
##' @aliases filter.traj filter.traj,ANY-method filter.traj,missing-method
##' @include pfilter.R pmcmc.R
##' @rdname filter_traj
##' @family particle filter methods
##' @inheritParams filter.mean
NULL

setGeneric(
  "filter.traj",
  function (object,...) standardGeneric("filter.traj")
)

setMethod(
  "filter.traj",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("filter.traj","object")
  }
)

setMethod(
  "filter.traj",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("filter.traj",object)
  }
)

##' @name filter.traj-pfilterd_pomp
##' @aliases filter.traj,pfilterd_pomp-method
##' @rdname filter_traj
##'
##' @export
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
##' @export
setMethod(
  "filter.traj",
  signature=signature(object="pfilterList"),
  definition=function (object, vars, ...) {
    fts <- lapply(object,filter.traj,vars=vars,...)
    d <- sapply(fts,dim)
    if (!all(apply(d,1L,function(x)x==x[1L])))
      pStop("filter.traj","incommensurate dimensions.")
    d <- d[,1L]
    nm <- dimnames(fts[[1L]])
    x <- do.call(c,fts)
    dim(x) <- c(d,length(fts))
    dimnames(x) <- c(nm,list(chain=names(fts)))
    x
  }
)

##' @name filter.traj-pmcmcd_pomp
##' @aliases filter.traj,pmcmcd_pomp-method
##' @rdname filter_traj
##' @export
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
##' @export
setMethod(
  "filter.traj",
  signature=signature(object="pmcmcList"),
  definition=function (object, vars, ...) {
    fts <- lapply(object,filter.traj,vars=vars,...)
    d <- dim(fts[[1]])
    nm <- dimnames(fts[[1]])
    x <- do.call(c,fts)
    dim(x) <- c(d,length(fts))
    dimnames(x) <- c(nm,list(chain=names(fts)))
    x
  }
)
