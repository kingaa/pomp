##' Filtering trajectories
##'
##' Drawing from the smoothing distribution
##'
##' The smoothing distribution is the distribution of
##' \deqn{X(t_k) | Y(t_1)=y^*_1, \dots, Y(t_n)=y^*_n,}{Xk | Y1=y1*, \dots, Yn=yn*,}
##' where \eqn{X(t_k)}{Xk} is the latent state process and \eqn{Y(t_k)}{Yk} is the observable process at time \eqn{t_k}{tk}, and \eqn{n} is the number of observations.
##'
##' To draw samples from this distribution, one can run a number of independent particle filter (\code{\link{pfilter}}) operations, sampling the full trajectory of \emph{one} randomly-drawn particle from each one.
##' One should view these as \emph{weighted} samples from the smoothing distribution, where the weights are the \emph{likelihoods} returned by each of the \code{\link{pfilter}} computations.
##'
##' One accomplishes this by setting \code{filter.traj = TRUE} in each \code{\link{pfilter}} computation and extracting the trajectory using the \code{filter.traj} command.
##'
##' In particle MCMC (\code{\link{pmcmc}}), the tracking of an individual trajectory is performed automatically.
##' 
##' @name filter.traj
##' @aliases filter.traj filter.traj,ANY-method filter.traj,missing-method
##' @include pfilter.R pmcmc.R
##' @rdname filter_traj
##' @family particle_filter_methods
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
