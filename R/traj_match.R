##' Trajectory matching
##'
##' Fitting a deterministic model to data.
##'
##' In trajectory matching, one attempts to minimize the discrepancy between a \acronym{POMP} model's predictions and data under the assumption that the process model is deterministic and all discrepancies between model and data are due to measurement error.
##' The measurement model likelihood (evaluated by \code{dmeasure}), or rather its negative, is a natural measure of the discrepancy.
##'
##' Thus trajectory matching is a generalization of the traditional nonlinear least squares approach.
##' In particular, if, on some scale, the measurement errors are normal with constant variance, then trajectory matching is equivalent to least squares (on that particular scale).
##'
##' This function attempts to match trajectories of a model's deterministic
##' skeleton to data.  Trajectory matching is equivalent to maximum likelihood
##' estimation under the assumption that process noise is entirely absent,
##' i.e., that all stochasticity is measurement error.  Accordingly, this
##' method uses only the \code{skeleton} and \code{dmeasure} components of a
##' \acronym{POMP} model.
##'
##' In \pkg{pomp}, trajectory matching is the term used for maximizing the
##' likelihood of the data under the assumption that there is no process noise.
##' For any value of the model parameters, this objective function is calculated by
##' \enumerate{
##' \item computing the deterministic trajectory of the model given the
##' parameters.  This is the trajectory returned by \code{\link{trajectory}},
##' which relies on the model's deterministic skeleton as specified in the
##' construction of the \sQuote{pomp} object \code{object}.
##' \item evaluating
##' the negative log likelihood of the data under the measurement model given
##' the deterministic trajectory and the model parameters.  This is
##' accomplished via the model's \code{dmeasure} slot.  The negative log
##' likelihood is the objective function's value.
##' }
##'
##' The objective function itself --- in a form suitable for use with
##' \code{\link{optim}}-like optimizers --- is created by a call to
##' \code{traj.match.objfun}.  Specifically, \code{traj.match.objfun} will
##' return a function that takes a single numeric-vector argument that is
##' assumed to cotain the parameters named in \code{est}, in that order.
##'
##' @name traj.match
##' @rdname traj_match
##' @include trajectory.R pomp_class.R workhorses.R
##' @aliases traj.match.objfun,missing-method traj.match.objfun,ANY-method
##'
##' @return
##' \code{traj.match.objfun} returns a function suitable for use as an
##' objective function in an \code{\link{optim}}-like optimizer.
##'
##' @seealso \code{\link{trajectory}}, \code{\link{optim}},
##' \code{\link[subplex]{subplex}}
NULL

setGeneric(
  "traj.match.objfun",
  function (object, ...)
    standardGeneric("traj.match.objfun")
)

setMethod(
  "traj.match.objfun",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("traj.match.objfun"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "traj.match.objfun",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("traj.match.objfun")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

##' @name traj.match.objfun-pomp
##' @aliases traj.match.objfun traj.match.objfun,pomp-method
##' @rdname traj_match
##'
##' @inheritParams probe.match
##'
setMethod(
  "traj.match.objfun",
  signature=signature(object="pomp"),
  function (object, params, est, transform = FALSE, ...) {

    tmof.internal(
      object=object,
      params=params,
      est=est,
      transform=transform,
      ...
    )

  }
)

tmof.internal <- function (object, params, est, transform, ...) {

  ep <- paste0("in ",sQuote("traj.match.objfun"),": ")

  object <- pomp(object,...)

  if (missing(est)) est <- character(0)
  est <- as.character(est)
  transform <- as.logical(transform)

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)
  if ((!is.numeric(params))||(is.null(names(params))))
    stop(ep,sQuote("params")," must be a named numeric vector",call.=FALSE)
  if (transform)
    params <- partrans(object,params,dir="toEst")
  par.est.idx <- match(est,names(params))
  if (any(is.na(par.est.idx)))
    stop(ep,"parameter(s): ",
      paste(sapply(est[is.na(par.est.idx)],sQuote),collapse=","),
      " not found in ",sQuote("params"),call.=FALSE)

  function (par) {
    pompLoad(object)
    on.exit(pompUnload(object))
    params[par.est.idx] <- par
    if (transform)
      tparams <- partrans(object,params,dir="fromEst")
    d <- dmeasure(
      object,
      y=object@data,
      x=trajectory(
        object,
        params=if (transform) tparams else params,
        ...
      ),
      times=time(object),
      params=if (transform) tparams else params,
      log=TRUE
    )
    -sum(d)
  }
}
