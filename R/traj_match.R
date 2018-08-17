##' Trajectory matching
##'
##' Fitting a deterministic model to data.
##'
##' @rdname traj_match
##' @name Trajectory matching
##' @include trajectory.R pomp_class.R workhorses.R
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
##' Specifically, \code{traj.match} calls an optimizer (\code{\link{optim}},
##' \code{\link[subplex]{subplex}}, and \code{\link{sannbox}} are the currently
##' supported options) to minimize an objective function.  For any value of the
##' model parameters, this objective function is calculated by \enumerate{
##' \item computing the deterministic trajectory of the model given the
##' parameters.  This is the trajectory returned by \code{\link{trajectory}},
##' which relies on the model's deterministic skeleton as specified in the
##' construction of the \sQuote{pomp} object \code{object}.  \item evaluating
##' the negative log likelihood of the data under the measurement model given
##' the deterministic trajectory and the model parameters.  This is
##' accomplished via the model's \code{dmeasure} slot.  The negative log
##' likelihood is the objective function's value.  }
##'
##' The objective function itself --- in a form suitable for use with
##' \code{\link{optim}}-like optimizers --- is created by a call to
##' \code{traj.match.objfun}.  Specifically, \code{traj.match.objfun} will
##' return a function that takes a single numeric-vector argument that is
##' assumed to cotain the parameters named in \code{est}, in that order.
##'
##' @param object A \sQuote{\link{pomp}} object.  If \code{object} has no
##' \code{skeleton} slot, an error will be generated.
##' @param start named numeric vector containing an initial guess for
##' parameters.  By default \code{start=coef(object)} if the latter exists.
##' @param params optional named numeric vector of parameters.  This should
##' contain all parameters needed by the \code{skeleton} and \code{dmeasure}
##' slots of \code{object}.  In particular, any parameters that are to be
##' treated as fixed should be present here.  Parameter values given in
##' \code{params} for parameters named in \code{est} will be ignored.  By
##' default, \code{params=coef(object)} if the latter exists.
##' @param est character vector containing the names of parameters to be
##' estimated.  In the case of \code{traj.match.objfun}, the objective function
##' that is constructed will assume that its argument contains the parameters
##' in this order.
##' @param method Optimization method.  Choices are
##' \code{\link[subplex]{subplex}}, \dQuote{sannbox}, any of the methods used
##' by \code{\link{optim}}, and \code{nloptr}.  The latter makes available all
##' the optimization algorithms of the \pkg{nloptr} package
##' (\url{https://cran.r-project.org/package=nloptr}).
##' @param transform logical; if \code{TRUE}, optimization is performed on the
##' transformed scale.
##' @param \dots Extra arguments that will be passed either to the optimizer
##' (\code{\link{optim}}, \code{\link[subplex]{subplex}},
##' \code{\link[nloptr]{nloptr}}, or \code{\link{sannbox}}, via their
##' \code{control} (\code{optim}, \code{subplex}, \code{sannbox}) or
##' \code{opts} (\code{nloptr}) lists) or to the ODE integrator.  In
##' \code{traj.match}, extra arguments will be passed to the optimizer.  In
##' \code{traj.match.objfun}, extra arguments are passed to
##' \code{\link{trajectory}}.  If extra arguments are needed by both optimizer
##' and \code{\link{trajectory}}, construct an objective function first using
##' \code{traj.match.objfun}, then give this objective function to the
##' optimizer.
##'
##' @return
##' \code{traj.match} returns an object of class
##' \sQuote{traj_matched_pomp}.  Available methods for objects of this class
##' include \code{summary} and \code{logLik}.
##'
##' \code{traj.match.objfun} returns a function suitable for use as an
##' objective function in an \code{\link{optim}}-like optimizer.
##'
##' @seealso \code{\link{trajectory}}, \code{\link{optim}},
##' \code{\link[subplex]{subplex}}
NULL

setClass(
  "traj_matched_pomp",
  contains="pomp",
  slots=c(
    transform="logical",
    est="character",
    evals="integer",
    convergence="integer",
    msg="character",
    value="numeric"
  )
)

setGeneric("traj.match.objfun",
  function(object,...)standardGeneric("traj.match.objfun"))

setGeneric("traj.match",
  function(object,...)standardGeneric("traj.match"))

##' @rdname traj_match
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

##' @rdname traj_match
setMethod(
  "traj.match",
  signature=signature(object="pomp"),
  function (object, start, est = character(0),
    method = c("Nelder-Mead","subplex","SANN","BFGS",
      "sannbox","nloptr"),
    transform = FALSE, ...)
  {

    if (missing(start)) start <- coef(object)
    if (is.list(start)) start <- unlist(start)

    method <- match.arg(method)
    est <- as.character(est)
    transform <- as.logical(transform)

    m <- minim.internal(
      objfun=traj.match.objfun(
        object=object,
        params=start,
        est=est,
        transform=transform
      ),
      start=start,
      est=est,
      object=object,
      method=method,
      transform=transform,
      ...
    )

    ## fill params slot appropriately
    coef(object) <- m$params

    ## fill states slot appropriately
    x <- trajectory(object)
    object@states <- array(data=x,dim=dim(x)[c(1L,3L)])
    rownames(object@states) <- rownames(x)

    new(
      "traj_matched_pomp",
      object,
      transform=transform,
      est=est,
      value=-m$value,
      evals=m$evals,
      convergence=m$convergence,
      msg=m$msg
    )
  }
)

##' @rdname traj_match
setMethod(
  "traj.match",
  signature=signature(object="traj_matched_pomp"),
  function (object, est, transform, ...)
  {
    if (missing(est)) est <- object@est
    if (missing(transform)) transform <- object@transform

    traj.match(as(object,"pomp"),est=est,transform=transform,...)
  }
)

setMethod(
  "traj.match",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("traj.match"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "traj.match",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("traj.match")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
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

tmof.internal <- function (object, params, est, transform, ...) {

  ep <- paste0("in ",sQuote("traj.match.objfun"),": ")
  object <- as(object,"pomp")
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
