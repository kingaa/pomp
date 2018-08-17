##' Trajectory matching
##'
##' Estimation of parameters for deterministic \acronym{POMP} models
##'
##' In trajectory matching, one attempts to minimize the discrepancy between a \acronym{POMP} model's predictions and data under the assumption that the latent state process is deterministic and all discrepancies between model and data are due to measurement error.
##' The measurement model likelihood (\code{dmeasure}), or rather its negative, is the natural measure of the discrepancy.
##'
##' Trajectory matching is a generalization of the traditional nonlinear least squares approach.
##' In particular, if, on some scale, measurement errors are normal with constant variance, then trajectory matching is equivalent to least squares on that particular scale.
##'
##' \code{traj.match.objfun} constructs an objective function that evaluates the likelihood function.
##' It can be passed to any one of a variety of numerical optimization routines, which will adjust model parameters to minimize the discrepancies between the power spectrum of model simulations and that of the data.
##'
##' @name traj.match
##' @docType methods
##' @rdname traj_match
##' @include trajectory.R pomp_class.R workhorses.R
##' @aliases traj.match.objfun,missing-method traj.match.objfun,ANY-method
##'
##' @return
##' \code{traj.match.objfun} constructs a stateful objective function for spectrum matching.
##' Specifically, \code{traj.match.objfun} returns an object of class \sQuote{traj_match_objfun}, which is a function suitable for use in an \code{\link{optim}}-like optimizer.
##' In particular, this function takes a single numeric-vector argument that is assumed to contain the parameters named in \code{est}, in that order.
##' When called, it will return the negative log likelihood.
##' It is a stateful function:
##' Each time it is called, it will remember the values of the parameters and its estimate of the log likelihood.
##'
##' @section Important Note:
##' Since \pkg{pomp} cannot guarantee that the \emph{final} call an optimizer makes to the function is a call \emph{at} the optimum, it cannot guarantee that the parameters stored in the function are the optimal ones.
##' One should check that the parameters agree with those that are returned by the optimizer.
##'
##' @seealso \code{\link{trajectory}}, \code{\link{optim}},
##' \code{\link[subplex]{subplex}}, \code{\link[nloptr]{nloptr}}
NULL

setClass(
  "traj_match_objfun",
  contains="function",
  slots=c(
    env="environment"
  )
)

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
##' @inheritParams trajectory
##'
##' @param ode_control optional list;
##' the elements of this list will be passed to \code{\link[=deSolve]{ode}}.
##' @param \dots additional arguments will modify the model structure
##'
setMethod(
  "traj.match.objfun",
  signature=signature(object="pomp"),
  function (object, params, est, transform = FALSE, fail.value = NA,
    ode_control = list(), ...) {

    tmof.internal(
      object=object,
      params=params,
      est=est,
      transform=transform,
      fail.value=fail.value,
      ode_control=ode_control,
      ...
    )

  }
)

tmof.internal <- function (object, params, est, transform, fail.value,
  ode_control, ...) {

  ep <- paste0("in ",sQuote("traj.match.objfun"),": ")

  object <- pomp(object,params=params,...)

  transform <- as.logical(transform)
  fail.value <- as.numeric(fail.value)

  if (missing(est)) est <- character(0)
  est <- as.character(est)
  est <- est[nzchar(est)]

  params <- coef(object,transform=transform)

  idx <- match(est,names(params))
  if (any(is.na(idx))) {
    missing <- est[is.na(idx)]
    stop(ep,ngettext(length(missing),"parameter","parameters")," ",
      paste(sQuote(missing),collapse=","),
      " not found in ",sQuote("params"),call.=FALSE)
  }

  pompLoad(object)

  loglik <- -traj.match.nll(object,ode_control=ode_control)

  ofun <- function (par) {
    params[idx] <- par
    coef(object,transform=transform) <<- params
    ll <- traj.match.nll(object,ode_control=ode_control)
    loglik <<- -ll
    if (is.finite(ll) || is.na(fail.value)) ll else fail.value
  }

  environment(ofun) <- list2env(
    list(object=object,transform=transform,fail.value=fail.value,
      params=params,idx=idx,loglik=loglik,ode_control=ode_control),
    parent=parent.frame(2)
  )

  new("traj_match_objfun",ofun,env=environment(ofun))

}

traj.match.nll <- function (object, seed, ode_control) {
  object@states <- do.call(trajectory,c(list(object),ode_control))
  -sum(dmeasure(object,y=obs(object),x=object@states,
    times=time(object),params=coef(object),
    log=TRUE))
}
