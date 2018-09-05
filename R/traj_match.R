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
##' @inheritParams probe.match
##' @inheritParams trajectory
##' @inheritParams pomp
##'
##' @param ode_control optional list;
##' the elements of this list will be passed to \code{\link[=deSolve]{ode}}.
##' @param \dots additional arguments will modify the model structure
##'
##' @return
##' \code{traj.match.objfun} constructs a stateful objective function for spectrum matching.
##' Specifically, \code{traj.match.objfun} returns an object of class \sQuote{traj_match_objfun}, which is a function suitable for use in an \code{\link{optim}}-like optimizer.
##' In particular, this function takes a single numeric-vector argument that is assumed to contain the parameters named in \code{est}, in that order.
##' When called, it will return the negative log likelihood.
##' It is a stateful function:
##' Each time it is called, it will remember the values of the parameters and its estimate of the log likelihood.
##'
##' @inheritSection objfun Important Note
##' @seealso \code{\link{trajectory}}, \code{\link{optim}},
##' \code{\link[subplex]{subplex}}, \code{\link[nloptr]{nloptr}}
NULL

setClass(
  "traj_match_objfun",
  contains="function",
  slots=c(
    env="environment",
    est="character"
  )
)

setGeneric(
  "traj.match.objfun",
  function (data, ...)
    standardGeneric("traj.match.objfun")
)

setMethod(
  "traj.match.objfun",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("traj.match.objfun","data")
  }
)

setMethod(
  "traj.match.objfun",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("traj.match.objfun",data)
  }
)

##' @name traj.match.objfun-data.frame
##' @aliases traj.match.objfun,data.frame-method
##' @rdname traj_match
##' @export
setMethod(
  "traj.match.objfun",
  signature=signature(data="data.frame"),
  definition=function(data,
    est = character(0), fail.value = NA,
    ode_control = list(),
    params, rinit, skeleton, dmeasure, partrans,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      tmof.internal(
        data,
        est=est,
        fail.value=fail.value,
        ode_control=ode_control,
        params=params,
        rinit=rinit,
        skeleton=skeleton,
        dmeasure=dmeasure,
        partrans=partrans,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("traj.match.objfun",conditionMessage(e))
    )

  }
)

##' @name traj.match.objfun-pomp
##' @aliases traj.match.objfun traj.match.objfun,pomp-method
##' @rdname traj_match
##' @export
setMethod(
  "traj.match.objfun",
  signature=signature(data="pomp"),
  function (data,
    est, fail.value = NA, ode_control = list(),
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      tmof.internal(
        data,
        est=est,
        fail.value=fail.value,
        ode_control=ode_control,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("traj.match.objfun",conditionMessage(e))
    )

  }
)

##' @name traj.match.objfun-traj_match_objfun
##' @aliases traj.match.objfun,traj_match_objfun-method
##' @rdname traj_match
##' @export
setMethod(
  "traj.match.objfun",
  signature=signature(data="traj_match_objfun"),
  definition=function (data,
    est, fail.value, ode_control,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(est)) est <- data@est
    if (missing(fail.value)) fail.value <- data@env$fail.value
    if (missing(ode_control)) ode_control <- data@env$ode_control

    traj.match.objfun(
      data@env$object,
      est=est,
      fail.value=fail.value,
      ode_control=ode_control,
      ...,
      verbose=verbose
    )

  }
)

tmof.internal <- function (object,
  est, fail.value, ode_control,
  ..., verbose) {

  verbose <- as.logical(verbose)
  object <- pomp(object,...,verbose=verbose)

  if (undefined(object@skeleton) || undefined(object@dmeasure))
    pStop_(paste(sQuote(c("skeleton","dmeasure")),collapse=", ")," are needed basic components.")

  fail.value <- as.numeric(fail.value)

  if (missing(est)) est <- character(0)
  est <- as.character(est)
  est <- est[nzchar(est)]

  params <- coef(object,transform=TRUE)

  idx <- match(est,names(params))
  if (any(is.na(idx))) {
    missing <- est[is.na(idx)]
    pStop_("parameter",ngettext(length(missing),"","s")," ",
      paste(sQuote(missing),collapse=",")," not found in ",sQuote("params"))
  }

  pompLoad(object,verbose=verbose)

  loglik <- traj.match.loglik(object,ode_control=ode_control)

  ofun <- function (par) {
    params[idx] <- par
    coef(object,transform=TRUE) <<- params
    loglik <<- traj.match.loglik(object,ode_control=ode_control)
    if (is.finite(loglik) || is.na(fail.value)) -loglik else fail.value
  }

  environment(ofun) <- list2env(
    list(object=object,fail.value=fail.value,params=params,
      idx=idx,loglik=loglik,ode_control=ode_control),
    parent=parent.frame(2)
  )

  new("traj_match_objfun",ofun,env=environment(ofun),est=est)

}

traj.match.loglik <- function (object, seed, ode_control) {
  object@states <- do.call(trajectory,c(list(object),ode_control))
  sum(dmeasure(object,y=obs(object),x=object@states,
    times=time(object),params=coef(object),
    log=TRUE))
}

##' @name trajectory-traj_match_objfun
##' @aliases trajectory,traj_match_objfun-method
##' @rdname trajectory
##' @export
setMethod(
  "trajectory",
  signature=signature(object="traj_match_objfun"),
  definition=function (object,
    ..., verbose = getOption("verbose", FALSE)) {

    trajectory(
      object@env$object,
      ...,
      verbose=verbose
    )

  }
)
