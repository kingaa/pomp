##' Trajectory matching
##'
##' Estimation of parameters for deterministic \acronym{POMP} models via trajectory matching.
##'
##' In trajectory matching, one attempts to minimize the discrepancy between a \acronym{POMP} model's predictions and data under the assumption that the latent state process is deterministic and all discrepancies between model and data are due to measurement error.
##' The measurement model likelihood (\code{dmeasure}), or rather its negative, is the natural measure of the discrepancy.
##'
##' Trajectory matching is a generalization of the traditional nonlinear least squares approach.
##' In particular, if, on some scale, measurement errors are normal with constant variance, then trajectory matching is equivalent to least squares on that particular scale.
##'
##' \code{traj_objfun} constructs an objective function that evaluates the likelihood function.
##' It can be passed to any one of a variety of numerical optimization routines, which will adjust model parameters to minimize the discrepancies between the power spectrum of model simulations and that of the data.
##'
##' @name trajectory matching
##' @rdname traj_match
##' @docType methods
##' @include trajectory.R pomp_class.R workhorses.R
##' @aliases traj_objfun traj_objfun,missing-method traj_objfun,ANY-method
##' @concept trajectory matching
##' @family deterministic methods
##' @family methods based on maximization
##'
##' @inheritParams probe matching
##' @inheritParams trajectory
##' @inheritParams pomp
##'
##' @param \dots additional arguments will modify the model structure
##'
##' @return
##' \code{traj_objfun} constructs a stateful objective function for spectrum matching.
##' Specifically, \code{traj_objfun} returns an object of class \sQuote{traj_match_objfun}, which is a function suitable for use in an \code{\link[stats]{optim}}-like optimizer.
##' In particular, this function takes a single numeric-vector argument that is assumed to contain the parameters named in \code{est}, in that order.
##' When called, it will return the negative log likelihood.
##' It is a stateful function:
##' Each time it is called, it will remember the values of the parameters and its estimate of the log likelihood.
##'
##' @inheritSection pomp Note for Windows users
##' @inheritSection objfun Important Note
##' @inheritSection objfun Warning! Objective functions based on C snippets
##' 
##' @seealso \code{\link[stats]{optim}}, \code{\link[subplex]{subplex}}, \code{\link[nloptr]{nloptr}}
##' 
##' @example examples/traj_match.R
##' 
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
  "traj_objfun",
  function (data, ...)
    standardGeneric("traj_objfun")
)

setMethod(
  "traj_objfun",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("traj_objfun","data")
  }
)

setMethod(
  "traj_objfun",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("traj_objfun",data)
  }
)

##' @rdname traj_match
##' @export
setMethod(
  "traj_objfun",
  signature=signature(data="data.frame"),
  definition=function(data,
    est = character(0), fail.value = NA,
    ode_control = list(),
    params, rinit, skeleton, dmeasure, partrans,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      tmof_internal(
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
      error = function (e) pStop(who="traj_objfun",conditionMessage(e))
    )

  }
)

##' @rdname traj_match
##' @export
setMethod(
  "traj_objfun",
  signature=signature(data="pomp"),
  function (data,
    est = character(0), fail.value = NA, ode_control = list(),
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      tmof_internal(
        data,
        est=est,
        fail.value=fail.value,
        ode_control=ode_control,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop(who="traj_objfun",conditionMessage(e))
    )

  }
)

##' @rdname traj_match
##' @export
setMethod(
  "traj_objfun",
  signature=signature(data="traj_match_objfun"),
  definition=function (data,
    est, fail.value, ode_control,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(est)) est <- data@est
    if (missing(fail.value)) fail.value <- data@env$fail.value
    if (missing(ode_control)) ode_control <- data@env$ode_control

    traj_objfun(
      data@env$object,
      est=est,
      fail.value=fail.value,
      ode_control=ode_control,
      ...,
      verbose=verbose
    )

  }
)

tmof_internal <- function (
  object, est, fail.value, ode_control, ..., verbose
) {

  verbose <- as.logical(verbose)
  object <- pomp(object,...,verbose=verbose)

  if (undefined(object@skeleton) || undefined(object@dmeasure))
    pStop_(paste(sQuote(c("skeleton","dmeasure")),collapse=", "),
      " are needed basic components.")

  fail.value <- as.numeric(fail.value)
  .gnsi <- TRUE

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

  loglik <- traj_match_logLik(object,ode_control=ode_control)

  ofun <- function (par = numeric(0)) {
    params[idx] <- par
    coef(object,transform=TRUE,.gnsi=.gnsi) <<- params
    loglik <<- traj_match_logLik(object,ode_control=ode_control,.gnsi=.gnsi)
    .gnsi <<- FALSE
    if (is.finite(loglik) || is.na(fail.value)) -loglik else fail.value
  }

  environment(ofun) <- list2env(
    list(object=object,fail.value=fail.value,params=params,
      idx=idx,loglik=loglik,ode_control=ode_control,.gnsi=.gnsi),
    parent=parent.frame(2)
  )

  new("traj_match_objfun",ofun,env=environment(ofun),est=est)

}

traj_match_logLik <- function (object, ode_control, .gnsi = TRUE) {
  object@states <- do.call(
    flow,
    c(list(object,x0=rinit(object),.gnsi=.gnsi),ode_control)
  )
  sum(
    dmeasure(
      object,
      y=object@data,
      x=object@states,
      times=object@times,
      params=object@params,
      log=TRUE,
      .gnsi=.gnsi
    )
  )
}

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
