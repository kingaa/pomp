##' Probe matching
##'
##' Estimation of parameters by maximum synthetic likelihood
##'
##' In probe-matching, one attempts to minimize the discrepancy between simulated and actual data, as measured by a set of summary statistics called \emph{probes}.
##' In \pkg{pomp}, this discrepancy is measured using the \dQuote{synthetic likelihood} as defined by Wood (2010).
##'
##' @name probe.match
##' @rdname probe_match
##' @aliases probe.match probe.objfun probe.objfun,missing-method
##' probe.objfun,ANY-method
##' @include probe.R
##' @author Aaron A. King
##' @family summary statistics methods
##' @family \pkg{pomp} parameter estimation methods
##' @seealso \code{\link{optim}} \code{\link[subplex]{subplex}} \code{\link[nloptr]{nloptr}}
##'
##' @param est character vector; the names of parameters to be estimated.
##'
##' @param fail.value optional numeric scalar;
##' if non-\code{NA}, this value is substituted for non-finite values of the objective function.
##' It should be a large number (i.e., bigger than any legitimate values the objective function is likely to take).
##'
##' @param seed  integer.
##' When fitting, it is often best to fix the seed of the random-number generator (RNG).
##' This is accomplished by setting \code{seed} to an integer.
##' By default, \code{seed = NULL}, which does not alter the RNG state.
##'
##' @inheritParams probe
##' @inheritParams pomp
##'
##' @return
##' \code{probe.objfun} constructs a stateful objective function for probe matching.
##' Specifically, \code{probe.objfun} returns an object of class \sQuote{probe_match_objfun}, which is a function suitable for use in an \code{\link{optim}}-like optimizer.
##' In particular, this function takes a single numeric-vector argument that is assumed to contain the parameters named in \code{est}, in that order.
##' When called, it will return the negative synthetic log likelihood for the probes specified.
##' It is a stateful function:
##' Each time it is called, it will remember the values of the parameters and its estimate of the synthetic likelihood.
##'
##' @inheritSection objfun Important Note
##'
NULL

setClass(
  "probe_match_objfun",
  contains="function",
  slots=c(
    env="environment",
    est="character"
  )
)

setGeneric(
  "probe.objfun",
  function (data, ...)
    standardGeneric("probe.objfun")
)

setMethod(
  "probe.objfun",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("probe.objfun","data")
  }
)

setMethod(
  "probe.objfun",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("probe.objfun",data)
  }
)

##' @name probe.objfun-data.frame
##' @aliases probe.objfun,data.frame-method
##' @rdname probe_match
##' @export
setMethod(
  "probe.objfun",
  signature=signature(data="data.frame"),
  definition=function (data,
    est = character(0), fail.value = NA,
    probes, nsim, seed = NULL,
    params, rinit, rprocess, rmeasure, partrans,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      pmof.internal(
        data,
        est=est,
        fail.value=fail.value,
        probes=probes,
        nsim=nsim,
        seed=seed,
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        rmeasure=rmeasure,
        partrans=partrans,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("probe.objfun",conditionMessage(e))
    )

  }
)

##' @name probe.objfun-pomp
##' @aliases probe.objfun,pomp-method
##' @rdname probe_match
##' @export
setMethod(
  "probe.objfun",
  signature=signature(data="pomp"),
  definition=function (data,
    est = character(0), fail.value = NA,
    probes, nsim, seed = NULL,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      pmof.internal(
        data,
        est=est,
        fail.value=fail.value,
        probes=probes,
        nsim=nsim,
        seed=seed,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("probe.objfun",conditionMessage(e))
    )

  }
)

##' @name probe.objfun-probed_pomp
##' @aliases probe.objfun,probed_pomp-method
##' @rdname probe_match
##' @export
setMethod(
  "probe.objfun",
  signature=signature(data="probed_pomp"),
  definition=function (data,
    est = character(0), fail.value = NA,
    probes, nsim, seed = NULL,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(probes)) probes <- data@probes
    if (missing(nsim)) nsim <- data@nsim

    probe.objfun(
      as(data,"pomp"),
      est=est,
      fail.value=fail.value,
      probes=probes,
      nsim=nsim,
      seed=seed,
      ...,
      verbose=verbose
    )

  }
)

##' @name probe.objfun-probe_match_objfun
##' @aliases probe.objfun,probe_match_objfun-method
##' @rdname probe_match
##' @export
setMethod(
  "probe.objfun",
  signature=signature(data="probe_match_objfun"),
  definition=function (data,
    est, fail.value, seed = NULL,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(est)) est <- data@est
    if (missing(fail.value)) fail.value <- data@env$fail.value

    probe.objfun(
      data@env$object,
      est=est,
      fail.value=fail.value,
      seed=seed,
      ...,
      verbose=verbose
    )

  }
)

pmof.internal <- function (object,
  est, fail.value = NA,
  probes, nsim, seed = NULL,
  ..., verbose) {

  verbose <- as.logical(verbose)

  object <- probe(object,probes=probes,nsim=nsim,seed=seed,...,verbose=verbose)

  fail.value <- as.numeric(fail.value)
  loglik <- logLik(object)

  est <- as.character(est)
  est <- est[nzchar(est)]

  params <- coef(object,transform=TRUE)

  idx <- match(est,names(params))
  if (any(is.na(idx))) {
    missing <- est[is.na(idx)]
    pStop_("parameter",ngettext(length(missing),"","s")," ",
      paste(sQuote(missing),collapse=",")," not found in ",sQuote("params"),".")
  }

  pompLoad(object,verbose=verbose)

  ofun <- function (par) {
    params[idx] <- par
    coef(object,transform=TRUE) <<- params
    loglik <<- probe.eval(object)
    if (is.finite(loglik) || is.na(fail.value)) -loglik else fail.value
  }

  environment(ofun) <- list2env(
    list(object=object,fail.value=fail.value,
      params=params,idx=idx,loglik=loglik,seed=seed),
    parent=parent.frame(2)
  )

  new("probe_match_objfun",ofun,env=environment(ofun),est=est)

}

probe.eval <- function (object) {

  ## apply probes to model simulations
  simvals <- tryCatch(
    freeze(
      .Call(P_apply_probe_sim,object=object,nsim=object@nsim,params=object@params,
        probes=object@probes,datval=object@datvals,.gnsi=TRUE),
      seed=object@seed
    ),
    error = function (e) pStop_("applying probes to simulated data: ",conditionMessage(e))
  )

  tryCatch(
    .Call(P_synth_loglik,simvals,object@datvals),
    error = function (e) pStop_("in synthetic likelihood computation: ",conditionMessage(e))
  )

}

##' @name probe-probe_match_obfjun
##' @rdname probe
##' @aliases probe,probe_match_objfun-method
##'
##' @details
##' When \code{probe} operates on a probe-matching objective function (a \sQuote{probe_match_objfun} object), by default, the
##' random-number generator seed is fixed at the value given when the objective function was constructed.
##' Specifying \code{NULL} or an integer for \code{seed} overrides this behavior.
##'
##' @export
setMethod(
  "probe",
  signature=signature(data="probe_match_objfun"),
  definition=function (data, seed,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(seed)) seed <- data@env$seed

    probe(
      data@env$object,
      seed=seed,
      ...,
      verbose=verbose
    )

  }
)

##' @name coerce-probe_match_objfun-probed_pomp
##' @aliases coerce,probe_match_objfun,probed_pomp-method
##' @rdname probe
##'
setAs(
  from="probe_match_objfun",
  to="probed_pomp",
  def = function (from) {
    from@env$object
  }
)
