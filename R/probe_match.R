##' Probe matching
##'
##' Estimation of parameters by maximum synthetic likelihood
##'
##' In probe-matching, one attempts to minimize the discrepancy between simulated and actual data, as measured by a set of summary statistics called \emph{probes}.
##' In \pkg{pomp}, this discrepancy is measured using the \dQuote{synthetic likelihood} as defined by Wood (2010).
##'
##'
##' @name probe.match
##' @rdname probe_match
##' @keywords optimize
##' @aliases probe.match probe.match.objfun probe.match.objfun,missing-method
##' probe.match.objfun,ANY-method
##' @include probe.R loglik.R summary.R
##' @family summary statistics methods
##'
##' @return
##' \code{probe.match.objfun} construct a stateful objective function for probe matching.
##' Specfically, \code{probe.match.objfun} returns an object of class \sQuote{probe_match_objfun}, which is a function suitable for use in an \code{\link{optim}}-like optimizer.
##' In particular, this function takes a single numeric-vector argument that is assumed to contain the parameters named in \code{est}, in that order.
##' When called, it wil return the negative synthetic log likelihood for the probes specified.
##' It is a stateful function:
##' Each time it is called, it will remember the values of the parameters and its estimate of the synthetic likelihood.
##'
##' @section Important Note:
##' Since \pkg{pomp} cannot guarantee that the \emph{final} call an optimizer makes to the function is a call \emph{at} the optimum, it cannot guarantee that the parameters stored in the function are the optimal ones.
##' One should check that the parameters agree with those that are returned by the optimizer.
##' The best practice is to call \code{\link[=probe,probe_match_objfun-method]{probe}} on the objective function after the optimization has been performed, thus obtaining a \sQuote{probed_pomp} object containing the (putative) optimal parameters and synthetic likelihood.
##' @seealso \code{\link{trajectory}}, \code{\link{optim}},
##' \code{\link[subplex]{subplex}}, \code{\link[nloptr]{nloptr}}
##'
NULL

setClass(
  "probe_match_objfun",
  contains="function",
  slots=c(
    env="environment"
  )
)

setGeneric(
  "probe.match.objfun",
  function (object, ...)
    standardGeneric("probe.match.objfun")
)

setMethod(
  "probe.match.objfun",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("probe.match.objfun"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "probe.match.objfun",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("probe.match.objfun")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' @name probe.match.objfun-pomp
##' @aliases probe.match.objfun,pomp-method
##' @rdname probe_match
##'
##' @inheritParams probe-pomp
##' @param est character vector; the names of parameters to be estimated.
##' @param fail.value optional numeric scalar;
##' if non-\code{NA}, this value is substituted for non-finite values of the objective function.
##' It should be a large number (i.e., bigger than any legitimate values the objective function is likely to take).
##' @param transform logical;
##' if \code{TRUE}, optimization is to be performed on the transformed scale.
##'
setMethod(
  "probe.match.objfun",
  signature=signature(object="pomp"),
  definition=function (object, params, est, probes,
    nsim, seed = NULL, fail.value = NA,
    transform = FALSE, ...) {

    object <- probe(object,probes=probes,nsim=nsim,
      seed=seed,params=params,...)

    probe.match.objfun(
      object=object,
      est=est,
      fail.value=fail.value,
      transform=transform
    )

  }
)

##' @name probe.match.objfun-probed_pomp
##' @aliases probe.match.objfun,probed_pomp-method
##' @rdname probe_match
setMethod(
  "probe.match.objfun",
  signature=signature(object="probed_pomp"),
  definition=function (object, params, est, probes, nsim, seed = NULL,
    fail.value = NA, transform = FALSE, ...) {

    if (missing(probes)) probes <- object@probes
    if (missing(nsim)) nsim <- nrow(object@simvals)
    if (missing(seed)) seed <- object@seed

    pmof.internal(
      object=object,
      params=params,
      est=est,
      probes=probes,
      nsim=nsim,
      seed=seed,
      fail.value=fail.value,
      transform=transform,
      ...
    )

  }
)

##' @name probe.match.objfun-probe_match_objfun
##' @aliases probe.match.objfun,probe_match_objfun-method
##' @rdname probe_match
setMethod(
  "probe.match.objfun",
  signature=signature(object="probe_match_objfun"),
  definition=function (object, ...) {
    probe.match.objfun(object@env$object,...)
  }
)

pmof.internal <- function (object, params, est, probes, nsim, seed = NULL,
  fail.value = NA, transform = FALSE, ...)
{

  ep <- paste0("in ",sQuote("probe.match.objfun"),": ")

  transform <- as.logical(transform)
  fail.value <- as.numeric(fail.value)
  loglik <- logLik(object)

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

  ofun <- function (par) {
    params[idx] <- par
    coef(object,transform=transform) <<- params
    loglik <<- probe.eval(object)
    if (is.finite(loglik) || is.na(fail.value)) -loglik else fail.value
  }

  environment(ofun) <- list2env(
    list(object=object,transform=transform,fail.value=fail.value,
      params=params,idx=idx,loglik=loglik,seed=seed),
    parent=parent.frame(2)
  )

  new("probe_match_objfun",ofun,env=environment(ofun))

}

probe.eval <- function (object) {

  ep <- "in objective function: "
  ## apply probes to model simulations
  simvals <- tryCatch(
    .Call(
      apply_probe_sim,
      object=object,
      nsim=object@nsim,
      params=object@params,
      seed=object@seed,
      probes=object@probes,
      datval=object@datvals,
      .getNativeSymbolInfo=TRUE
    ),
    error = function (e) {
      stop(ep,"applying probes to simulated data: ",
        conditionMessage(e),call.=FALSE)
    }
  )

  loglik <- tryCatch(
    .Call(synth_loglik,simvals,object@datvals),
    error = function (e) {
      stop(ep,"in synthetic likelihood computation: ",
        conditionMessage(e),call.=FALSE)
    }
  )

  loglik
}

##' @name probe-probe_match_obfjun
##' @rdname probe
##' @aliases probe,probe_match_objfun-method
setMethod(
  "probe",
  signature=signature(object="probe_match_objfun"),
  definition=function (object, ...) {
    probe(object@env$object,...)
  }
)

##' @name summary-probe_match_obfjun
##' @rdname summary
##' @aliases summary,probe_match_objfun-method
setMethod(
  "summary",
  signature=signature(object="probe_match_objfun"),
  definition=function (object) {
    summary(object@env$object)
  }
)

##' @name logLik-probe_match_obfjun
##' @rdname loglik
##' @aliases logLik,probe_match_objfun-method
setMethod(
  "logLik",
  signature=signature(object="probe_match_objfun"),
  definition=function (object) {
    logLik(object@env$object)
  }
)

##' @name coef-probe_match_objfun
##' @rdname coef
##' @aliases coef,probe_match_objfun-method
setMethod(
  "coef",
  signature=signature(object="probe_match_objfun"),
  definition=function (object, ...) {
    coef(object@env$object,...)
  }
)
