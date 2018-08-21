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
##' @include probe.R
##' @author Aaron A. King
##' @family summary statistics methods
##' @family \pkg{pomp} parameter estimation methods
##' @seealso \code{\link{optim}} \code{\link[subplex]{subplex}} \code{\link[nloptr]{nloptr}}
##'
##' @inheritParams probe
##' @inheritParams pomp
##' @param est character vector; the names of parameters to be estimated.
##' @param fail.value optional numeric scalar;
##' if non-\code{NA}, this value is substituted for non-finite values of the objective function.
##' It should be a large number (i.e., bigger than any legitimate values the objective function is likely to take).
##'
##' @return
##' \code{probe.match.objfun} construct a stateful objective function for probe matching.
##' Specfically, \code{probe.match.objfun} returns an object of class \sQuote{probe_match_objfun}, which is a function suitable for use in an \code{\link{optim}}-like optimizer.
##' In particular, this function takes a single numeric-vector argument that is assumed to contain the parameters named in \code{est}, in that order.
##' When called, it wil return the negative synthetic log likelihood for the probes specified.
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
    env="environment"
  )
)

setGeneric(
  "probe.match.objfun",
  function (data, ...)
    standardGeneric("probe.match.objfun")
)

setMethod(
  "probe.match.objfun",
  signature=signature(data="missing"),
  definition=function (...) {
    stop("in ",sQuote("probe.match.objfun"),": ",sQuote("data"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "probe.match.objfun",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    stop(sQuote("probe.match.objfun")," is not defined for objects of class ",
      sQuote(class(data)),call.=FALSE)
  }
)

##' @name probe.match.objfun-data.frame
##' @aliases probe.match.objfun,data.frame-method
##' @rdname probe_match
##' @export
setMethod(
  "probe.match.objfun",
  signature=signature(data="data.frame"),
  definition=function (data,
    rinit, rprocess, rmeasure, partrans, params,
    est = character(0), probes, nsim, seed = NULL, fail.value = NA, ...,
    verbose = getOption("verbose", FALSE)) {

    data <- tryCatch(
      pomp(data,rinit=rinit,rprocess=rprocess,rmeasure=rmeasure,
        partrans=partrans,params=params,...,verbose=verbose),
      error = function (e) pomp_stop(conditionMessage(e))
    )

    probe.match.objfun(data,est=est,probes=probes,nsim=nsim,seed=seed,
      fail.value=fail.value,verbose=verbose)

  }
)

##' @name probe.match.objfun-pomp
##' @aliases probe.match.objfun,pomp-method
##' @rdname probe_match
##' @export
setMethod(
  "probe.match.objfun",
  signature=signature(data="pomp"),
  definition=function (data, est = character(0), probes, nsim, seed = NULL,
    fail.value = NA, ..., verbose = getOption("verbose", FALSE)) {

    pmof.internal(data,est=est,probes=probes,nsim=nsim,seed=seed,
      fail.value=fail.value,...,verbose=verbose)

  }
)

##' @name probe.match.objfun-probed_pomp
##' @aliases probe.match.objfun,probed_pomp-method
##' @rdname probe_match
##' @export
setMethod(
  "probe.match.objfun",
  signature=signature(data="probed_pomp"),
  definition=function (data, est = character(0), probes, nsim, seed = NULL,
    fail.value = NA, ..., verbose = getOption("verbose", FALSE)) {

    if (missing(probes)) probes <- data@probes
    if (missing(nsim)) nsim <- data@nsim

    probe.match.objfun(as(data,"pomp"),est=est,probes=probes,nsim=nsim,
      seed=seed,fail.value=fail.value,...,verbose=verbose)

  }
)

##' @name probe.match.objfun-probe_match_objfun
##' @aliases probe.match.objfun,probe_match_objfun-method
##' @rdname probe_match
##' @export
setMethod(
  "probe.match.objfun",
  signature=signature(data="probe_match_objfun"),
  definition=function (data, ..., verbose = getOption("verbose", FALSE)) {
    probe.match.objfun(data@env$object,...,verbose=verbose)
  }
)

pmof.internal <- function (object, est, probes, nsim, seed = NULL,
  fail.value = NA, ..., verbose) {

  object <- tryCatch(
    probe(object,probes=probes,nsim=nsim,seed=seed,...),
    error = function (e) pomp_stop(conditionMessage(e))
  )

  fail.value <- as.numeric(fail.value)
  loglik <- logLik(object)

  est <- as.character(est)
  est <- est[nzchar(est)]

  params <- coef(object,transform=TRUE)

  idx <- match(est,names(params))
  if (any(is.na(idx))) {
    missing <- est[is.na(idx)]
    pomp_stop(ngettext(length(missing),"parameter","parameters")," ",
      paste(sQuote(missing),collapse=",")," not found in ",sQuote("params"),".")
  }

  pompLoad(object)

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
      .getnativesymbolinfo=TRUE
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
##' @export
setMethod(
  "probe",
  signature=signature(data="probe_match_objfun"),
  definition=function (data, ...,
    verbose = getOption("verbose", FALSE)) {
    probe(data@env$object,...,seed=data@env$seed,verbose=verbose)
  }
)
