##' Probes (AKA summary statistics)
##'
##' Probe a partially-observed Markov process by computing summary statistics
##' and the synthetic likelihood.
##'
##' \code{probe} applies one or more \dQuote{probes} to time series data and
##' model simulations and compares the results.  It can be used to diagnose
##' goodness of fit and/or as the basis for \dQuote{probe-matching}, a
##' generalized method-of-moments approach to parameter estimation.
##'
##' A call to \code{probe} results in the evaluation of the probe(s) in
##' \code{probes} on the data.  Additionally, \code{nsim} simulated data sets
##' are generated (via a call to \code{\link[=simulate-pomp]{simulate}}) and
##' the probe(s) are applied to each of these.  The results of the probe
##' computations on real and simulated data are stored in an object of class
##' \sQuote{probed_pomp}.
##'
##' @name Probes
##' @docType methods
##' @rdname probe
##' @include pomp_class.R pomp_fun.R
##' @aliases probe probe,missing-method probe,ANY-method
##' @family summary statistics
##'
##' @section Methods:
##' \describe{
##' \item{plot}{ displays diagnostic plots.  }
##' \item{summary}{ displays summary information.  The summary includes
##' quantiles (fractions of simulations with probe values less than those
##' realized on the data) and the corresponding two-sided p-values.  In
##' addition, the \dQuote{synthetic likelihood} (Wood 2010) is computed, under
##' the assumption that the probe values are multivariate-normally distributed.
##' }
##' \item{probevals}{ extracts the realized values of the probes on the data
##' and on the simulations.  These are returned in a list of two elements,
##' \code{datvals} and \code{simvals}.  }
##' \item{logLik}{ returns the synthetic
##' likelihood for the probes.  NB: in general, this is not the same as the
##' likelihood.  }
##' \item{as.data.frame, as(object,"data.frame")}{ coercing a
##' \sQuote{probed_pomp} to a \sQuote{data.frame}, gives the realized values of
##' the probes on the data and on the simulations.  The variable \code{.id}
##' indicates whether the probes are from the data or simulations.  }
##' }
##'
##' @return
##' \code{probe} returns an object of class \sQuote{probed_pomp}, which
##' is derived from the \sQuote{pomp} class and contains additional information
##' about the \code{probe} calculation.
##'
##' This information can be summarized via a call to \code{summary} and
##' displayed graphically via a call to \code{plot}.
##'
##' @author Daniel C. Reuman, Aaron A. King
##'
##' @seealso \link{Basic probes}, \code{\link{spect}}, and the
##' tutorials on the \href{https://kingaa.github.io/pomp/}{package website}.
##'
##' @references
##' B. E. Kendall, C. J. Briggs, W. M. Murdoch, P. Turchin, S. P.
##' Ellner, E. McCauley, R. M. Nisbet, S. N. Wood Why do populations cycle? A
##' synthesis of statistical and mechanistic modeling approaches, Ecology,
##' 80:1789--1805, 1999.
##'
##' S. N. Wood Statistical inference for noisy nonlinear ecological dynamic
##' systems, Nature, 466: 1102--1104, 2010.
NULL

setClass(
  "probed_pomp",
  contains="pomp",
  slots=c(
    probes="list",
    nsim="integer",
    datvals="numeric",
    simvals="array",
    quantiles="numeric",
    pvals="numeric",
    synth.loglik="numeric",
    seed="integer"
  )
)

setGeneric(
  "probe",
  function (object, ...)
    standardGeneric("probe")
)

##' @name probe-pomp
##' @aliases probe,pomp-method
##' @rdname probe
##'
##' @param object an object of class \sQuote{pomp}, or of a class the extends \sQuote{pomp}
##' @param probes a single probe or a list of one or more probes.
##' A probe is simply a scalar- or vector-valued function of one argument that can be applied to the data array of a \sQuote{pomp}.
##' A vector-valued probe must always return a vector of the same size.
##' A number of useful probes are provided with the package:
##' see \link{Basic probes}.
##' @param params optional named numeric vector of model parameters.
##' By default, \code{params=coef(object)}.
##' @param nsim the number of model simulations to be computed.
##' @param seed optional integer;
##' if non-\code{NULL}, the random number generator will be initialized with this seed for simulations.
##' See \code{\link[=simulate-pomp]{simulate}}.
##' @param verbose logical; print diagnostic messages?
##' @param \dots Additional arguments.  In the case of \code{probe} and
##' \code{probe.match.objfun}, additional arguments are passed to
##' \code{\link{pomp}}, allowing one to supply new or modify existing model
##' characteristics or components.
##'
setMethod(
  "probe",
  signature=signature(object="pomp"),
  definition=function (object, probes, params, nsim, seed = NULL, ...,
    verbose = getOption("verbose", FALSE))
  {

    if (missing(probes)) probes <- NULL
    if (missing(params)) params <- coef(object)
    if (missing(nsim)) nsim <- NULL

    probe.internal(
      object=object,
      probes=probes,
      params=params,
      nsim=nsim,
      seed=seed,
      ...,
      verbose=verbose
    )
  }
)

##' @name probe-probed_pomp
##' @aliases probe,probed_pomp-method
##' @rdname probe
setMethod(
  "probe",
  signature=signature(object="probed_pomp"),
  definition=function (object, probes, nsim, seed = NULL, ...,
    verbose = getOption("verbose", FALSE)) {

    if (missing(probes)) probes <- object@probes
    if (missing(nsim)) nsim <- object@nsim

    probe(
      as(object,"pomp"),
      probes=probes,
      nsim=nsim,
      seed=seed,
      ...,
      verbose=verbose
    )
  }
)

setMethod(
  "probe",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("probe"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "probe",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("probe")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

probe.internal <- function (object, probes, params, nsim, seed,
  .getnativesymbolinfo = TRUE, ..., verbose) {

  ep <- paste0("in ",sQuote("probe"),": ")
  verbose <- as.logical(verbose)

  object <- pomp(object,...)

  if (is.null(probes))
    stop(ep,sQuote("probes")," must be furnished.",call.=FALSE)
  if (!is.list(probes)) probes <- list(probes)
  if (!all(sapply(probes,is.function)))
    stop(ep,sQuote("probes")," must be a function or a list of functions.",
      call.=FALSE)
  if (!all(sapply(probes,function(f)length(formals(f))==1)))
    stop(ep,"each probe must be a function of a single argument.",call.=FALSE)

  if (is.list(params)) params <- unlist(params)
  if (is.null(params)) params <- numeric(0)
  if (length(params)==0)
    stop(ep,sQuote("params")," must be specified.",call.=FALSE)
  if (!is.numeric(params) || is.null(names(params)))
    stop(ep,sQuote("params")," must be furnished as a named numeric vector.",
      call.=FALSE)

  nsim <- as.integer(nsim)
  if (length(nsim) < 1)
    stop(ep,sQuote("nsim")," must be specified.",call.=FALSE)
  if (length(nsim) > 1 || !is.finite(nsim) || nsim <= 0)
    stop(ep,"number of simulations, ",sQuote("nsim"),
      ", must be a single positive integer.",call.=FALSE)

  seed <- as.integer(seed)

  gnsi <- as.logical(.getnativesymbolinfo)

  ## apply probes to data
  datval <- tryCatch(
    .Call(apply_probe_data,object,probes),
    error = function (e) {
      stop(ep,"applying probes to actual data: ",
        conditionMessage(e),call.=FALSE)
    }
  )

  nprobes <- length(datval)

  if (nprobes >= nsim)
    stop(ep,sQuote("nsim")," (=",nsim,"), should be (much) larger than the ",
      "number of probes (=",nprobes,").",call.=FALSE)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  ## apply probes to model simulations
  simval <- tryCatch(
    .Call(
      apply_probe_sim,
      object=object,
      nsim=nsim,
      params=params,
      seed=seed,
      probes=probes,
      datval=datval,
      gnsi=gnsi
    ),
    error = function (e) {
      stop(ep,"applying probes to simulated data: ",
        conditionMessage(e),call.=FALSE)
    }
  )

  pvals <- numeric(nprobes)
  names(pvals) <- names(datval)
  quants <- numeric(nprobes)
  names(quants) <- names(datval)

  for (k in seq_len(nprobes)) {
    r <- min(sum(simval[,k]>datval[k]),sum(simval[,k]<datval[k]))
    tails <- (r+1)/(nsim+1)
    pvals[k] <- min(2*tails,1)
    quants[k] <- sum(simval[,k]<datval[k])/nsim
  }

  ll <- tryCatch(
    .Call(synth_loglik,simval,datval),
    error = function (e) {
      stop(ep,"in synthetic likelihood computation: ",  # nocov
        conditionMessage(e),call.=FALSE)                # nocov
    }
  )

  coef(object) <- params

  names(dimnames(simval)) <- c("rep","probe")

  new(
    "probed_pomp",
    object,
    probes=probes,
    nsim=nsim,
    datvals=datval,
    simvals=simval,
    quantiles=quants,
    pvals=pvals,
    synth.loglik=ll,
    seed=seed
  )
}
