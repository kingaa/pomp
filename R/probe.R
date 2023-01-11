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
##' are generated (via a call to \code{\link{simulate}}) and
##' the probe(s) are applied to each of these.  The results of the probe
##' computations on real and simulated data are stored in an object of class
##' \sQuote{probed_pomp}.
##'
##' @docType methods
##' @name probe
##' @rdname probe
##' @aliases probe,missing-method probe,ANY-method
##' @author Daniel C. Reuman, Aaron A. King
##' @concept synthetic likelihood
##' @family elementary algorithms
##' @family summary statistic-based methods
##'
##' @inheritSection pomp Note for Windows users
##' 
##' @include pomp_class.R pomp_fun.R pomp.R
##' @importFrom stats quantile
##'
##' @param probes a single probe or a list of one or more probes.
##' A probe is simply a scalar- or vector-valued function of one argument that can be applied to the data array of a \sQuote{pomp}.
##' A vector-valued probe must always return a vector of the same size.
##' A number of useful probes are provided with the package:
##' see \link[=basic probes]{basic probes}.
##' @param nsim the number of model simulations to be computed.
##' @param seed optional integer;
##' if non-\code{NULL}, the random number generator will be initialized with this seed for simulations.
##' See \code{\link{simulate}}.
##' @inheritParams pomp
##'
##' @return
##' \code{probe} returns an object of class \sQuote{probed_pomp}, which contains the data and the model, together with the results of the \code{probe} calculation.
##'
##' @section Methods:
##' The following methods are available.
##' \describe{
##' \item{\code{plot}}{ displays diagnostic plots.  }
##' \item{\code{summary}}{ displays summary information.
##' The summary includes quantiles (fractions of simulations with probe values less than those realized on the data) and the corresponding two-sided p-values.
##' In addition, the \dQuote{synthetic likelihood} (Wood 2010) is computed,
##' under the assumption that the probe values are multivariate-normally distributed.  }
##' \item{\code{logLik}}{ returns the synthetic likelihood for the probes.
##' NB: in general, this is not the same as the likelihood.  }
##' \item{\code{as.data.frame}}{
##'  coerces a \sQuote{probed_pomp} to a \sQuote{data.frame}.
##'  The latter contains the realized values of the probes on the data and on the simulations.
##' The variable \code{.id} indicates whether the probes are from the data or simulations.  }
##' }
##'
##' @references
##'
##' \Kendall1999
##'
##' \Wood2010
##' 
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
  function (data, ...)
    standardGeneric("probe")
)

setMethod(
  "probe",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("probe","data")
  }
)

setMethod(
  "probe",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("probe",data)
  }
)

##' @rdname probe
##' @export
setMethod(
  "probe",
  signature=signature(data="data.frame"),
  definition=function (data,
    probes, nsim, seed = NULL,
    params, rinit, rprocess, rmeasure,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      probe_internal(data,probes=probes,nsim=nsim,seed=seed,
        rinit=rinit,rprocess=rprocess,rmeasure=rmeasure,params=params,
        ...,verbose=verbose),
      error = function (e) pStop("probe",conditionMessage(e))
    )


  }
)

##' @rdname probe
##' @export
setMethod(
  "probe",
  signature=signature(data="pomp"),
  definition=function (data,
    probes, nsim, seed = NULL, ...,
    verbose = getOption("verbose", FALSE))
  {

    tryCatch(
      probe_internal(data,probes=probes,nsim=nsim,seed=seed,
        ...,verbose=verbose),
      error = function (e) pStop("probe",conditionMessage(e))
    )

  }
)

##' @rdname probe
##' @export
setMethod(
  "probe",
  signature=signature(data="probed_pomp"),
  definition=function (data,
    probes, nsim, seed = NULL, ...,
    verbose = getOption("verbose", FALSE)) {

    if (missing(probes)) probes <- data@probes
    if (missing(nsim)) nsim <- data@nsim

    probe(as(data,"pomp"),probes=probes,nsim=nsim,seed=seed,
      ...,verbose=verbose)

  }
)

probe_internal <- function (object, probes, nsim, seed, ...,
  .gnsi = TRUE, verbose) {

  verbose <- as.logical(verbose)

  object <- pomp(object,...,verbose=verbose)

  if (undefined(object@rprocess) || undefined(object@rmeasure))
    pStop_(paste(sQuote(c("rprocess","rmeasure")),collapse=", ")," are needed basic components.")

  if (missing(probes)) probes <- NULL
  if (missing(nsim)) nsim <- NULL

  if (is.null(probes)) pStop_(sQuote("probes")," must be furnished.")
  if (!is.list(probes)) probes <- list(probes)
  if (!all(vapply(probes,is.function,logical(1L))))
    pStop_(sQuote("probes")," must be a function or a list of functions.")
  if (!all(vapply(probes,\(f)length(formals(f))==1L,logical(1L))))
    pStop_("each probe must be a function of a single argument.")

  nsim <- as.integer(nsim)
  if (length(nsim) < 1) pStop_(sQuote("nsim")," must be specified.")
  if (length(nsim) > 1 || !is.finite(nsim) || nsim <= 0)
    pStop_("number of simulations, ",sQuote("nsim"),", must be a single positive integer.")

  seed <- as.integer(seed)
  gnsi <- as.logical(.gnsi)

  params <- coef(object)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  ## apply probes to data
  datval <- tryCatch(
    .Call(P_apply_probe_data,object,probes),
    error = function (e)
      pStop_("applying probes to actual data: ",conditionMessage(e))
  )

  nprobes <- length(datval)

  if (nprobes >= nsim)
    pStop_(sQuote("nsim")," (=",nsim,"), should be (much) larger than the ", "number of probes (=",nprobes,").")

  ## apply probes to model simulations
  simval <- tryCatch(
    freeze(
      .Call(P_apply_probe_sim,object=object,nsim=nsim,params=params,probes=probes,
        datval=datval,gnsi=gnsi),
      seed=seed
    ),
    error = function (e)
      pStop_("applying probes to simulated data: ",conditionMessage(e))
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
    .Call(P_synth_loglik,simval,datval),
    error = function (e)
      pStop_("in synthetic likelihood computation: ",conditionMessage(e))
  )

  names(dimnames(simval)) <- c(".id","probe")

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

##' @rdname summary
##' @include summary.R
##' @export
setMethod(
  "summary",
  signature=signature(object="probed_pomp"),
  definition=function (object, ...) {
    list(
      coef=coef(object),
      nsim=nrow(object@simvals),
      quantiles=object@quantiles,
      pvals=object@pvals,
      synth.loglik=object@synth.loglik
    )
  }
)
