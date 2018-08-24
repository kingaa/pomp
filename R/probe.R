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
##' @name probe
##' @docType methods
##' @rdname probe
##' @aliases probe probe,missing-method probe,ANY-method
##' @author Daniel C. Reuman, Aaron A. King
##' @family elementary POMP methods
##' @family summary statistics methods
##'
##' @include pomp_class.R pomp_fun.R pomp.R
##' @importFrom stats quantile
##'
##' @inheritParams pomp
##' @param probes a single probe or a list of one or more probes.
##' A probe is simply a scalar- or vector-valued function of one argument that can be applied to the data array of a \sQuote{pomp}.
##' A vector-valued probe must always return a vector of the same size.
##' A number of useful probes are provided with the package:
##' see \link[=basic_probes]{basic probes}.
##' @param nsim the number of model simulations to be computed.
##' @param seed optional integer;
##' if non-\code{NULL}, the random number generator will be initialized with this seed for simulations.
##' See \code{\link[=simulate-pomp]{simulate}}.
##'
##' @return
##' \code{probe} returns an object of class \sQuote{probed_pomp}, which contains the data and the model, together with the results of the \code{probe} calculation.
##'
##' @section Methods:
##' The following methods are available.
##' \describe{
##' \item{plot}{ displays diagnostic plots.  }
##' \item{summary}{ displays summary information.
##' The summary includes quantiles (fractions of simulations with probe values less than those realized on the data) and the corresponding two-sided p-values.
##' In addition, the \dQuote{synthetic likelihood} (Wood 2010) is computed,
##' under the assumption that the probe values are multivariate-normally distributed.  }
##' \item{logLik}{ returns the synthetic likelihood for the probes.
##' NB: in general, this is not the same as the likelihood.  }
##' \item{as.data.frame}{
##'  coerces a \sQuote{probed_pomp} to a \sQuote{data.frame}.
##'  The latter contains the realized values of the probes on the data and on the simulations.
##' The variable \code{.id} indicates whether the probes are from the data or simulations.  }
##' }
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
  function (data, ...)
    standardGeneric("probe")
)

setMethod(
  "probe",
  signature=signature(data="missing"),
  definition=function (...) {
    stop("in ",sQuote("probe"),": ",sQuote("data"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "probe",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    stop(sQuote("probe")," is not defined when ",sQuote("data")," is of class ",
      sQuote(class(data)),call.=FALSE)
  }
)

##' @name probe-data.frame
##' @aliases probe,data.frame-method
##' @rdname probe
##' @export
setMethod(
  "probe",
  signature=signature(data="data.frame"),
  definition=function (data, rinit, rprocess, rmeasure, params,
    probes, nsim, seed = NULL, ...,
    verbose = getOption("verbose", FALSE)) {

    if (missing(rprocess) || missing(rmeasure) || missing(params))
      pStop("probe",paste(sQuote(c("rprocess","rmeasure","params")),
        collapse=", ")," are required arguments.")

    object <- tryCatch(
      pomp(data,rinit=rinit,rprocess=rprocess,rmeasure=rmeasure,
        params=params,...,verbose=verbose),
      error = function (e)
        pStop("probe",conditionMessage(e))
    )

    probe(object,probes=probes,nsim=nsim,seed=seed,verbose=verbose)

  }
)

##' @name probe-pomp
##' @aliases probe,pomp-method
##' @rdname probe
##' @export
setMethod(
  "probe",
  signature=signature(data="pomp"),
  definition=function (data, probes, nsim, seed = NULL, ...,
    verbose = getOption("verbose", FALSE))
  {

    if (missing(probes)) probes <- NULL
    if (missing(nsim)) nsim <- NULL

    probe.internal(data,probes=probes,nsim=nsim,seed=seed,...,verbose=verbose)

  }
)

##' @name probe-probed_pomp
##' @aliases probe,probed_pomp-method
##' @rdname probe
##' @export
setMethod(
  "probe",
  signature=signature(data="probed_pomp"),
  definition=function (data, probes, nsim, seed = NULL, ...,
    verbose = getOption("verbose", FALSE)) {

    if (missing(probes)) probes <- data@probes
    if (missing(nsim)) nsim <- data@nsim

    probe(as(data,"pomp"),probes=probes,nsim=nsim,seed=seed,...,verbose=verbose)

  }
)

probe.internal <- function (object, probes, nsim, seed, ...,
  .getnativesymbolinfo = TRUE, verbose) {

  verbose <- as.logical(verbose)

  object <- tryCatch(
    pomp(object,...),
    error = function (e) pStop("probe",conditionMessage(e))
  )

  if (is.null(probes))
    pStop("probe",sQuote("probes")," must be furnished.")
  if (!is.list(probes)) probes <- list(probes)
  if (!all(sapply(probes,is.function)))
    pStop("probe",sQuote("probes")," must be a function or a list of functions.")
  if (!all(sapply(probes,function(f)length(formals(f))==1)))
    pStop("probe","each probe must be a function of a single argument.")

  nsim <- as.integer(nsim)
  if (length(nsim) < 1)
    pStop("probe",sQuote("nsim")," must be specified.")
  if (length(nsim) > 1 || !is.finite(nsim) || nsim <= 0)
    pStop("probe","number of simulations, ",sQuote("nsim"),
      ", must be a single positive integer.")

  seed <- as.integer(seed)
  gnsi <- as.logical(.getnativesymbolinfo)

  params <- coef(object)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  ## apply probes to data
  datval <- tryCatch(
    .Call(apply_probe_data,object,probes),
    error = function (e)
      pStop("probe","applying probes to actual data: ",conditionMessage(e))
  )

  nprobes <- length(datval)

  if (nprobes >= nsim)
    pStop("probe",sQuote("nsim")," (=",nsim,"), should be (much) larger than the ",
      "number of probes (=",nprobes,").")

  ## apply probes to model simulations
  simval <- tryCatch(
    freeze(
      .Call(apply_probe_sim,object=object,nsim=nsim,params=params,probes=probes,
        datval=datval,gnsi=gnsi),
      seed=seed
    ),
    error = function (e)
      pStop("probe","applying probes to simulated data: ",conditionMessage(e))
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
    error = function (e)
      pStop("probe","in synthetic likelihood computation: ",conditionMessage(e))
  )

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

##' @name summary-probed_pomp
##' @aliases summary,probed_pomp-method
##' @rdname summary
setMethod(
  "summary",
  signature=signature(object="probed_pomp"),
  definition=function (object) {
    list(
      coef=coef(object),
      nsim=nrow(object@simvals),
      quantiles=object@quantiles,
      pvals=object@pvals,
      synth.loglik=object@synth.loglik
    )
  }
)

##' @name logLik-probed_pomp
##' @aliases logLik,probed_pomp-method
##' @rdname loglik
##'
##' @return
##' When \code{object} is of \sQuote{probed_pomp} class (i.e., the result of a \code{probe} computation), \code{logLik} retrieves the \dQuote{synthetic likelihood} (see \code{\link{probe}}).
##'
setMethod(
  "logLik",
  signature=signature(object="probed_pomp"),
  definition=function(object)object@synth.loglik
)
