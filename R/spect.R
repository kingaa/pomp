##' Power spectrum
##'
##' Power spectrum computation and spectrum-matching for partially-observed
##' Markov processes.
##'
##' \code{spect} estimates the power spectrum of time series data and model
##' simulations and compares the results.  It can be used to diagnose goodness
##' of fit and/or as the basis for frequency-domain parameter estimation
##' (\code{spect.match}).
##'
##' A call to \code{spect} results in the estimation of the power spectrum for
##' the (transformed, detrended) data and \code{nsim} model simulations.  The
##' results of these computations are stored in an object of class
##' \sQuote{spectd_pomp}.
##'
##' @name spect
##' @docType methods
##' @rdname spect
##' @aliases spect,missing-method spect,ANY-method
##' @family summary statistic-based methods
##' @family elementary algorithms
##'
##' @inheritSection pomp Note for Windows users
##'
##' @include simulate.R pomp.R
##' @importFrom stats spec.pgram kernel .lm.fit
##'
##' @param vars optional; names of observed variables for which the power spectrum will be computed.
##' By default, the spectrum will be computed for all observables.
##' @param kernel.width width parameter for the smoothing kernel used for
##' calculating the estimate of the spectrum.
##' @param nsim number of model simulations to be computed.
##' @param seed optional; if non-\code{NULL}, the random number generator will
##' be initialized with this seed for simulations.
##' See \code{\link{simulate}}.
##' @param transform.data function; this transformation will be applied to the
##' observables prior to estimation of the spectrum, and prior to any
##' detrending.
##' @param detrend de-trending operation to perform.  Options include no
##' detrending, and subtraction of constant, linear, and quadratic trends from
##' the data.  Detrending is applied to each data series and to each model
##' simulation independently.
##' @inheritParams probe
##' @inheritParams pomp
##'
##' @return
##' An object of class \sQuote{spectd_pomp}, which contains the model, the data, and the results of the \code{spect} computation.
##' The following methods are available:
##' \describe{
##' \item{plot}{produces some diagnostic plots}
##' \item{summary}{displays a summary}
##' \item{logLik}{gives a measure of the agreement of the power spectra}
##' }
##'
##' @author Daniel C. Reuman, Cai GoGwilt, Aaron A. King
##'
##' @references
##'
##' \Reuman2006
##'
##' \Reuman2008
##'
NULL

setClass(
  "spectd_pomp",
  contains="pomp",
  slots=c(
    nsim="integer",
    seed="integer",
    kernel.width="numeric",
    transform.data="function",
    vars="character",
    freq="numeric",
    datspec="array",
    simspec="array",
    pvals="numeric",
    detrend="character"
  )
)

setGeneric(
  "spect",
  function (data, ...)
    standardGeneric("spect")
)

setMethod(
  "spect",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("spect","data")
  }
)

setMethod(
  "spect",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("spect",data)
  }
)

##' @rdname spect
##' @export
setMethod(
  "spect",
  signature(data="data.frame"),
  function (
    data,
    ...,
    vars, kernel.width, nsim, seed = NULL,
    transform.data = identity,
    detrend = c("none","mean","linear","quadratic"),
    params, rinit, rprocess, rmeasure,
    verbose = getOption("verbose", FALSE)
  ) {

    tryCatch(
      spect_internal(
        data,
        ...,
        vars=vars,
        kernel.width=kernel.width,
        nsim=nsim,
        seed=seed,
        transform.data=match.fun(transform.data),
        detrend=match.arg(detrend),
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        rmeasure=rmeasure,
        verbose=verbose
      ),
      error = function (e) pStop(who="spect",conditionMessage(e))
    )

  }
)
##' @rdname spect
##' @export
setMethod(
  "spect",
  signature(data="pomp"),
  function (
    data,
    ...,
    vars, kernel.width, nsim, seed = NULL,
    transform.data = identity,
    detrend = c("none","mean","linear","quadratic"),
    verbose = getOption("verbose", FALSE)
  ) {

    tryCatch(
      spect_internal(
        data,
        ...,
        vars=vars,
        kernel.width=kernel.width,
        nsim=nsim,
        seed=seed,
        transform.data=match.fun(transform.data),
        detrend=match.arg(detrend),
        verbose=verbose
      ),
      error = function (e) pStop(who="spect",conditionMessage(e))
    )

  }
)

##' @rdname spect
##' @export
setMethod(
  "spect",
  signature=signature(data="spectd_pomp"),
  definition=function (
    data,
    ...,
    vars, kernel.width, nsim, seed = NULL,
    transform.data,
    detrend,
    verbose = getOption("verbose", FALSE)
  ) {

    if (missing(vars)) vars <- colnames(data@datspec)
    if (missing(kernel.width)) kernel.width <- data@kernel.width
    if (missing(nsim)) nsim <- data@nsim
    if (missing(transform.data)) transform.data <- data@transform.data
    if (missing(detrend)) detrend <- data@detrend

    spect(
      as(data,"pomp"),
      ...,
      vars=vars,
      kernel.width=kernel.width,
      nsim=nsim,
      seed=seed,
      transform.data=transform.data,
      detrend=detrend,
      verbose=verbose
    )

  }
)

spect_internal <- function (
  object,
  ...,
  vars, kernel.width, nsim, seed = NULL,
  transform.data, detrend,
  .gnsi = TRUE, verbose
) {

  verbose <- as.logical(verbose)

  object <- pomp(object,...,verbose=verbose)

  if (undefined(object@rprocess) || undefined(object@rmeasure))
    pStop_(paste(sQuote(c("rprocess","rmeasure")),collapse=", ")," are needed basic components.")

  if (missing(vars)) vars <- rownames(object@data)

  if (missing(kernel.width) || length(kernel.width) > 1 ||
        !is.numeric(kernel.width) ||
        !is.finite(kernel.width) || kernel.width < 0)
    pStop_(sQuote("kernel.width")," must be a positive integer.")

  if (missing(nsim) || length(nsim) > 1 || !is.numeric(nsim)||
        !is.finite(nsim) || (nsim<1))
    pStop_(sQuote("nsim")," must be a positive integer.")

  nsim <- as.integer(nsim)
  seed <- as.integer(seed)

  ker <- reuman_kernel(kernel.width)

  params <- coef(object)

  pompLoad(object)
  on.exit(pompUnload(object))

  ds <- compute_spect_data(
    object,
    vars=vars,
    transform.data=transform.data,
    detrend=detrend,
    ker=ker
  )
  freq <- ds$freq
  datspec <- ds$spec

  simspec <- compute_spect_sim(
    object,
    params=params,
    vars=vars,
    nsim=nsim,
    seed=seed,
    transform.data=transform.data,
    detrend=detrend,
    ker=ker,
    .gnsi=.gnsi
  )
  .gnsi <- FALSE

  pvals <- numeric(length(vars)+1)
  names(pvals) <- c(vars,"all")
  mean.simspec <- colMeans(simspec) # mean spectrum of simulations
  totdatdist <- 0
  totsimdist <- 0
  for (j in seq_along(vars)) {
    ## L-2 distance between data and mean simulated spectrum
    datdist <- sum((datspec[,j]-mean.simspec[,j])^2)
    ## L-2 distance betw. each sim. and mean simulated spectrum
    simdist <- vapply(
      seq_len(nsim),
      \(k)sum((simspec[k,,j]-mean.simspec[,j])^2),
      numeric(1L)
    )
    pvals[j] <- (nsim+1-sum(simdist<datdist))/(nsim+1)
    totdatdist <- totdatdist+datdist
    totsimdist <- totsimdist+simdist
  }
  pvals[length(vars)+1] <- (nsim+1-sum(totsimdist<totdatdist))/(nsim+1)

  coef(object) <- params

  new(
    "spectd_pomp",
    object,
    nsim=nsim,
    seed=seed,
    kernel.width=kernel.width,
    transform.data=transform.data,
    vars=vars,
    detrend=detrend,
    freq=freq,
    datspec=datspec,
    simspec=simspec,
    pvals=pvals
  )
}

compute_spect_data <- function (object, vars, transform.data, detrend, ker) {

  dat <- obs(object,vars)
  if (any(!is.finite(dat)))
    pStop_("missing or infinite values in the data.")

  dt <- diff(time(object,t0=FALSE))
  base.freq <- 1/mean(dt)
  dt.tol <- 0.025
  if (max(dt)-min(dt)>dt.tol*mean(dt))
    pStop_(sQuote("spect")," assumes evenly spaced times.")

  for (j in seq_along(vars)) {
    sp <- spec.pgram(
      pomp_detrend(transform.data(dat[j,]),type=detrend),spans=ker,taper=0,
      pad=0,fast=FALSE,detrend=FALSE,plot=FALSE
    )
    if (j==1) {
      freq <- base.freq*sp$freq
      datspec <- array(
        dim=c(length(freq),nrow(dat)),
        dimnames=list(NULL,vars)
      )
    }
    datspec[,j] <- log10(sp$spec)
  }
  list(freq=freq,spec=datspec)
}

compute_spect_sim <- function (object, params, vars, nsim, seed,
  transform.data, detrend, ker, .gnsi = TRUE) {

  sims <- tryCatch(
  {
    s <- freeze(
      .Call(P_do_simulate,object,params,nsim,0L,gnsi=.gnsi),
      seed=seed
    )
    s$obs[vars,,,drop=FALSE]
  },
  error = function (e) pStop_("in simulation: ",
    conditionMessage(e))
  )

  if (any(!is.finite(sims)))
    pStop_("missing or infinite values in simulated data.")

  nobs <- length(vars)
  for (j in seq_len(nobs)) {
    for (k in seq_len(nsim)) {
      sp <- tryCatch(
        spec.pgram(pomp_detrend(transform.data(sims[j,k,]),type=detrend),
          spans=ker,taper=0,pad=0,fast=FALSE,detrend=FALSE,plot=FALSE),
        error = function (e) pStop(who="spec.pgram",conditionMessage(e))
      )
      if ((j==1)&&(k==1)) {
        simspec <- array(
          dim=c(nsim,length(sp$freq),nobs),
          dimnames=list(NULL,NULL,vars)
        )
      }
      simspec[k,,j] <- log10(sp$spec)
    }
  }
  simspec
}

## detrends in one of several ways, according to type.
## tseries is a numeric vector,
pomp_detrend <- function (tseries, type) {
  switch(
    type,
    mean=tseries-mean(tseries),
    linear={
      m <- cbind(1,seq_along(tseries))
      .lm.fit(m,tseries)$residuals
    },
    quadratic={
      x <- seq_along(tseries)
      m <- cbind(1,x,x*x)
      .lm.fit(m,tseries)$residuals
    },
    tseries
  )
}

## The default smoothing kernel for the R spec.pgram function is weird.
## This function creates a better one.
reuman_kernel <- function (kernel.width) {
  ker <- kernel("modified.daniell",m=kernel.width,r=NA)
  x <- seq.int(from=0,to=kernel.width,by=1)/kernel.width
  ker[[1L]] <- (15/(16*2*pi))*((x-1)^2)*((x+1)^2)
  ker[[1L]] <- ker[[1L]]/(2*sum(ker[[1L]][-1])+ker[[1L]][1L])
  attr(ker,"name") <- NULL
  ker
}

##' @rdname summary
##' @export
setMethod(
  "summary",
  signature=signature(object="spectd_pomp"),
  definition=function (object, ...) {
    list(
      coef=coef(object),
      nsim=nrow(object@simspec),
      pvals=object@pvals
    )
  }
)
