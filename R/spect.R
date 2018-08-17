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
##' \code{spect.match} tries to match the power spectrum of the model to that
##' of the data.  It calls an optimizer to adjust model parameters to minimize
##' the discrepancy between simulated and actual data.
##'
##' A call to \code{spect} results in the estimation of the power spectrum for
##' the (transformed, detrended) data and \code{nsim} model simulations.  The
##' results of these computations are stored in an object of class
##' \sQuote{spectd_pomp}.
##'
##' A call to \code{spect.match} results in an attempt to optimize the
##' agreement between model and data spectrum over the parameters named in
##' \code{est}.  The results, including coefficients of the fitted model and
##' power spectra of fitted model and data, are stored in an object of class
##' \sQuote{spect_matched_pomp}.
##'
##' @name Power spectrum
##' @rdname spect
##' @include simulate_pomp.R
##' @aliases spect spect,missing-method spect,ANY-method
##'
##' @param object An object of class \sQuote{pomp}.
##' @param params optional named numeric vector of model parameters.  By
##' default, \code{params=coef(object)}.
##' @param vars optional; names of observed variables for which the power
##' spectrum will be computed.  This must be a subset of
##' \code{rownames(obs(object))}.  By default, the spectrum will be computed
##' for all observables.
##' @param kernel.width width parameter for the smoothing kernel used for
##' calculating the estimate of the spectrum.
##' @param nsim number of model simulations to be computed.
##' @param seed optional; if non-\code{NULL}, the random number generator will
##' be initialized with this seed for simulations.  See
##' \code{\link[=simulate-pomp]{simulate}}.
##' @param transform.data function; this transformation will be applied to the
##' observables prior to estimation of the spectrum, and prior to any
##' detrending.
##' @param detrend de-trending operation to perform.  Options include no
##' detrending, and subtraction of constant, linear, and quadratic trends from
##' the data.  Detrending is applied to each data series and to each model
##' simulation independently.
##' @param weights optional.  The mismatch between model and data is measured
##' by a weighted average of mismatch at each frequency.  By default, all
##' frequencies are weighted equally.  \code{weights} can be specified either
##' as a vector (which must have length equal to the number of frequencies) or
##' as a function of frequency.  If the latter, \code{weights(freq)} must
##' return a nonnegative weight for each frequency.
##' @param start named numeric vector; the initial guess of parameters.
##' @param est character vector; the names of parameters to be estimated.
##' @param method Optimization method.  Choices are
##' \code{\link[subplex]{subplex}}, \code{\link{sannbox}}, and any of the
##' methods used by \code{\link{optim}}.
##' @param verbose logical; print diagnostic messages?
##' @param fail.value optional scalar; if non-\code{NA}, this value is
##' substituted for non-finite values of the objective function.
##' @param \dots Additional arguments.  In the case of \code{spect}, these are
##' currently ignored.  In the case of \code{spect.match}, these are passed to
##' \code{optim} or \code{subplex} in the \code{control} list.
##'
##' @return
##' \code{spect} returns an object of class \sQuote{spectd_pomp}.
##'
##' \code{spect.match} returns an object of class \sQuote{spect_matched_pomp},
##' which is derived from class \sQuote{spectd_pomp} and therefore has all the
##' slots of that class.  In addition, \sQuote{spect_matched_pomp} objects have
##' the following slots: \describe{ \item{est, weights, fail.value}{values of
##' the corresponding arguments in the call to \code{spect.match}.}
##' \item{evals}{ number of function and gradient evaluations by the optimizer.
##' See \code{\link{optim}}.  } \item{value}{Value of the objective function.}
##' \item{convergence, msg}{ Convergence code and message from the optimizer.
##' See \code{\link{optim}}.  } }
##'
##' @author Daniel C. Reuman, Cai GoGwilt, Aaron A. King
##'
##' @seealso \code{\link{simulate}}, \code{\link{probe}}
##'
##' @references
##' D.C. Reuman, R.A. Desharnais, R.F. Costantino, O. Ahmad, J.E.
##' Cohen (2006) Power spectra reveal the influence of stochasticity on
##' nonlinear population dynamics.  \emph{Proceedings of the National Academy
##' of Sciences} \bold{103}, 18860-18865.
##'
##' D.C. Reuman, R.F. Costantino, R.A. Desharnais, J.E. Cohen (2008) Color of
##' environmental noise affects the nonlinear dynamics of cycling,
##' stage-structured populations.  \emph{Ecology Letters}, \bold{11}, 820-830.
NULL

## power spectrum
## Authors:
## Cai GoGwilt, Daniel Reuman, Aaron A. King

setClass(
  "spectd_pomp",
  contains="pomp",
  slots=c(
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
  function (object, ...)
    standardGeneric("spect")
)

##' @name spect-pomp
##' @aliases spect,pomp-method
##' @rdname spect
setMethod(
  "spect",
  signature(object="pomp"),
  function (object, params, vars, kernel.width, nsim, seed = NULL,
    transform.data = identity, detrend = c("none","mean","linear","quadratic"),
    ...) {

    detrend <- match.arg(detrend)
    transform.data <- match.fun(transform.data)

    spect.internal(
      object,
      params=params,
      vars=vars,
      kernel.width=kernel.width,
      nsim=nsim,
      seed=seed,
      transform.data=transform.data,
      detrend=detrend,
      ...
    )

  }
)

##' @name spect-spectd_pomp
##' @aliases spect,spectd_pomp-method
##' @rdname spect
setMethod(
  "spect",
  signature=signature(object="spectd_pomp"),
  definition=function (object, params, vars, kernel.width,
    nsim, seed = NULL, transform.data, detrend, ...) {

    if (missing(params)) params <- coef(object)
    if (missing(vars)) vars <- colnames(object@datspec)
    if (missing(kernel.width)) kernel.width <- object@kernel.width
    if (missing(nsim)) nsim <- nrow(object@simspec)
    if (missing(transform.data)) transform.data <- object@transform.data
    if (missing(detrend)) detrend <- object@detrend

    spect(
      as(object,"pomp"),
      params=params,
      vars=vars,
      kernel.width=kernel.width,
      nsim=nsim,
      seed=seed,
      transform.data=transform.data,
      detrend=detrend,
      ...
    )

  }
)

setMethod(
  "spect",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("spect"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "spect",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("spect")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

spect.internal <- function (object, params, vars, kernel.width, nsim,
  seed = NULL, transform.data, detrend, ...) {

  ep <- paste0("in ",sQuote("spect"),": ")

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)

  if (missing(vars)) vars <- rownames(object@data)

  if (missing(kernel.width) || length(kernel.width) > 1 ||
      !is.numeric(kernel.width) ||
      !is.finite(kernel.width) || kernel.width < 0)
    stop(ep,sQuote("kernel.width"),
      " must be specified as a single positive integer.",call.=FALSE)

  if (missing(nsim) || length(nsim) > 1 || !is.numeric(nsim)||
      !is.finite(nsim) || (nsim<1))
    stop(ep,sQuote("nsim")," must be specified as a positive integer.",
      call.=FALSE)

  ker <- reuman.kernel(kernel.width)

  object <- pomp(object,...)

  pompLoad(object)
  on.exit(pompUnload(object))

  ds <- compute.spect.data(
    object,
    vars=vars,
    transform.data=transform.data,
    detrend=detrend,
    ker=ker
  )
  freq <- ds$freq
  datspec <- ds$spec

  simspec <- compute.spect.sim(
    object,
    params=params,
    vars=vars,
    nsim=nsim,
    seed=seed,
    transform.data=transform.data,
    detrend=detrend,
    ker=ker
  )

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
      function(k)sum((simspec[k,,j]-mean.simspec[,j])^2),
      numeric(1)
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

## detrends in one of several ways, according to type.
## tseries is a numeric vector,
pomp.detrend <- function (tseries, type) {
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
reuman.kernel <- function (kernel.width) {
  ker <- kernel("modified.daniell",m=kernel.width)
  x <- seq.int(from=0,to=kernel.width,by=1)/kernel.width
  ker[[1L]] <- (15/(16*2*pi))*((x-1)^2)*((x+1)^2)
  ker[[1L]] <- ker[[1L]]/(2*sum(ker[[1L]][-1])+ker[[1L]][1L])
  attr(ker,"name") <- NULL
  ker
}

compute.spect.data <- function (object, vars, transform.data, detrend, ker) {

  ep <- paste0("in ",sQuote("spect"),": ")

  dat <- obs(object,vars)
  if (any(!is.finite(dat)))
    stop(ep,"missing or infinite values in the data.",call.=FALSE)

  dt <- diff(time(object,t0=FALSE))
  base.freq <- 1/mean(dt)
  dt.tol <- 0.025
  if (max(dt)-min(dt)>dt.tol*mean(dt))
    stop(ep,sQuote("spect")," assumes evenly spaced times.",call.=FALSE)

  for (j in seq_along(vars)) {
    sp <- spec.pgram(
      pomp.detrend(transform.data(dat[j,]),type=detrend),spans=ker,taper=0,
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

compute.spect.sim <- function (object, params, vars, nsim, seed,
  transform.data, detrend, ker) {

  ep <- paste0("in ",sQuote("spect"),": ")

  sims <- tryCatch(
    simulate(object,nsim=nsim,seed=seed,params=params,states=FALSE,obs=TRUE),
    error = function (e) {
      stop(ep,conditionMessage(e),call.=FALSE)
    }
  )

  sims <- sims[vars,,,drop=FALSE]

  if (any(!is.finite(sims)))
    stop(ep,"missing or infinite values in simulated data.",call.=FALSE)

  nobs <- length(vars)
  for (j in seq_len(nobs)) {
    for (k in seq_len(nsim)) {
      sp <- spec.pgram(pomp.detrend(transform.data(sims[j,k,]),type=detrend),
        spans=ker,taper=0,pad=0,fast=FALSE,detrend=FALSE,plot=FALSE)
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
