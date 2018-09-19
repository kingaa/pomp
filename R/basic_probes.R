##' Useful probes for partially-observed Markov processes
##'
##' Several simple and configurable probes are provided with in the package.
##' These can be used directly and as templates for custom probes.
##'
##' @name basic_probes
##' @rdname basic_probes
##' @importFrom stats median spec.pgram kernel quantile sd var
##'
##' @param var,vars character; the name(s) of the observed variable(s).
##'
##' @param trim the fraction of observations to be trimmed (see \code{\link{mean}}).
##'
##' @param transform transformation to be applied to the data before the probe is computed.
##'
##' @param na.rm if \code{TRUE}, remove all NA observations prior to computing the probe.
##'
##' @param kernel.width width of modified Daniell smoothing kernel to be used
##' in power-spectrum computation: see \code{\link{kernel}}.
##'
##' @param lags In \code{probe.ccf}, a vector of lags between time series.
##' Positive lags correspond to \code{x} advanced relative to \code{y};
##' negative lags, to the reverse.
##'
##' In \code{probe.nlar}, a vector of lags present in the nonlinear
##' autoregressive model that will be fit to the actual and simulated data.
##' See Details, below, for a precise description.
##'
##' @param powers the powers of each term (corresponding to \code{lags}) in the
##' the nonlinear autoregressive model that will be fit to the actual and
##' simulated data.  See Details, below, for a precise description.
##'
##' @param type Compute autocorrelation or autocovariance?
##'
##' @param ... additional arguments passed to the underlying algorithms.
##'
##' @return
##' A call to any one of these functions returns a probe function,
##' suitable for use in \code{\link{probe}} or \code{\link{probe.match}}.  That
##' is, the function returned by each of these takes a data array (such as
##' comes from a call to \code{\link{obs}}) as input and returns a single
##' numerical value.
##'
##' @author Daniel C. Reuman, Aaron A. King
##' @family summary statistics methods
##' @concept diagnostics
##'
##' @references
##' B. E. Kendall, C. J. Briggs, W. M. Murdoch, P. Turchin, S. P.
##' Ellner, E. McCauley, R. M. Nisbet, S. N. Wood Why do populations cycle? A
##' synthesis of statistical and mechanistic modeling approaches, Ecology,
##' 80:1789--1805, 1999.
##'
##' S. N. Wood Statistical inference for noisy nonlinear ecological dynamic
##' systems, Nature, 466: 1102--1104, 2010.
##'
NULL

##'@rdname basic_probes
##' @export
probe.mean <- function (var, trim = 0, transform = identity, na.rm = TRUE) {
  if (length(var)>1)
    pStop_(sQuote("probe.mean")," is a univariate probe.")
  transform <- match.fun(transform)
  function(y) mean(x=transform(y[var,]),trim=trim,na.rm=na.rm)
}

##'@rdname basic_probes
##' @export
probe.median <- function (var, na.rm = TRUE) {
  if (length(var)>1) pStop_(sQuote("probe.median")," is a univariate probe.")
  function(y) median(x=as.numeric(y[var,]),na.rm=na.rm)
}

##'@rdname basic_probes
##' @export
probe.var <- function (var, transform = identity, na.rm = TRUE) {
  if (length(var)>1) pStop_(sQuote("probe.var")," is a univariate probe.")
  transform <- match.fun(transform)
  function(y) var(x=transform(y[var,]),na.rm=na.rm)
}

##'@rdname basic_probes
##' @export
probe.sd <- function (var, transform = identity, na.rm = TRUE) {
  if (length(var)>1) pStop_(sQuote("probe.sd")," is a univariate probe.")
  transform <- match.fun(transform)
  function(y) sd(x=transform(y[var,]),na.rm=na.rm)
}

##'@rdname basic_probes
##' @export
probe.period <- function (var, kernel.width, transform = identity) {
  if (length(var)>1) pStop_(sQuote("probe.period")," is a univariate probe.")
  transform <- match.fun(transform)
  function (y) {
    zz <- spec.pgram(
      x=transform(y[var,]),
      kernel=kernel("modified.daniell",m=kernel.width,r=NA),
      taper=0,
      fast=FALSE,
      pad=0,
      detrend=FALSE,
      plot=FALSE
    )
    1/zz$freq[which.max(zz$spec)]
  }
}

##' @rdname basic_probes
##' @param probs the quantile or quantiles to compute: see \code{\link{quantile}}.
##' @export
##'
probe.quantile <- function (var, probs, ...) {
  if (length(var)>1) pStop_(sQuote("probe.quantile")," is a univariate probe.")
  function (y) quantile(y[var,],probs=probs, ...)
}

##' @rdname basic_probes
##' @export
probe.acf <- function (var, lags, type = c("covariance", "correlation"),
  transform = identity) {
  ep <- "probe.acf"
  type <- match.arg(type)
  corr <- type=="correlation"
  transform <- match.fun(transform)
  if (corr && any(lags<=0)) pStop(ep,"lags must be positive integers.")
  lags <- as.integer(lags)
  function (y) tryCatch(
    .Call(P_probe_acf,x=transform(y[var,,drop=FALSE]),lags=lags,corr=corr),
    error = function (e) pStop(ep,conditionMessage(e))
  )
}

##' @rdname basic_probes
##' @export
probe.ccf <- function (vars, lags, type = c("covariance", "correlation"),
  transform = identity) {
  ep <- "probe.ccf"
  type <- match.arg(type)
  corr <- type=="correlation"
  transform <- match.fun(transform)
  if (length(vars)!=2)
    pStop(ep,sQuote("vars")," must name two variables.")
  lags <- as.integer(lags)
  function (y) tryCatch(
    .Call(P_probe_ccf,x=transform(y[vars[1L],,drop=TRUE]),y=transform(y[vars[2L],,drop=TRUE]),lags=lags,corr=corr),
    error = function (e) pStop(ep,conditionMessage(e))
  )
}

##' @rdname basic_probes
##' @param ref empirical reference distribution.  Simulated data will be
##' regressed against the values of \code{ref}, sorted and, optionally,
##' differenced.  The resulting regression coefficients capture information
##' about the shape of the marginal distribution.  A good choice for \code{ref}
##' is the data itself.
##'
##' @param order order of polynomial regression.
##'
##' @param diff order of differencing to perform.
##'
##' @export
probe.marginal <- function (var, ref, order = 3, diff = 1,
  transform = identity) {
  if (length(var)>1) pStop_(sQuote("probe.marginal")," is a univariate probe.")
  transform <- match.fun(transform)
  setup <- .Call(P_probe_marginal_setup,transform(ref),order,diff)
  function (y) tryCatch(
    .Call(P_probe_marginal_solve,x=transform(y[var,,drop=TRUE]),setup=setup,diff=diff),
    error = function (e) pStop("probe.marginal",conditionMessage(e))
  )
}

##' @rdname basic_probes
##' @export
##'
probe.nlar <- function (var, lags, powers, transform = identity) {
  ep <- "probe.nlar"
  if (length(var)>1) pStop_(sQuote(ep)," is a univariate probe.")
  transform <- match.fun(transform)
  if (missing(lags) || missing(powers))
    pStop(ep,sQuote("lags")," and ",sQuote("powers")," are required arguments.")
  lags <- as.integer(lags)
  powers <- as.integer(powers)
  if (any(lags<1)||any(powers<1))
    pStop(ep,sQuote("lags")," and ",sQuote("powers")," must be positive integers.")
  if (length(lags)<length(powers)) {
    if (length(lags)>1)
      pStop(ep,sQuote("lags")," must match ",sQuote("powers")," in length, or have length 1.")
    lags <- rep(lags,length(powers))
  } else if (length(lags)>length(powers)) {
    if (length(powers)>1)
      pStop(ep,sQuote("powers")," must match ",sQuote("lags")," in length, or have length 1.")
    powers <- rep(powers,length(lags))
  }
  lags <- as.integer(lags)
  powers <- as.integer(powers)
  function (y) tryCatch(
    .Call(P_probe_nlar,x=transform(y[var,,drop=TRUE]),lags=lags,powers=powers),
    error = function (e) pStop(ep,conditionMessage(e))
  )
}
