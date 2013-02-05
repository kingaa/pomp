probe.mean <- function (var, trim = 0, transform = identity, na.rm = TRUE) {
  if (length(var)>1) stop(sQuote("probe.mean")," is a univariate probe")
  transform <- match.fun(transform)
  function(y) mean(x=transform(y[var,]),trim=trim,na.rm=na.rm)
}

probe.median <- function (var, na.rm = TRUE) {
  if (length(var)>1) stop(sQuote("probe.median")," is a univariate probe")
  function(y) median(x=as.numeric(y[var,]),na.rm=na.rm)
}

probe.var <- function (var, transform = identity, na.rm = TRUE) {
  if (length(var)>1) stop(sQuote("probe.var")," is a univariate probe")
  transform <- match.fun(transform)
  function(y) var(x=transform(y[var,]),na.rm=na.rm)
}

probe.sd <- function (var, transform = identity, na.rm = TRUE) {
  if (length(var)>1) stop(sQuote("probe.sd")," is a univariate probe")
  transform <- match.fun(transform)
  function(y) sd(x=transform(y[var,]),na.rm=na.rm)
}

probe.period <- function (var, kernel.width, transform = identity) {
  if (length(var)>1) stop(sQuote("probe.period")," is a univariate probe")
  transform <- match.fun(transform)
  function (y) {
    zz <- spec.pgram(
                     x=transform(y[var,]),
                     kernel=kernel("modified.daniell",m=kernel.width),
                     taper=0,
                     fast=FALSE,
                     pad=0,
                     detrend=FALSE,
                     plot=FALSE
                     )
    1/zz$freq[which.max(zz$spec)]
  }
}

probe.quantile <- function (var, prob, transform = identity) {
  if (length(var)>1) stop(sQuote("probe.quantile")," is a univariate probe")
  transform <- match.fun(transform)
  function (y) quantile(transform(y[var,]),probs=prob)
}

probe.cov <- function (
                       vars,
                       lag,
                       method = c("pearson", "kendall", "spearman"),
                       transform = identity
                       ) {
  method <- match.arg(method)
  lag <- as.integer(lag)
  transform <- match.fun(transform)
  var1 <- vars[1L]
  if (length(vars)>1)
    var2 <- vars[2L]
  else
    var2 <- var1
  function (y) {
    if (lag>=0) {
      val <- cov(
                 x=transform(y[var1,seq(from=1+lag,to=ncol(y),by=1)]),
                 y=transform(y[var2,seq(from=1,to=ncol(y)-lag,by=1)]),
                 method=method
                 )
    } else {
      val <- cov(
                 x=transform(y[var1,seq(from=1,to=ncol(y)+lag,by=1)]),
                 y=transform(y[var2,seq(from=-lag,to=ncol(y),by=1)]),
                 method=method
                 )
    }
    val
  }
}

probe.cor <- function (
                       vars,
                       lag,
                       method = c("pearson", "kendall", "spearman"),
                       transform = identity
                       ) {
  method <- match.arg(method)
  lag <- as.integer(lag)
  transform <- match.fun(transform)
  var1 <- vars[1L]
  if (length(vars)>1)
    var2 <- vars[2L]
  else
    var2 <- var1
  function (y) {
    if (lag>=0) {
      val <- cor(
                 x=transform(y[var1,seq(from=1+lag,to=ncol(y),by=1)]),
                 y=transform(y[var2,seq(from=1,to=ncol(y)-lag,by=1)]),
                 method=method
                 )
    } else {
      val <- cor(
                 x=transform(y[var1,seq(from=1,to=ncol(y)+lag,by=1)]),
                 y=transform(y[var2,seq(from=-lag,to=ncol(y),by=1)]),
                 method=method
                 )
    }
    val
  }
}

probe.acf <- function (var, lags, type = c("covariance", "correlation"), transform = identity) {
  type <- match.arg(type)
  corr <- type=="correlation"
  transform <- match.fun(transform)
  if (corr && any(lags==0)) {
    warning("useless zero lag discarded in ",sQuote("probe.acf"))
    lags <- lags[lags!=0]
  }
  lags <- as.integer(lags)
  function (y) .Call(
                     probe_acf,
                     x=transform(y[var,,drop=FALSE]),
                     lags=lags,
                     corr=corr
                     )
}

probe.ccf <- function (vars, lags, type = c("covariance", "correlation"), transform = identity) {
  type <- match.arg(type)
  corr <- type=="correlation"
  transform <- match.fun(transform)
  if (length(vars)!=2)
    stop(sQuote("vars")," must name two variables")
  lags <- as.integer(lags)
  function (y) .Call(
                     probe_ccf,
                     x=transform(y[vars[1L],,drop=TRUE]),
                     y=transform(y[vars[2L],,drop=TRUE]),
                     lags=lags,
                     corr=corr
                     )
}

probe.marginal <- function (var, ref, order = 3, diff = 1, transform = identity) {
  if (length(var)>1) stop(sQuote("probe.marginal")," is a univariate probe")
  transform <- match.fun(transform)
  setup <- .Call(probe_marginal_setup,transform(ref),order,diff)
  function (y) .Call(
                     probe_marginal_solve,
                     x=transform(y[var,,drop=TRUE]),
                     setup=setup,
                     diff=diff
                     )
}

probe.nlar <- function (var, lags, powers, transform = identity) {
  if (length(var)>1) stop(sQuote("probe.nlar")," is a univariate probe")
  transform <- match.fun(transform)
  if (any(lags<1)||any(powers<1))
    stop(sQuote("lags")," and ",sQuote("powers")," must be positive integers")
  if (length(lags)<length(powers)) {
    if (length(lags)>1) stop(sQuote("lags")," must match ",sQuote("powers")," in length, or have length 1")
    lags <- rep(lags,length(powers))
  } else if (length(lags)>length(powers)) {
    if (length(powers)>1) stop(sQuote("powers")," must match ",sQuote("lags")," in length, or have length 1")
    powers <- rep(powers,length(lags))
  }
  lags <- as.integer(lags)
  powers <- as.integer(powers)
  function (y) .Call(
                     probe_nlar,
                     x=transform(y[var,,drop=TRUE]),
                     lags=lags,
                     powers=powers
                     )
}
