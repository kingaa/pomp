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
  var1 <- vars[1]
  if (length(vars)>1)
    var2 <- vars[2]
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
  var1 <- vars[1]
  if (length(vars)>1)
    var2 <- vars[2]
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

probe.marginal <- function (var, ref, order = 3, diff = 1, transform = identity) {
  if (length(var)>1) stop(sQuote("probe.marginal")," is a univariate probe")
  transform <- match.fun(transform)
  setup <- .Call(probe_marginal_setup,transform(ref),order,diff)
  function (y) .Call(
                     probe_marginal_solve,
                     x=transform(y[var,]),
                     setup=setup,
                     diff=diff
                     )
}

probe.acf <- function (var, lag.max, type = c("covariance", "correlation"), transform = identity) {
  type <- match.arg(type)
  transform <- match.fun(transform)
  corr <- type=="correlation"
  function (y) .Call(
                     probe_acf,
                     x=transform(y[var,,drop=FALSE]),
                     lag_max=lag.max,
                     corr=corr
                     )
}
