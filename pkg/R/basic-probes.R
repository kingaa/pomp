probe.mean <- function (var, trim = 0, transform = identity, na.rm = TRUE) {
  function(y) mean(x=transform(y[var,]),trim=trim,na.rm=na.rm)
}

probe.var <- function (var, transform = identity, na.rm = TRUE) {
  function(y) var(x=transform(y[var,]),na.rm=na.rm)
}

probe.sd <- function (var, transform = identity, na.rm = TRUE) {
  function(y) sd(x=transform(y[var,]),na.rm=na.rm)
}

probe.period <- function (var, kernel.width, transform = identity) {
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
  function (y) quantile(transform(y[var,]),probs=prob)
}

probe.acf <- function (
                       var,
                       lag = 0,
                       type = c("correlation", "covariance", "partial"),
                       transform = identity,
                       ...
                       ) {
  args <- list(...)
  type <- match.arg(type)
  function (y) {
    zz <- do.call(acf,c(list(x=transform(y[var,]),lag.max=lag,plot=FALSE),args))
    if (type=="partial")
      val <- zz$acf[lag]
    else
      val <- zz$acf[lag+1]
    val
  }
}

probe.cov <- function (
                       vars,
                       lag,
                       method = c("pearson", "kendall", "spearman"),
                       transform = identity
                       ) {
  method <- match.arg(method)
  lag <- as.integer(lag)
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
