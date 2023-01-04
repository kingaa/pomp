##' The log-mean-exp trick
##'
##' \code{logmeanexp} computes \deqn{\log\frac{1}{N}\sum_{n=1}^N\!e^x_i,}{log
##' mean exp(x_i),} avoiding over- and under-flow in doing so.  It can
##' optionally return an estimate of the standard error in this quantity.
##'
##' When \code{se = TRUE}, \code{logmeanexp} uses a jackknife estimate of the
##' variance in \eqn{log(x)}.
##'
##' @importFrom stats sd
##' @param x numeric
##' @param se logical; give approximate standard error?
##' @return
##' \code{log(mean(exp(x)))} computed so as to avoid over- or
##' underflow.  If \code{se = FALSE}, the approximate standard error is
##' returned as well.
##' @author Aaron A. King
##' @example examples/logmeanexp.R
##' @export
logmeanexp <- function (x, se = FALSE) {
  se <- isTRUE(se)
  mx <- max(x)
  mean <- mx+log(mean(exp(x-mx)))
  if (se) {
    n <- length(x)
    jk <- vapply(
      seq_len(n),
      \(k) logmeanexp(x[-k]),
      numeric(1L)
    )
    c(mean,se=(n-1)*sd(jk)/sqrt(n))
  } else {
    mean
  }
}
