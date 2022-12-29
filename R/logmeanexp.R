##' The log-mean-exp trick
##'
##' \code{logmeanexp} computes \deqn{\log\frac{1}{N}\sum_{n=1}^N\!e^x_i,}{log
##' mean exp(x_i),} avoiding over- and under-flow in doing so.  It can
##' optionally return an estimate of the standard error in this quantity.
##'
##' When \code{se = TRUE}, \code{logmeanexp} uses a jackknife estimate of the
##' variance in \eqn{log(x)}.
##' 
##' When \code{ess = TRUE}, \code{logmeanexp} returns an estimate of the effective sample size.
##'
##' @importFrom stats sd
##' @param x numeric
##' @param se logical; give approximate standard error?
##' @param ess logical; give effective sample size?
##' @return
##' \code{log(mean(exp(x)))} computed so as to avoid over- or underflow.
##' If \code{se = TRUE}, the approximate standard error is returned as well.
##' If \code{ess = TRUE}, the effective sample size is returned also.
##' @author Aaron A. King
##' @example examples/logmeanexp.R
##' @export
logmeanexp <- function (x, se = FALSE, ess = FALSE) {
  se <- isTRUE(se)
  ess <- isTRUE(ess)
  mx <- max(x)
  w <- exp(x-mx)
  rv <- mx+log(mean(w))
  if (se) {
    n <- length(x)
    jk <- vapply(
      seq_len(n),
      \(k) logmeanexp(x[-k]),
      numeric(1L)
    )
    rv <- c(rv,se=(n-1)*sd(jk)/sqrt(n))
  }
  if (ess) {
    rv <- c(rv,ess=sum(w)^2/sum(w^2))
  }
  rv
}
