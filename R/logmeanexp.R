##' The log-mean-exp trick
##'
##' \code{logmeanexp} computes
##' \deqn{\log\frac{1}{n}\sum_{i=1}^n\!e^{x_i},}{log(mean(exp(x))),}
##' avoiding over- and under-flow in doing so.
##' It can optionally return an estimate of the standard error in this quantity.
##'
##' When \code{se = TRUE}, \code{logmeanexp} uses a jackknife estimate of the variance in \eqn{log(x)}.
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
  x <- as.numeric(x)
  est <- .Call(P_logmeanexp,x,-1L)
  if (se || ess) {
    if (se) {
      n <- length(x)
      jk <- vapply(
        seq_len(n),
        \(k) .Call(P_logmeanexp,x,k),
        double(1L)
      )
      xse <- (n-1)*sd(jk)/sqrt(n)
    }
    if (ess) {
      w <- exp(x-max(x))
      xss <- sum(w)^2/sum(w^2)
      if (se) {
        c(est=est,se=xse,ess=xss)
      } else {
        c(est=est,ess=xss)
      }
    } else {
      c(est=est,se=xse)
    }
  } else {
    est
  }
}
