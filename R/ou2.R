##' Two-dimensional discrete-time Ornstein-Uhlenbeck process
##'
##' \code{ou2} is a \sQuote{pomp} object encoding a bivariate discrete-time
##' Ornstein-Uhlenbeck process.
##'
##' If the state process is \eqn{X(t) = (x_{1}(t),x_{2}(t))}, then \deqn{X(t+1)
##' = \alpha X(t) + \sigma \epsilon(t),} where \eqn{\alpha} and \eqn{\sigma}
##' are 2x2 matrices, \eqn{\sigma} is lower-triangular, and \eqn{\epsilon(t)}
##' is standard bivariate normal.  The observation process is \eqn{Y(t) =
##' (y_1(t),y_2(t))}, where \eqn{y_i(t) \sim \mathrm{normal}(x_i(t),\tau)}.
##' The functions \code{rprocess}, \code{dprocess}, \code{rmeasure},
##' \code{dmeasure}, and \code{skeleton} are implemented using compiled C code
##' for computational speed: see the source code for details.
##'
##' @name ou2
##' @rdname ou2
##' @docType data
##' @seealso \code{\link{pomp}}
##' @keywords models datasets
##' @examples
##'
##' pompExample(ou2)
##' plot(ou2)
##' coef(ou2)
##' x <- simulate(ou2)
##' plot(x)
##' pf <- pfilter(ou2,Np=1000)
##' logLik(pf)
##'
NULL
