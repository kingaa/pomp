##' Gompertz model with log-normal observations.
##'
##' \code{gompertz} is a \sQuote{pomp} object encoding a stochastic Gompertz
##' population model with log-normal measurement error.
##'
##' The state process is \eqn{X_{t+1} = K^{1-S} X_{t}^S
##' \epsilon_{t}}{X[t+1]=K^(1-S) X[t]^S eps[t]}, where \eqn{S=e^{-r}}{S=e^{-r}}
##' and the \eqn{\epsilon_t}{eps[t]} are i.i.d. lognormal random deviates with
##' variance \eqn{\sigma^2}{sigma^2}.  The observed variables \eqn{Y_t} are
##' distributed as
##' \eqn{\mathrm{lognormal}(\log{X_t},\tau)}{lognormal(log(X[t]),tau)}.
##' Parameters include the per-capita growth rate \eqn{r}, the carrying
##' capacity \eqn{K}, the process noise s.d. \eqn{\sigma}{sigma}, the
##' measurement error s.d. \eqn{\tau}{tau}, and the initial condition
##' \eqn{X_0}{X[0]}.  The \sQuote{pomp} object includes parameter
##' transformations that log-transform the parameters for estimation purposes.
##'
##' @name gompertz
##' @docType data
##' @keywords models datasets
##' @family pomp examples
##' @examples
##'
##' pompExample(gompertz)
##' plot(gompertz)
##' coef(gompertz)
##'
NULL
