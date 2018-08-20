##' Ricker model with Poisson observations.
##'
##' \code{ricker} is a \sQuote{pomp} object encoding a stochastic Ricker model
##' with Poisson measurement error.
##'
##' The state process is \eqn{N_{t+1} = r N_{t} \exp(-c N_{t}+e_{t})}{N[t+1] =
##' r N[t] exp(-c N[t]+e[t])}, where the \eqn{e_t}{e[t]} are i.i.d. normal
##' random deviates with zero mean and variance \eqn{\sigma^2}{sigma^2}.  The
##' observed variables \eqn{y_t}{y[t]} are distributed as
##' \eqn{\mathrm{Poisson}(\phi N_t)}{Poisson(phi N[t])}.
##'
##' @name ricker
##' @docType data
##' @keywords datasets models
##' @family pomp examples
##' @examples
##'
##' pompExample(ricker)
##' plot(ricker)
##' coef(ricker)
##'
NULL
