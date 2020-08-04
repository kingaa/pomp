##' Gompertz model with log-normal observations.
##'
##' \code{gompertz()} constructs a \sQuote{pomp} object encoding a stochastic Gompertz population model with log-normal measurement error.
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
##' @keywords models pomp_datasets
##' @include simulate.R
##' @family pomp_examples
##'
##' @return
##' A \sQuote{pomp} object with simulated data.
##'
##' @examples
##'
##' plot(gompertz())
##' plot(gompertz(K=2,r=0.01))
##'
NULL

##' @rdname gompertz
##'
##' @param r growth rate
##' @param K carrying capacity
##' @param sigma process noise intensity
##' @param tau measurement error s.d.
##' @param X_0 value of the latent state variable \code{X} at the zero time
##' @param t0 zero time
##' @param times observation times
##'
##' @export
gompertz <- function (K = 1, r = 0.1, sigma = 0.1, tau = 0.1, X_0 = 1,
  times = 1:100, t0 = 0)
{

  simulate(
    times=times, t0=t0,
    params=c(K=K,r=r,sigma=sigma,tau=tau,X_0=X_0),
    partrans=parameter_trans(
      toEst="_gompertz_to_trans",
      fromEst="_gompertz_from_trans"
    ),
    rprocess=discrete_time(
      step.fun="_gompertz_simulator"
    ),
    rmeasure="_gompertz_normal_rmeasure",
    dmeasure="_gompertz_normal_dmeasure",
    skeleton=map("_gompertz_skeleton",delta.t=1),
    PACKAGE="pomp",
    paramnames=c("r","K","sigma","tau","X_0"),
    obsnames=c("Y"),
    statenames=c("X"),
    seed=299438676L
  )

}
