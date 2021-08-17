##' Two-dimensional discrete-time Ornstein-Uhlenbeck process
##'
##' \code{ou2()} constructs a \sQuote{pomp} object encoding a bivariate discrete-time Ornstein-Uhlenbeck process with noisy observations.
##'
##' If the state process is \eqn{X(t) = (x_{1}(t),x_{2}(t))}, then
##' \deqn{X(t+1) = \alpha X(t) + \sigma \epsilon(t),}
##' where \eqn{\alpha} and \eqn{\sigma} are 2x2 matrices,
##' \eqn{\sigma} is lower-triangular, and
##' \eqn{\epsilon(t)} is standard bivariate normal.
##' The observation process is \eqn{Y(t) = (y_1(t),y_2(t))}, where
##' \eqn{y_i(t) \sim \mathrm{normal}(x_i(t),\tau)}.
##'
##' @name ou2
##' @rdname ou2
##' @include simulate.R
##' @docType data
##' @keywords models pomp_datasets
##' @family pomp examples
##' @return A \sQuote{pomp} object with simulated data.
##' @example examples/ou2.R
##'
NULL

##' @rdname ou2
##' @param alpha_1,alpha_2,alpha_3,alpha_4 entries of the \eqn{alpha} matrix, in column-major order.
##' That is, \code{alpha_2} is in the lower-left position.
##' @param sigma_1,sigma_2,sigma_3 entries of the lower-triangular \eqn{sigma} matrix.
##' \code{sigma_2} is the entry in the lower-left position.
##' @param tau measurement error s.d.
##' @param x1_0,x2_0 latent variable values at time \code{t0}
##' @param t0 the zero time
##' @param times vector of observation times
##'
##' @export
ou2 <- function (
  alpha_1 = 0.8, alpha_2 = -0.5, alpha_3 = 0.3, alpha_4 = 0.9,
  sigma_1 = 3, sigma_2 = -0.5, sigma_3 = 2,
  tau = 1,
  x1_0 = -3, x2_0 = 4,
  times = 1:100, t0 = 0
)
{

  simulate(
    times=times, t0=t0,
    seed=787832394,
    rprocess=discrete_time("_ou2_step"),
    dprocess="_ou2_pdf",
    dmeasure = "_ou2_dmeasure",
    rmeasure = "_ou2_rmeasure",
    skeleton = map("_ou2_skel",delta.t=1),
    PACKAGE="pomp",
    paramnames = c(
      "alpha_1","alpha_2","alpha_3","alpha_4",
      "sigma_1","sigma_2","sigma_3",
      "tau"
    ),
    statenames = c("x1","x2"),
    obsnames = c("y1","y2"),
    params=c(
      alpha_1=alpha_1, alpha_2=alpha_2, alpha_3=alpha_3, alpha_4=alpha_4,
      sigma_1=sigma_1, sigma_2=sigma_2, sigma_3=sigma_3,
      tau=tau,
      x1_0=x1_0, x2_0=x2_0
    )
  )

}
