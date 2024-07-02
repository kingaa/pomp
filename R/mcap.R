##' Monte Carlo adjusted profile
##'
##' Given a collection of points maximizing the likelihood over a range
##' of fixed values of a focal parameter, this function constructs
##' a profile likelihood confidence interval accommodating both
##' Monte Carlo error in the profile and statistical uncertainty present
##' in the likelihood function.
##'
##' @param logLik numeric; a vector of profile log likelihood evaluations.
##' @param parameter numeric; the corresponding values of the focal parameter.
##' @param level numeric; the confidence level required.
##' @param span numeric; the \code{\link[stats]{loess}} smoothing parameter.
##' @param Ngrid integer; the number of points to evaluate the smoothed profile.
##'
##' @references
##'
##' \Ionides2017
##'
##' @return
##' \code{mcap} returns a list including the \code{\link[stats]{loess}}-smoothed
##' profile, a quadratic approximation, and the constructed confidence interval.
##'
##' @author Edward L. Ionides
##' @concept profile likelihood
##' @importFrom stats loess predict lm vcov qchisq
##' @include package.R design.R
##'
##' @rdname mcap
##' @export
mcap <- function (
  logLik, parameter,
  level = 0.95,
  span = 0.75, Ngrid = 1000
)
{
  smooth_fit <- loess(logLik ~ parameter,span=span)
  parameter_grid <- seq(min(parameter),max(parameter),length.out=Ngrid)
  smoothed_logLik <- predict(smooth_fit,newdata=parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_logLik)]
  dist <- abs(parameter-smooth_arg_max)
  included <- dist < sort(dist)[trunc(span*length(dist))]
  maxdist <- max(dist[included])
  weights <- numeric(length(parameter))
  weights[included] <- (1-(dist[included]/maxdist)^3)^3
  quadratic_fit <- lm(
    logLik ~ a + b,
    weights=weights,
    data = data.frame(logLik=logLik,b=parameter,a=-parameter^2)
  )
  b <- unname(coef(quadratic_fit)["b"] )
  a <- unname(coef(quadratic_fit)["a"] )
  m <- vcov(quadratic_fit)
  var_b <- m["b","b"]
  var_a <- m["a","a"]
  cov_ab <- m["a","b"]
  se_mc_squared <- (1/(4*a*a))*(var_b-(2*b/a)*cov_ab+(b*b/a/a)*var_a)
  se_stat_squared <- 1/2/a
  se_total_squared <- se_mc_squared + se_stat_squared
  delta <- qchisq(level,df=1) * (a*se_mc_squared+0.5)
  logLik_diff <- max(smoothed_logLik)-smoothed_logLik
  ci <- range(parameter_grid[logLik_diff < delta])
  list(
    logLik=logLik,
    parameter=parameter,
    level=level,
    span=span,
    quadratic_fit=quadratic_fit,
    quadratic_max=b/(2*a),
    smooth_fit=smooth_fit,
    fit=data.frame(
      parameter=parameter_grid,
      smoothed=smoothed_logLik,
      quadratic=predict(
        quadratic_fit,
        newdata=list(
          b=parameter_grid,
          a=-parameter_grid^2
        )
      )
    ),
    mle=smooth_arg_max,
    ci=ci,
    delta=delta,
    se_stat=sqrt(se_stat_squared),
    se_mc=sqrt(se_mc_squared),
    se=sqrt(se_total_squared)
  )
}
