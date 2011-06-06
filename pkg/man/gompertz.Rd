\name{gompertz}
\alias{gompertz}
\docType{data}
\title{Gompertz model with log-normal observations.}
\description{
  \code{gompertz} is a \code{pomp} object encoding a stochastic Gompertz population model with log-normal measurement error.
}
\usage{data(gompertz)}
\details{
  The state process is \eqn{X_{t+1} = K^{(1-S)} X_{t}^S \varepsilon_{t}}, where \eqn{S=e^{-r}} and the \eqn{\varepsilon_t} are i.i.d. lognormal random deviates with variance \eqn{\sigma^2}.
  The observed variables \eqn{Y_t} are distributed as \eqn{\mathrm{lognormal}(\log{X_t},\tau)}.
  The model is parameterized by the logarithms of \eqn{r}, \eqn{K}, \eqn{\sigma}, and \eqn{\tau};
  the initial condition is parameterized directly.
}
\examples{
data(gompertz)
plot(gompertz)
coef(gompertz)
}
\seealso{
  \code{\link{pomp-class}} and the introductory vignette \code{vignette("intro_to_pomp")}.
}
\keyword{datasets}
