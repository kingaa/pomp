\name{ricker}
\alias{ricker}
\docType{data}
\title{Ricker model with Poisson observations.}
\description{
  \code{ricker} is a \code{pomp} object encoding a stochastic Ricker model with Poisson measurement error.
}
\details{
  The state process is \eqn{N_{t+1} = r N_{t} \exp(-N_{t}+e_{t})}{N[t+1] = r N[t] exp(-N[t]+e[t])}, where the \eqn{e_t}{e[t]} are i.i.d. normal random deviates with zero mean and variance \eqn{\sigma^2}{sigma^2}.
  The observed variables \eqn{y_t}{y[t]} are distributed as \eqn{\mathrm{Poisson}(\phi N_t)}{Poisson(phi N[t])}.
}
\examples{
pompExample(ricker)
plot(ricker)
coef(ricker)
}
\seealso{\code{\link{pomp}}, \code{\link{gompertz}}, and the tutorials on the \href{https://kingaa.github.io/pomp}{package website}.}
\keyword{datasets}
\keyword{models}
