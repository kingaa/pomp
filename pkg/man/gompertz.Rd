\name{gompertz}
\alias{gompertz}
\docType{data}
\title{Gompertz model with normal observations.}
\description{
  \code{gompertz} is a \code{pomp} object encoding a stochastic Gompertz population model with normal measurement error.
}
\usage{data(gompertz)}
\details{
  The state process is \eqn{X_{t+1} = K^(1-R) X_{t}^R eps_{t}}, where the \eqn{eps_t} are i.i.d. lognormal random deviates with variance \eqn{sigma^2}.
  The observed variables \eqn{Y_t} are distributed as \eqn{normal(X_t,tau)}.
}
\examples{
data(gompertz)
plot(gompertz)
coef(gompertz)
}
\seealso{\code{\link{pomp-class}} and the vignettes}
\keyword{datasets}
