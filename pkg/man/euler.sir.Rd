\name{euler.sir}
\alias{euler.sir}
\docType{data}
\title{Seasonal SIR model implemented as an Euler-multinomial model}
\description{
  \code{euler.sir} is a \code{pomp} object encoding a simple seasonal SIR model.
}
\usage{data(euler.sir)}
\examples{
data(euler.sir)
plot(euler.sir)
x <- simulate(euler.sir,nsim=10,seed=20348585)
plot(x[[1]])
}
\seealso{\code{\link{pomp-class}} and the vignettes}
\keyword{datasets}
