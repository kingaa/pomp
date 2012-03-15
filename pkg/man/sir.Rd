\name{sir}
\alias{sir}
\alias{euler.sir}
\alias{gillespie.sir}
\docType{data}
\title{Seasonal SIR model implemented using two stochastic simulation algorithms.}
\description{
  \code{euler.sir} is a \code{pomp} object encoding a simple seasonal SIR model.
  Simulation is performed using an Euler multinomial approximation (AKA tau leap method).
  \code{gillespie.sir} has the same model implemented using Gillespie's algorithm.
}
\usage{
data(euler.sir)
data(gillespie.sir)
}
\examples{
data(euler.sir)
plot(euler.sir)
x <- simulate(euler.sir,nsim=10,seed=20348585)
plot(x[[1]])

data(gillespie.sir)
plot(gillespie.sir)
x <- simulate(gillespie.sir,nsim=1,seed=20348585)
plot(x)
}
\seealso{\code{\link{pomp-class}} and the vignettes}
\keyword{datasets}
