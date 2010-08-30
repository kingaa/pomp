\name{dacca}
\alias{dacca}
\docType{data}
\title{Model of cholera transmission for historic Bengal.}
\description{
  \code{dacca} contains census and cholera mortality data from the Dacca district of the former British province of Bengal over the years 1891 to 1940 together with a stochastic differential equation transmission model.
  The model is that of King et al. (2008).
  It also has the MLE for the SIRS model with seasonal reservoir.

  Data are provided courtesy of Dr. Menno J. Bouma, London School of Tropical Medicine and Hygiene.
}
\usage{
data(dacca)
}
\details{
  \code{dacca} is a list consisting of
  \describe{
    \item{po}{
      A \code{pomp} object containing the model, data, and MLE parameters.
    }
    \item{transform}{
      A function for transforming the parameters.
    }
  }
}
\examples{
data(dacca)
plot(dacca$po)
theta <- dacca$transform(coef(dacca$po),dir="inverse") #MLEs
plot(simulate(dacca$po))
theta["eps"] <- 1
plot(simulate(dacca$po,params=dacca$transform(theta,dir="forward")))
}
\references{
  King, A. A., Ionides, E. L., Pascual, M., and Bouma, M. J.
  Inapparent infections and cholera dynamics.
  Nature 454:877-880 (2008)
}
\seealso{
  \code{\link{euler.sir}},
  \code{\link{pomp}}
}
\keyword{datasets}
