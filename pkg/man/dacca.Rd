\name{dacca}
\alias{dacca}
\alias{dacca.transform}
\docType{data}
\title{Model of cholera transmission for historic Bengal.}
\description{
  \code{dacca} is a \code{pomp} object containing census and cholera mortality data from the Dacca district of the former British province of Bengal over the years 1891 to 1940 together with a stochastic differential equation transmission model.
  The model is that of King et al. (2008).
  It also has the MLE for the SIRS model with seasonal reservoir.

  Data are provided courtesy of Dr. Menno J. Bouma, London School of Tropical Medicine and Hygiene.
}
\usage{
data(dacca)
}
\details{
  \code{dacca} is a \code{pomp} object containing the model, data, and MLE parameters.
  \code{dacca.transform} is a function that transforms the parameters.
}
\examples{
data(dacca)
plot(dacca)
theta <- dacca.transform(coef(dacca),dir="inverse") #MLEs
plot(simulate(dacca))
theta["eps"] <- 1
plot(simulate(dacca,params=dacca.transform(theta,dir="forward")))
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
