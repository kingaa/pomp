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
dacca.transform(params,dir=c("forward","inverse"))
}
\details{
  \code{dacca} is a \code{pomp} object containing the model, data, and MLE parameters.
  \code{dacca.transform} is a function that transforms the parameters.
  If \code{params} is a matrix or vector of parameters on the natural scale for this model, \code{dacca.transform(params,dir="forward")} is the corresponding matrix or vector of parameters on the scale that they are used internally.
  If \code{params} is a matrix or vector of parameters on the internal scale for this model, \code{dacca.transform(params,dir="inverse")} is the corresponding matrix or vector of parameters on the natural scale.
  Parameters that naturally range over the positive reals are log-transformed by \code{dacca.transform};
  parameters that range over the unit interval are logit-transformed;
  parameters that are naturally unbounded or take integer values are not transformed.
}
\examples{
data(dacca)
plot(dacca)
#MLEs on the natural scale
theta <- dacca.transform(coef(dacca),dir="inverse") 
plot(simulate(dacca))
# now change 'eps' and simulate again
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
