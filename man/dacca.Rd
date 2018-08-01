\name{dacca}
\alias{dacca}
\docType{data}
\title{Model of cholera transmission for historic Bengal.}
\description{
  \code{dacca} is a \sQuote{pomp} object containing census and cholera mortality data from the Dacca district of the former British province of Bengal over the years 1891 to 1940 together with a stochastic differential equation transmission model.
  The model is that of King et al. (2008).
  The parameters are the MLE for the SIRS model with seasonal reservoir.

  Data are provided courtesy of Dr. Menno J. Bouma, London School of Tropical Medicine and Hygiene.
}
\details{
  \code{dacca} is a \sQuote{pomp} object containing the model, data, and MLE parameters.
  Parameters that naturally range over the positive reals are log-transformed;
  parameters that range over the unit interval are logit-transformed;
  parameters that are naturally unbounded or take integer values are not transformed.
}
\examples{
pompExample(dacca)
plot(dacca)
#MLEs on the natural scale
coef(dacca)
#MLEs on the transformed scale
coef(dacca,transform=TRUE)
plot(simulate(dacca))
}
\references{
  A. A. King, E. L. Ionides, M. Pascual, and M. J. Bouma,
  Inapparent infections and cholera dynamics,
  Nature, 454:877-880, 2008
}
\seealso{
  \code{\link{sir}},
  \code{\link{pomp}}
}
\keyword{models}
\keyword{datasets}
