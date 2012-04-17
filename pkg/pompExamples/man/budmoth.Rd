\name{budmoth}
\alias{budmoth.sim}
\docType{data}
\title{Larch budmoth model POMPs with real and simulated data.}
\description{
  \code{budmoth.sim} is a list of \code{pomp} objects containing the larch budmoth model, budmoth density, parasitism rate, and food quality (needle-length) simulated data.
  Four parameter regimes are represented.
}
\usage{
data(budmoth.sim)
}
\examples{
data(budmoth.sim)
names(budmoth.sim)
## three regimes, high and low noise regimes for parasitism and tritrophic
plot(budmoth.sim[["food"]])
plot(budmoth.sim[["para1"]])
plot(budmoth.sim[["para2"]])
plot(budmoth.sim[["tri"]])
}
\seealso{the vignettes}
\keyword{datasets}
