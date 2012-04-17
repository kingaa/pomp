\name{budmoth}
\alias{budmoth.eng}
\alias{budmoth.sim}
\docType{data}
\title{Larch budmoth model POMPs with real and simulated data.}
\description{
  \code{budmoth.eng} is a \code{pomp} object containing the larch budmoth model, budmoth density, parasitism rate, and food quality (needle-length) data from the Engadine valley of Switzerland.
  \code{budmoth.sim} is a list of \code{pomp} objects:
  each encodes the larch budmoth model and contains simulated data.
  Three parameter regimes are represented, with a number of replicates for each.
}
\usage{
data(budmoth.eng)
data(budmoth.sim)
}
\examples{
data(budmoth.eng)
plot(budmoth.eng)

data(budmoth.sim)
names(budmoth.sim)
## three regimes, high and low noise regimes for parasitism and tritrophic
plot(budmoth.sim[["food.1"]])
plot(budmoth.sim[["para.1"]])
plot(budmoth.sim[["para.2"]])
plot(budmoth.sim[["tri.1"]])
}
\seealso{the vignettes}
\keyword{datasets}
