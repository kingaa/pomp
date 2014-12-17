\name{budmoth}
\docType{data}
\alias{budmoth.sim}
\title{Larch budmoth model POMPs with real and simulated data.}
\description{
  \code{pompExample(budmoth.sim)} constructs a \code{pomp} object containing the larch budmoth model and simulated budmoth density, parasitism rate, and food quality (needle-length) data.
  Four datasets, representing four distinct parameter regimes, are avaiable.
}
\examples{
## three regimes, high and low noise regimes for parasitism and tritrophic
bm <- pompExample(budmoth,envir=NULL)
plot(bm$food)
plot(bm$para1)
plot(bm$para2)
plot(bm$tri)
}
\keyword{datasets}
