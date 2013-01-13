\name{budmoth}
\alias{budmoth.sim}
\title{Larch budmoth model POMPs with real and simulated data.}
\description{
  \code{budmoth.sim} constructs a \code{pomp} object containing the larch budmoth model and simulated budmoth density, parasitism rate, and food quality (needle-length) data.
  Four datasets, representing four distinct parameter regimes, are avaiable.
}
\usage{
budmoth.sim(which)
}
\arguments{
  \item{which}{
    dataset to load given as a name or literal character string.
    Evoked without an argument, \code{budmoth.sim} lists all available datasets.
  }
}
\examples{
budmoth.sim()			## print a list of all available datasets
## three regimes, high and low noise regimes for parasitism and tritrophic
plot(budmoth.sim(food))
plot(budmoth.sim(para1))
plot(budmoth.sim(para2))
plot(budmoth.sim(tri))
}
\seealso{the dQuote{budmoth-model} vignette}
\keyword{datasets}
