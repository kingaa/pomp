\name{pertussis}
\alias{pertussis.sim}
\title{Pertussis models with simulated data.}
\description{
  \code{pertussis.sim} constructs a \code{pomp} object containing one of the pertussis models and simulated incidence data.
}
\usage{
pertussis.sim(which)
}
\arguments{
  \item{which}{
    dataset to load given as a name or literal character string.
    Evoked without an argument, \code{pertussis.sim} lists all available datasets.
  }
}
\examples{
pertussis.sim()
plot(pertussis.sim(SEIR.small))
plot(pertussis.sim(SEIRS.big))
plot(pertussis.sim(SEIRR.small))
plot(pertussis.sim("full.big"))
}
\seealso{the \dQuote{pertussis-model} vignette}
\keyword{models}
