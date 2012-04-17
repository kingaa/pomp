\name{pertussis}
\alias{pertussis.sim}
\docType{data}
\title{Pertussis models with simulated data.}
\description{
  \code{pertussis.sim} is a list of \code{pomp} objects:
  each encodes one of the pertussis models and contains simulated data.
}
\usage{
data(pertussis.sim)
}
\examples{
data(pertussis.sim)
names(pertussis.sim)

plot(pertussis.sim[["SEIR.small"]])
plot(pertussis.sim[["SEIRS.big"]])
plot(pertussis.sim[["SEIRR.small"]])
}
\seealso{the vignettes}
\keyword{datasets}
