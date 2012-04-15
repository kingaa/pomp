\name{sir}
\alias{euler.sir}
\alias{gillespie.sir}
\alias{bbs}
\docType{data}
\title{SIR models.}
\description{
  \code{euler.sir} is a \code{pomp} object encoding a simple seasonal SIR model.
  Simulation is performed using an Euler multinomial approximation.
  \code{gillespie.sir} has the same model implemented using Gillespie's algorithm.
  \code{bbs} is a nonseasonal SIR model together with data from a 1978 outbreak of influenza in a British boarding school.
}
\usage{
data(euler.sir)
data(gillespie.sir)
data(bbs)
}
\details{
  This example is discussed extensively in the \dQuote{Introduction to \pkg{pomp}} and \dQuote{Advanced topics in \pkg{pomp}} vignettes.

  The codes that construct these \code{pomp} objects can be found in the \dQuote{data-R} directory in the installed package.
  Do \code{file.show(system.file("data-R/sir.R",package="pomp"))} to view these codes.
  For the basic \code{rprocess}, \code{dmeasure}, \code{rmeasure}, and \code{skeleton} functions, these codes use compiled native routines built into the package's library.
  View \dQuote{src/sir.c} in the package source or \code{file.show("examples/sir.c")} from an \R session to view these codes.

  The boarding school influenza outbreak is described in Anonymous (1978).
}
\examples{
data(euler.sir)
plot(euler.sir)
plot(simulate(euler.sir,seed=20348585))
coef(euler.sir)

data(gillespie.sir)
plot(gillespie.sir)
plot(simulate(gillespie.sir,seed=20348585))
coef(gillespie.sir)

data(bbs)
plot(bbs)
coef(bbs)
}
\references{
  Anonymous (1978).
  Influenza in a boarding school.
  British Medical Journal 1:587
}
\seealso{\code{\link{pomp-class}} and the vignettes}
\keyword{datasets}
