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
\details{
  This example is discussed tutorials available on the \href{http://pomp.r-forge.r-project.org}{package website}.

  The codes that construct these \code{pomp} objects can be found in the \dQuote{examples} directory in the installed package.
  Do \code{system.file("examples",package="pomp"))} to find this directory.
  For the basic \code{rprocess}, \code{dmeasure}, \code{rmeasure}, and \code{skeleton} functions, these codes use compiled native routines built into the package's library.
  View \dQuote{src/sir.c} in the package source or \code{file.show("examples/sir.c")} from an \R session to view these codes.

  The boarding school influenza outbreak is described in Anonymous (1978).
}
\examples{
pompExample(euler.sir)
plot(euler.sir)
plot(simulate(euler.sir,seed=20348585))
coef(euler.sir)

pompExample(gillespie.sir)
plot(gillespie.sir)
plot(simulate(gillespie.sir,seed=20348585))
coef(gillespie.sir)

pompExample(bbs)
plot(bbs)
coef(bbs)
}
\references{
  Anonymous (1978).
  Influenza in a boarding school.
  British Medical Journal 1:587
}
\seealso{\code{\link{pomp}} and the tutorials on the \href{http://pomp.r-forge.r-project.org}{package website}.}
\keyword{datasets}
\keyword{models}
