\name{sir}
\title{Compartmental epidemiological models}
\alias{euler.sir}
\alias{gillespie.sir}
\alias{bbs}
\docType{data}
\description{
  \code{euler.sir} is a \code{pomp} object encoding a simple seasonal SIR model.
  Simulation is performed using an Euler multinomial approximation.
  \code{gillespie.sir} has the same model implemented using Gillespie's algorithm.
  \code{bbs} is a nonseasonal SIR model together with data from a 1978 outbreak of influenza in a British boarding school.
}
\details{
  This and similar examples are discussed and constructed in tutorials available on the \href{https://kingaa.github.io/pomp}{package website}.

  The boarding school influenza outbreak is described in Anonymous (1978).
}
\examples{
pompExample(euler.sir)
plot(euler.sir)
plot(simulate(euler.sir,seed=20348585))
coef(euler.sir)

pompExample(gillespie.sir)
plot(gillespie.sir)
plot(simulate(window(gillespie.sir,end=3),seed=20348585))
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
\seealso{\code{\link{pomp}} and the tutorials on the \href{https://kingaa.github.io/pomp}{package website}.}
\keyword{datasets}
\keyword{models}
