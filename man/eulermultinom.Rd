\name{Euler-multinomial models}
\alias{eulermultinom}
\alias{reulermultinom}
\alias{deulermultinom}
\title{Euler-multinomial models}
\description{
  Density and random generation for the Euler-multinomial death process with parameters \code{size}, \code{rate}, and \code{dt}.
}
\usage{
reulermultinom(n = 1, size, rate, dt)
deulermultinom(x, size, rate, dt, log = FALSE)
}
\arguments{
  \item{n}{Number of random variates to generate.}
  \item{size}{Number of individuals at risk.}
  \item{rate}{Hazard rates.}
  \item{dt}{Duration of Euler step.}
  \item{x}{Matrix or vector containing number of individuals that have succumbed to each death process.}
  \item{log}{logical; if TRUE, return logarithm of probabilities.}
}
\value{
  \item{reulermultinom}{
    Returns a \code{length(rate)} by \code{n} matrix.
    Each column is a different random draw.
    Each row contains the numbers of individuals succumbed to the corresponding process.
  }
  \item{deulermultinom}{
    Returns a vector (of length equal to the number of columns of \code{x}) containing the probabilities of observing \code{x} given the specified parameters (\code{size}, \code{rate}, \code{dt}).
  }
}
\details{
  Direct access to the underlying C routines is available: see the header file "euler.h", included with the package.
}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{euler}}, \code{\link{pomp}}}
\examples{
x <- reulermultinom(5,size=100,rate=c(a=1,b=2,c=3),dt=0.1)
deulermultinom(x,size=100,rate=c(1,2,3),dt=0.1)
}
\keyword{distribution}
