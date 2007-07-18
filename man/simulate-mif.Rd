\name{simulate-mif}
\docType{methods}
\alias{simulate,mif-method}
\alias{simulate-mif}
\title{Simulating a fitted "mif" object}
\description{\code{simulate} can be used to generate simulated
  data sets and/or to simulate the state process.}
\usage{
\S4method{simulate}{mif}(object, nsim = 1, seed = NULL, \dots)
}
\arguments{
  \item{object}{An object of class \code{pomp}.}
  \item{nsim}{The number of simulations to perform.  Note that the
    number of replicates will be \code{nsim} times
    \code{ncol(xstart)}.}
  \item{seed}{The random seed to use.}
  \item{\dots}{further arguments are passed to the \code{simulate}
    method for the \code{pomp} object.  See \link{simulate-pomp}.}
}
\value{
  See \link{simulate-pomp}.
}
\details{
  When applied to a \code{mif} object, the underlying \code{pomp} object
  is simulated using the estimated parameters.  See \link{simulate-pomp}
  for details.
}
\references{}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{simulate}}, \link{simulate-pomp}, \link{mif-class}}
\keyword{models}
\keyword{ts}
