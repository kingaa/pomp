\name{rmeasure-pomp}
\docType{methods}
\alias{rmeasure}
\alias{rmeasure,pomp-method}
\alias{rmeasure-pomp}
\title{Simulate the measurement model of a partially-observed Markov process}
\description{
  The method \code{rmeasure} draws from the distribution of measurements given the state of the system.
  This function is part of the low-level interface to \code{pomp} objects.
  This help page does not give instructions on the implementation of models: see \code{\link{pomp}} for instructions.
}
\usage{
rmeasure(object, x, times, params, \dots)
\S4method{rmeasure}{pomp}(object, x, times, params, \dots)
}
\arguments{
  \item{object}{an object of class \code{pomp}.}
  \item{x}{
    a rank-3 array containing the states of the unobserved process.
    The dimensions of \code{x} are \code{nvars} x \code{nreps} x \code{ntimes}, where \code{nvars} is the number of state variables, \code{nreps} is the number of simulations, and \code{ntimes} is the number of distinct times.
  }
  \item{times}{
    a numerical vector containing the times at which the measurements are to be made.
  }
  \item{params}{
    a rank-2 array (dimensions \code{npars} x \code{nreps}) of parameters with the parameters corresponding to the columns of \code{x}.
  }
  \item{\dots}{
    at present, these are ignored.
  }
}
\value{
  Returns a rank-3 array of dimensions \code{nobs} x \code{nreps} x \code{ntimes}, where \code{nobs} is the number of observed variables.
}
\details{
  This function is essentially a wrapper around the user-supplied \code{rmeasure} slot of the \code{pomp} object.
  For specifications on writing such a function, see \code{\link{pomp}}.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{pomp-class}}, \code{\link{pomp}}}
\keyword{models}
\keyword{ts}
