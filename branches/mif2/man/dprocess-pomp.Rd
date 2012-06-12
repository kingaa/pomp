\name{dprocess-pomp}
\docType{methods}
\alias{dprocess}
\alias{dprocess,pomp-method}
\alias{dprocess-pomp}
\keyword{internal}
\title{
  Evaluate the probability density of state transitions in a Markov process
}
\description{
  The method \code{dprocess} evaluates the probability density of a set of consecutive state transitions.
  This function is part of the low-level interface to \code{pomp} objects.
  This help page does not give instructions on the implementation of models: see \code{\link{pomp}} for instructions.
}
\usage{
\S4method{dprocess}{pomp}(object, x, times, params, log = FALSE, \dots)
}
\arguments{
  \item{object}{an object of class \code{pomp}.}
  \item{x}{
    a rank-3 array containing the states of the unobserved process.
    The dimensions of \code{x} are \code{nvars} x \code{nreps} x \code{ntimes}, where \code{nvars} is the number of state variables, \code{nreps} is the number of replicates, and \code{ntimes} is the length of \code{times}.
  }
  \item{times}{
    a numeric vector containing the times corresponding to the given states.
  }
  \item{params}{
    a rank-2 array of parameters with columns corresponding to the columns of \code{x}.
    Note that the \code{x} and \code{params} must agree in the number of their columns.
  }
  \item{log}{if TRUE, log probabilities are returned.}
  \item{\dots}{at present, these are ignored.}
}
\value{
  Returns a matrix of dimensions \code{nreps} x \code{ntimes-1}.
  If \code{d} is the returned matrix, \code{d[j,k]} is the likelihood of the transition from state \code{x[,j,k-1]} at time \code{times[k-1]} to state \code{x[,j,k]} at time \code{times[k]}.
}
\details{
  This function is essentially a wrapper around the user-supplied \code{dprocess} slot of the \code{pomp} object.
  For specifications on writing such a function, see \code{\link{pomp}}.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{pomp-class}}, \code{\link{pomp}}}
\keyword{models}
\keyword{ts}
