\name{skeleton-pomp}
\docType{methods}
\alias{skeleton}
\alias{skeleton,pomp-method}
\alias{skeleton-pomp}

\title{Evaluate the deterministic skeleton at the given points in state space.}
\description{
  The method \code{skeleton} computes the deterministic skeleton.
  In the case of a discrete-time system, this is the one-step map.
  In the case of a continuous-time system, this is the vector-field.
  NB: \code{skeleton} just evaluates the deterministic skeleton;
  it does not iterate or integrate.
}
\usage{
skeleton(object, x, t, params, \dots)
\S4method{skeleton}{pomp}(object, x, t, params, \dots)
}
\arguments{
  \item{object}{an object of class \code{pomp}.}
  \item{x}{
    a rank-3 array containing the states of the unobserved process at which the deterministic skeleton is to be evaluated.
    The dimensions of \code{x} are \code{nvars} x \code{nreps} x \code{ntimes}, where \code{nvars} is the number of state variables, \code{nreps} is the number of replicates, and \code{ntimes} is the length of \code{times}.
    Note that if \code{ntimes != length(times)}, an error is generated.
  }
  \item{t}{
    a numeric vector containing the times at which the deterministic skeleton is to be evaluated.
  }
  \item{params}{
    a rank-2 array of parameters with columns corresponding to the columns of \code{x}.
    Note that the \code{x} and \code{params} must agree in the number of their columns.
  }
  \item{\dots}{at present, these are ignored.}
}
\value{
  Returns an array of dimensions \code{nvar} x \code{nreps} x \code{ntimes}.
  If \code{f} is the returned matrix, \code{f[i,j,k]} is the i-th component of the deterministic skeleton at time \code{times[k]} given the state \code{x[,j,k]} and parameters \code{params[,j]}.
}
\details{
  This function makes repeated calls to the user-supplied \code{skeleton} of the \code{pomp} object.
  For specifications on supplying these, see \code{\link{pomp}}.
}
\references{}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{pomp-class}}, \code{\link{pomp}}}
\keyword{models}
\keyword{ts}
