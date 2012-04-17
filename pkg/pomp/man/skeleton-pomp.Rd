\name{skeleton-pomp}
\docType{methods}
\alias{skeleton}
\alias{skeleton,pomp-method}
\alias{skeleton-pomp}
\keyword{internal}
\title{Evaluate the deterministic skeleton at the given points in state space.}
\description{
  The method \code{skeleton} computes the deterministic skeleton.
  In the case of a discrete-time system, this is the one-step map.
  In the case of a continuous-time system, this is the vector-field.
  NB: \code{skeleton} just evaluates the deterministic skeleton;
  it does not iterate or integrate.
  This function is part of the low-level interface to \code{pomp} objects.
  This help page does not give instructions on the implementation of models: see \code{\link{pomp}} for instructions.
}
\usage{
\S4method{skeleton}{pomp}(object, x, t, params, \dots)
}
\arguments{
  \item{object}{an object of class \code{pomp}.}
  \item{x}{
    a rank-3 array containing the states of the unobserved process at which the deterministic skeleton is to be evaluated.
    The dimensions of \code{x} are \code{nvars} x \code{nreps} x \code{ntimes}, where \code{nvars} is the number of state variables, \code{nreps} is the number of replicates, and \code{ntimes} is the length of \code{times}.
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
  For specifications on supplying this, see \code{\link{pomp}}.

  The function \code{\link{trajectory}} iterates the skeleton (in case it is a map) or integrates it using an ODE solver (in case it is a vectorfield).
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{pomp-class}}, \code{\link{pomp}}, \code{\link{trajectory}}}
\keyword{models}
\keyword{ts}
