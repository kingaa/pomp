\name{trajectory-pomp}
\docType{methods}
\alias{trajectory}
\alias{trajectory,pomp-method}
\alias{trajectory-pomp}

\title{Compute trajectories of the determinstic skeleton.}
\description{
  The method \code{trajectory} computes a trajectory of the deterministic skeleton of a Markov process.
  In the case of a discrete-time system, the deterministic skeleton is a map and a trajectory is obtained by iterating the map.
  In the case of a continuous-time system, the deterministic skeleton is a vector-field; \code{trajectory} integrates the vectorfield to obtain a trajectory.
}
\usage{
trajectory(object, params, times, \dots)
\S4method{trajectory}{pomp}(object, params, times, \dots)
}
\arguments{
  \item{object}{an object of class \code{pomp}.}
  \item{times}{
    a numeric vector containing the times at which a trajectory is desired.
    The first of these will be the initial time.
    By default, \code{times=time(object,t0=TRUE)}.
  }
  \item{params}{
    a rank-2 array of parameters.
    Each column of \code{params} is a distinct parameter vector.
  }
  \item{\dots}{at present, these are ignored.}
}
\value{
  Returns an array of dimensions \code{nvar} x \code{nreps} x \code{ntimes}.
  If \code{x} is the returned matrix, \code{x[i,j,k]} is the i-th component of the state vector at time \code{times[k]} given parameters \code{params[,j]}.
}
\details{
  This function makes repeated calls to the user-supplied \code{skeleton} of the \code{pomp} object.
  For specifications on supplying this, see \code{\link{pomp}}.

  What is the behavior if \code{type="map"} and \code{times} are non-integral?
}
\references{}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{pomp-class}}, \code{\link{pomp}}}
\keyword{models}
\keyword{ts}
