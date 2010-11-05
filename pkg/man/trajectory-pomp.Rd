\name{trajectory}
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
\S4method{trajectory}{pomp}(object, params, times, t0, \dots)
}
\arguments{
  \item{object}{an object of class \code{pomp}.}
  \item{params}{
    a rank-2 array of parameters.
    Each column of \code{params} is a distinct parameter vector.
  }
  \item{times, t0}{
    \code{times} is a numeric vector specifying the times at which a trajectory is desired.
    \code{t0} specifies the start time (the time at which the initial conditions hold).
    The default for \code{times} is \code{times=time(object,t0=FALSE)} and \code{t0=timezero(object)}, respectively.
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

  When the skeleton is a vectorfield, \code{trajectory} integrates it using \code{\link[deSolve]{lsoda}}.

  Unresolved issue: What is the behavior if \code{type="map"} and \code{times} are non-integral?
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
data(euler.sir)
x <- trajectory(euler.sir)
plot(time(euler.sir),x["I",1,],type='l',xlab='time',ylab='I')
lines(time(euler.sir)[-1],diff(x["cases",1,]),col='red')

coef(euler.sir,c("gamma")) <- log(12)
x <- trajectory(euler.sir)
plot(time(euler.sir),x["I",1,],type='l',xlab='time',ylab='I')
lines(time(euler.sir)[-1],diff(x["cases",1,]),col='red')
}
\seealso{\code{\link{pomp-class}}, \code{\link{pomp}}, \code{\link[deSolve]{lsoda}}}
\keyword{models}
\keyword{ts}
