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
\S4method{trajectory}{pomp}(object, params, times, t0, as.data.frame = FALSE, \dots)
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
  \item{as.data.frame}{
    logical; if \code{TRUE}, return the result as a data-frame.
  }
  \item{\dots}{
    additional arguments are passed to the ODE integrator if the skeleton is a vectorfield and ignored if it is a map.
    See \code{\link[deSolve]{ode}} for a description of the additional arguments accepted.
  }
}
\value{
  Returns an array of dimensions \code{nvar} x \code{nreps} x \code{ntimes}.
  If \code{x} is the returned matrix, \code{x[i,j,k]} is the i-th component of the state vector at time \code{times[k]} given parameters \code{params[,j]}.
}
\details{
  This function makes repeated calls to the user-supplied \code{skeleton} of the \code{pomp} object.
  For specifications on supplying this, see \code{\link{pomp}}.

  When the skeleton is a vectorfield, \code{trajectory} integrates it using \code{\link[deSolve]{ode}}.

  When the skeleton is a map, \code{trajectory} iterates it.
  By default, time is advanced 1 unit per iteration.
  The user can change this behavior by specifying the desired timestep using the argument \code{skelmap.delta.t} in the construction of the \code{pomp} object.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
data(euler.sir)
x <- trajectory(euler.sir)
plot(time(euler.sir),x["I",1,],type='l',xlab='time',ylab='I')
lines(time(euler.sir),x["cases",1,],col='red')

coef(euler.sir,c("gamma")) <- log(12)
x <- trajectory(euler.sir)
plot(time(euler.sir),x["I",1,],type='l',xlab='time',ylab='I')
lines(time(euler.sir),x["cases",1,],col='red')

x <- trajectory(euler.sir,as.data.frame=TRUE)
}
\seealso{\code{\link{pomp}}, \code{\link{traj.match}}, \code{\link[deSolve]{ode}}}
\keyword{models}
\keyword{ts}
