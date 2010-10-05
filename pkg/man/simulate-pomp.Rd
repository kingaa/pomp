\name{simulate-pomp}
\docType{methods}
\alias{simulate,pomp-method}
\alias{simulate-pomp}
\title{Running simulations of a partially-observed Markov process}
\description{
  \code{simulate} can be used to generate simulated data sets and/or to simulate the state process.
}
\usage{
\S4method{simulate}{pomp}(object, nsim = 1, seed = NULL, params,
         states = FALSE, obs = FALSE, times, t0, \dots)
}
\arguments{
  \item{object}{An object of class \code{pomp}.}
  \item{nsim}{
    The number of simulations to perform.
    Note that the number of replicates will be \code{nsim} times \code{ncol(xstart)}.
  }
  \item{seed}{
    optional; if set, the pseudorandom number generator (RNG) will be initialized with \code{seed}.
    the random seed to use.
    The RNG will be restored to its original state afterward.
  }
  \item{params}{
    The parameters to use in simulating the model.
    If \code{params} is not given, then the contents of the \code{params} slot of \code{object} will be used, if they exist.
  }
  \item{states}{Do we want the state trajectories?}
  \item{obs}{Do we want data-frames of the simulated observations?}
  \item{times, t0}{
    \code{times} specifies the times at which simulated observations will be made.
    \code{t0} specifies the start time (the time at which the initial conditions hold).
    The default for \code{times} is is \code{times=time(object,t0=FALSE)} and \code{t0=timezero(object)}, respectively.
  }
  \item{\dots}{further arguments that are currently ignored.}
}
\value{
  If \code{states=FALSE} and \code{obs=FALSE} (the default), a list of \code{nsim} \code{pomp} objects is returned.
  Each has a simulated data set, together with the parameters used (in slot \code{params}) and the state trajectories also (in slot \code{states}).
  If \code{times} is specified, then the simulated observations will be at times \code{times}.
  
  If \code{nsim=1}, then a single \code{pomp} object is returned (and not a singleton list).

  If \code{states=TRUE} and \code{obs=FALSE}, simulated state trajectories are returned as a rank-3 array with dimensions \code{nvar} x \code{(ncol(params)*nsim)} x \code{ntimes}.
  Here, \code{nvar} is the number of state variables and \code{ntimes} the length of the argument \code{times}.
  The measurement process is not simulated in this case.   

  If \code{states=FALSE} and \code{obs=TRUE}, simulated observations are returned as a rank-3 array with dimensions \code{nobs} x \code{(ncol(params)*nsim)} x \code{ntimes}.
  Here, \code{nobs} is the number of observables.

  If both \code{states=TRUE} and \code{obs=TRUE}, then a named list is returned.
  It contains the state trajectories and simulated observations as above.
}
\details{
  Simulation of the state process and of the measurement process are each accomplished by a single call to the user-supplied \code{rprocess} and \code{rmeasure} functions, respectively.
  This makes it possible for the user to write highly optimized code for these potentially expensive computations.
}
\examples{
data(ou2)
x <- simulate(ou2,seed=3495485,nsim=10)
x <- simulate(ou2,seed=3495485,nsim=10,states=TRUE,obs=TRUE)
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\link{pomp-class}}
\keyword{models}
\keyword{ts}
