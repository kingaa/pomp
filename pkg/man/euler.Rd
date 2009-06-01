\name{euler}
\alias{euler}
\alias{euler.simulate}
\alias{onestep.simulate}
\alias{onestep.density}
\title{Plug-ins for dynamical models based on stochastic Euler algorithms}
\description{
  Facilities for implementing discrete-time Markov processes and continuous-time Markov processes using the Euler algorithm.
}
\usage{
euler.simulate(xstart, times, params, step.fun, delta.t, \dots,
               statenames = character(0), paramnames = character(0),
               covarnames = character(0), zeronames = character(0),
               tcovar, covar, PACKAGE)
onestep.simulate(xstart, times, params, step.fun, \dots,
                 statenames = character(0), paramnames = character(0),
                 covarnames = character(0), zeronames = character(0),
                 tcovar, covar, PACKAGE)
onestep.density(x, times, params, dens.fun, \dots,
                statenames = character(0), paramnames = character(0),
                covarnames = character(0), tcovar, covar, log = FALSE,
                PACKAGE)
}
\arguments{
  \item{xstart}{
    Matrix (dimensions \code{nvar} x \code{nrep}) of states at initial time \code{times[1]}.
  }
  \item{x}{
    Matrix (dimensions \code{nvar} x \code{nrep} x \code{ntimes}) of states at times \code{times}.
  }
  \item{times}{
    Vector of times (length \code{ntimes}) at which states are required or given.
  }
  \item{params}{
    Matrix containing parameters of the model.
    The \code{nrep} columns of \code{params} correspond to those of \code{xstart}.
  }
  \item{step.fun}{
    This can be either an R function or a compiled, dynamically loaded native function containing the model simulator.
    It should be written to take a single Euler step from a single point in state space.
    If it is a native function, it must be of type "pomp\_onestep\_sim" as defined in the header "pomp.h", which is included with the package.
    For details on how to write such codes, see Details.
  }
  \item{dens.fun}{
    This can be either an R function or a compiled, dynamically loaded native function containing the model transition log probability density function.
    This function will be called to compute the log likelihood of the actual Euler steps.
    It must be of type "pomp\_onestep\_pdf" as defined in the header "pomp.h", which is included with the package.
    For details on how to write such codes, see Details.
  }
  \item{delta.t}{
    Time interval of Euler steps.
  }
  \item{statenames}{
    Names of state variables, in the order they will be expected by the routine named in \code{step.fun} and \code{dens.fun}.
  }
  \item{paramnames}{
    Names of parameters, in the order they will be expected by the routine named in \code{step.fun} and \code{dens.fun}.
  }
  \item{covarnames}{
    Names of columns of the matrix of covariates \code{covar}, in the order they will be expected by the routine named in \code{step.fun} and \code{dens.fun}.
  }
  \item{zeronames}{
    Names of additional variables which will be zeroed before each time in \code{times}.
    These are useful, e.g., for storing accumulations of state variables.
  }
  \item{tcovar}{
    Times at which covariates are measured.
  }
  \item{covar}{
    Matrix of covariates.
    This should have dimensions \code{length(tcovar)} x \code{ncovar}, where \code{ncovar} is the number of covariates.
  }
  \item{log}{
    logical; if TRUE, probabilities p are given as log(p).
  }
  \item{\dots}{
    if \code{step.fun} (or \code{dens.fun}) is an R function, then additional arguments will be passed to it.
    If \code{step.fun} (or \code{dens.fun}) is a native routine, then additional arguments are ignored.
  }
  \item{PACKAGE}{
    an optional argument that specifies to which dynamically loaded library we restrict the search for the native routines.
    If this is '"base"', we search in the R executable itself.
  }
}
\details{
  \code{onestep.simulate} assumes that a single call to \code{step.fun} will advance the state process from one time to the next.
  \code{euler.simulate} will take multiple Euler steps, each of size at most \code{delta.t} (see below for information on how the actual Euler step size is chosen) to get from one time to the next.

  \code{onestep.density} assumes that no state transitions occure between consecutive times.

  If \code{step.fun} is written as an R function, it must have at least the arguments \code{x}, \code{t}, \code{params}, \code{delta.t}, and \code{\dots}.
  On a call to this function, \code{x} will be a named vector of state variables, \code{t} a scalar time, and \code{params} a named vector of parameters.
  The length of the Euler step will be \code{delta.t}.
  If the argument \code{covars} is included and a covariate table has been included in the \code{pomp} object, then on a call to this function, \code{covars} will be filled with the values, at time \code{t}, of the covariates.
  This is accomplished via interpolation of the covariate table.
  Additional arguments may be given: these will be filled by the correspondingly-named elements in the \code{userdata} slot of the \code{pomp} object (see \code{\link{pomp}}).

  If \code{step.fun} is written in a native language, it must be a function of type "pomp\_onestep\_sim" as specified in the header "pomp.h" included with the package (see the directory "include" in the installed package directory).

  If \code{dens.fun} is written as an R function, it must have at least the arguments \code{x1}, \code{x2}, \code{t1}, \code{t2}, \code{params}, and \code{\dots}.
  On a call to this function, \code{x1} and \code{x2} will be named vectors of state variables at times \code{t1} and \code{t2}, respectively.
  The named vector \code{params} contains the parameters.
  If the argument \code{covars} is included and a covariate table has been included in the \code{pomp} object, then on a call to this function, \code{covars} will be filled with the values, at time \code{t1}, of the covariates.
  This is accomplished via interpolation of the covariate table.
  As above, any additional arguments will be filled by the correspondingly-named elements in the \code{userdata} slot of the \code{pomp} object (see \code{\link{pomp}}).

  If \code{dens.fun} is written in a native language, it must be a function of type "pomp\_onestep\_pdf" as defined in the header "pomp.h" included with the package (see the directory "include" in the installed package directory).
}
\value{
  \code{euler.simulate} and \code{onestep.simulate} each return a \code{nvar} x \code{nrep} x \code{ntimes} array, where \code{nvar} is the number of state variables, \code{nrep} is the number of replicate simulations (= number of columns of \code{xstart} and \code{params}), and \code{ntimes} is the length of \code{times}.
  If \code{x} is this array, \code{x[,,1]} will be identical to \code{xstart}; the rownames of \code{x} and \code{xstart} will also coincide.

  \code{onestep.density} returns a \code{nrep} x \code{ntimes-1} array.
  If \code{f} is this array, \code{f[i,j]} is the likelihood of a transition from \code{x[,i,j]} to \code{x[,i,j+1]} in exactly one Euler step of duration \code{times[j+1]-times[j]}.
}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{eulermultinom}}, \code{\link{pomp}}}
\examples{
## an example showing how to use these functions to implement a seasonal SIR model is contained
## in the 'examples' directory
\dontrun{
edit(file=system.file("examples/euler_sir.R",package="pomp"))
}
}
\keyword{models}
