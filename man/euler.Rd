\name{euler}
\alias{euler}
\alias{euler.simulate}
\title{Dynamical models based on stochastic Euler algorithms}
\description{
  Facilities for simulating stochastic Markov processes based on the Euler algorithm.
}
\usage{
euler.simulate(xstart, times, params, euler.step.fun, delta.t, statenames, paramnames, zeronames, tcovar, covar, PACKAGE)
}
\arguments{
  \item{xstart}{
    Matrix of states at initial time \code{times[1]}.
  }
  \item{times}{
    Times at which states are required.
  }
  \item{params}{
    Matrix containing parameters of the model.
    The columns of \code{params} correspond to those of \code{xstart}.
  }
  \item{euler.step.fun}{
    Name of a compiled, loaded native function containing the model itself.
    This function will be called to take the actual Euler steps.
    It must be of type "estep" as defined in the header "euler.h", which is included with the package.
  }
  \item{delta.t}{
    Time interval of Euler steps.
  }
  \item{statenames}{
    Names of state variables, in the order they will be expected by the routine named in \code{euler.step.fun}.
  }
  \item{paramnames}{
    Names of parameters, in the order they will be expected by the routine named in \code{euler.step.fun}.
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
  \item{PACKAGE}{
    an optional argument that specifies to which dynamically loaded library we restrict the search for \code{euler.step.fun}.
    If this is '"base"', we search in the R executable itself.
  }
}
\value{
  Returns a \code{nvar} x \code{nsim} x \code{ntimes} array, where \code{nvar} is the number of state variables, \code{nsim} is the number of replicate simulations (= number of columns of \code{xstart} and \code{params}), and \code{ntimes} is the length of \code{times}.
  If \code{x} is this array, \code{x[,,1]} will be identical to \code{xstart}; the rownames of \code{x} and \code{xstart} will also coincide.
}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{eulermultinom}}, \code{\link{pomp}}}
\examples{}
\keyword{smooth}
