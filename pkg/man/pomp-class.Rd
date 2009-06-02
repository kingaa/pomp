\name{pomp-class}
\docType{class}
\alias{pomp-class}
\title{Partially-observed Markov process}
\description{
  The class \code{pomp} encodes a partially-observed Markov process.
  This page documents the structure of the class:
  see the documentation for \code{\link{pomp}} for usage instructions.
}
\section{Objects from the Class}{
  Objects should be created by calls of the function \code{pomp}.
}
\section{Slots}{
  \describe{
    \item{data}{
      An array holding the data.
      This array is of dimensions \code{nobs} x \code{ntimes}, where \code{nobs} is the number of observed variables and \code{ntimes} is the number of times at which observations were made.
    }
    \item{times}{
      The times corresponding to the observations.
      \code{times} must be a strictly increasing numeric vector.
    }
    \item{t0}{The zero-time.}
    \item{rprocess}{
      Function of prototype \code{rprocess(xstart,times,params,\dots)} which simulates from the unobserved process.
    }
    \item{dprocess}{
      Function of prototype \code{dprocess(x,times,params,log=FALSE,\dots)} which evaluates the likelihood of a sequence of consecutive state transitions.
    }
    \item{dmeasure}{
      an object of class \dQuote{pomp.fun} which encodes the measurement model density.
    }
    \item{rmeasure}{
      an object of class \dQuote{pomp.fun} which encodes the measurement model simulator.
    }
    \item{skeleton.type}{
      a character variable specifying whether the deterministic skeleton is a map or a vectorfield.
    }
    \item{skeleton}{
      an object of class \dQuote{pomp.fun} which encodes the deterministic skeleton.
    }
    \item{initializer}{
      Function of prototype \code{initializer(params,t0,\dots)} which gives a vector of initial conditions when given a vector of parameters, \code{params}, and a time \code{t0}.
    }
    \item{states}{
      An array to hold a trajectory of the unobserved process.
    }
    \item{params}{
      A named numeric vector to hold model parameters.
    }
    \item{covar, tcovar}{
      An optional table of covariates.
    }
    \item{statenames, paramnames, covarnames}{
      Names state variables, parameters, and covariates, respectively.
      These are used to make the interface with native routines more robust.
    }
    \item{PACKAGE}{
      Character variable giving the name of the dynamically loadable library containing native routines.
    }
    \item{userdata}{
      A list containing any objects the user desires.
      Using this mechanism, the user can store additional information necessary for the definition of the model.
    }
  }
}
\section{Methods}{
  See the pomp methods documentation: \link{pomp-methods}.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{
  \code{\link{pomp}},
  \link{pomp-methods},
  \code{\link{rprocess}},
  \code{\link{dprocess}},
  \code{\link{rmeasure}},
  \code{\link{dmeasure}},
  \code{\link{init.state}},
  \link{simulate-pomp}
}
\keyword{models}
\keyword{ts}
