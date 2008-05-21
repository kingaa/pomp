\name{pomp-class}
\docType{class}
\alias{pomp-class}
\title{Partially-observed Markov process}
\description{The class \code{pomp} encodes a partially-observed Markov process.}
\section{Objects from the Class}{
  Objects should be created by calls of the function
  \code{pomp}.  See the documentation for \code{\link{pomp}} for
  usage instructions and important warnings.
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
    \item{rmeasure}{
      Function of prototype \code{rmeasure(x,times,params,\dots)} which simulates from the observation process.
    }
    \item{dmeasure}{
      Function of prototype \code{dmeasure(y,x,params,log=FALSE,\dots)} which gives the likelihood of \code{y} given \code{x}.
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
    \item{userdata}{
      A list containing any objects the user desires.
      Using this mechanism, the user can store additional information necessary for the definition of the model.
    }
  }
}
\section{Methods}{
  See the pomp methods documentation: \link{pomp-methods}.
}
\author{Aaron A. King (kingaa at umich dot edu)}
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
