\name{spect.pomp-class}
\docType{class}
\alias{spect.pomp-class}
\alias{spect.matched.pomp-class}
\keyword{internal}
\title{The "spect.pomp" and "spect.matched.pomp" classes}
\description{
  The \code{spect.pomp} class holds a \code{\link{pomp}} object together with power spectrum estimates of the data and of some number of model simulations.
  It is created by a call to \code{\link{spect}}.
  See \code{\link{spect}} for usage.
}
\section{Objects from the Class}{
  Objects can be created by calls to the \code{\link{spect}} method or to the \code{\link{spect.match}} function on a \code{\link{pomp}} object.
}
\section{Slots}{
  The \code{spect.pomp} and \code{spect.matched.pomp} classes are derived from the \code{pomp} class and therefore have all the slots of \code{pomp}:
  see \code{\link{pomp-class}} for details.
  A full description of slots in a \code{spect.pomp} or \code{spect.matched.pomp} object follows.
  \describe{
    \item{seed}{the seed of the RNG used.}
    \item{kernel.width}{width parameter of the smoothing kernel used.}
    \item{transform}{transformation function used.}
    \item{detrend}{detrending option used.}
    \item{fail.value}{value used to replace non-finite values of the objective function.}
    \item{freq}{numeric vector of the frequencies at which the power spectrum is estimated.}
    \item{datspec, simspec}{estimated power spectra for data and simulations, respectively.}
    \item{pvals}{
      one-sided p-values:
      fraction of the simulated spectra that differ more from the mean simulated spectrum than does the data.
      The metric used is \eqn{L^2}{L2} distance.
    }
    \item{evals}{
      number of function and gradient evaluations by the optimizer.
      See \code{\link{optim}}.
    }
    \item{value}{Value of the objective function.}
    \item{convergence, msg}{
      Convergence code and message from the optimizer.
      See \code{\link{optim}}.
    }
    \item{data, times, t0, rprocess, dprocess, dmeasure, rmeasure,
      skeleton.type, skeleton, initializer, states, params,
      statenames, paramnames, covarnames, tcovar, covar,
      obsnames, statenames, paramnames, covarnames, PACKAGE,
      userdata, call}{
      Inherited from the \code{pomp} class.
    }
  }
}
\section{Extends}{
  Class \code{pomp}, directly.
  See \code{\link{pomp-class}}.
}
\section{Methods}{
  \describe{
    \item{plot}{
      plots the power spectrum of the data and a representation of the distribution of the simulated power spectra.
    }
  }
}
\author{
  Daniel C. Reuman \email{d.reuman at imperial dot ac dot uk}
  
  Aaron A. King \email{kingaa at umich dot edu}
}
\seealso{
  \code{\link{spect}},
  \code{\link{spect.match}},
  \code{\link{pomp}},
  \link{pomp-class}
}
\keyword{models}
\keyword{ts}
