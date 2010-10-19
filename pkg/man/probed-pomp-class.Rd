\name{probed.pomp-class}
\docType{class}
\alias{probed.pomp-class}
\alias{probe.matched.pomp-class}
\keyword{internal}
\title{The "probed.pomp" and "probe.matched.pomp" classes}
\description{
  The \code{probed.pomp} class holds a \code{\link{pomp}} object together with a list of \dQuote{probes} and the results of applying them to the real and model-simulated data.
  It is created by a call to \code{\link{probe}}.
  See \code{\link{probe}} for usage.
}
\section{Objects from the Class}{
  Objects can be created by calls to the \code{\link{probe}} method on a \code{\link{pomp}} object.
}
\section{Slots}{
  The \code{probed.pomp} and \code{probe.matched.pomp} classes are derived from the \code{pomp} class and therefore have all the slots of \code{pomp}:
  see \code{\link{pomp-class}} for details.
  A full description of slots in a \code{probed.pomp} or \code{probe.matched.pomp} object follows.
  \describe{
    \item{probes}{list of the probes applied.}
    \item{seed}{the seed of the RNG used.}
    \item{datvals, simvals}{values of each of the probes applied to the real and simulated data, respectively.}
    \item{quantiles}{fraction of simulations with probe values less than the value of the probe of the data.}
    \item{pvals}{two-sided p-values: fraction of the \code{simvals} that deviate more extremely from the mean of the \code{simvals} than does \code{datavals}.}
    \item{synth.loglik}{the log synthetic likelihood (Wood 2010).  This is the likelihood assuming that the probes are multivariate-normally distributed.}
    \item{weights}{relative weights applied to the discrepancies in computing the probe-matching objective function.}
    \item{fail.value}{value to use to replace non-finite values of the objective function.}
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
      PACKAGE, userdata}{
      Inherited from the \code{pomp} class.
    }
  }
}
\section{Extends}{
  Class \code{pomp}, directly.
  See \code{\link{pomp-class}}.
}
\section{Methods}{
  See \code{\link{probe}}, \code{\link[=summary-probed.pomp]{summary}}, \code{\link[=plot-probed.pomp]{plot}}.
}
\references{
  B. E. Kendall, C. J. Briggs, W. M. Murdoch, P. Turchin, S. P. Ellner, E. McCauley, R. M. Nisbet, S. N. Wood
  Why do populations cycle? A synthesis of statistical and mechanistic modeling approaches,
  Ecology, 80:1789--1805, 1999.

  S. N. Wood
  Statistical inference for noisy nonlinear ecological dynamic systems,
  Nature, 466: 1102--1104, 2010.
}
\author{
  Daniel C. Reuman \email{d.reuman at imperial dot ac dot uk}
  
  Aaron A. King \email{kingaa at umich dot edu}
}
\seealso{\code{\link{probe}}, \link{probed.pomp-methods}, \code{\link{pomp}}, \link{pomp-class}}
\keyword{models}
\keyword{ts}
