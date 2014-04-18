\name{mif2d.pomp-class}
\docType{class}
\alias{mif2d.pomp-class}
\keyword{internal}
\title{The "mif2d.pomp" class}
\description{
  The \code{mif2d.pomp} class holds a fitted model and is created by a call to \code{\link{mif2}}.
  See \code{\link{mif2}} for usage.
}
\section{Objects from the Class}{
  Objects can be created by calls to the \code{\link{mif2}} method on a \code{\link{pomp}} or \code{\link[pomp]{pfilterd.pomp-class}} object.
  Such a call uses the MIF2 algorithm to fit the model parameters.
}
\section{Slots}{
  A \code{mif2d.pomp} object is derived from a \code{pfilterd.pomp} object and therefore has all the slots of such an object.
  The \code{mif2d.pomp} specific slots are described below.
  \describe{
    \item{transform}{
      logical scalar; was the estimation done on the transformed scale?
    }
    \item{Nmif}{
      integer; number of MIF2 iterations that have been completed.
    }
    \item{perturb.fn}{
      function.
    }
    \item{conv.rec}{
      The \dQuote{convergence record}: a matrix containing a record of the parameter values, log likelihoods, and other pertinent information, with one row for each MIF2 iteration.
    }
    \item{pred.mean, pred.var, filter.mean,
      eff.sample.size, cond.loglik, saved.states, saved.params, seed,
      Np, tol, nfail, loglik,
      data, times, t0, rprocess, dprocess, dmeasure, rmeasure,
      skeleton.type, skeleton, initializer, states, params,
      statenames, paramnames, covarnames, obsnames, tcovar, covar,
      PACKAGE, userdata}{
      Inherited from the \code{pfilterd.pomp} class.
    }
  }
}
\section{Extends}{
  Class \code{pfilterd.pomp}, directly.
  See \code{\link{pfilter}}.
}
\section{Methods}{
  See \code{\link{mif2}}, \link{mif2-methods}.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link[pomp]{mif}}, \link{mif2-methods}, \code{\link[pomp]{pomp}}, \link[pomp]{pomp-class}}
\keyword{models}
\keyword{ts}
