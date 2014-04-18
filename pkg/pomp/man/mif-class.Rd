\name{mif-class}
\docType{class}
\alias{mif-class}
\keyword{internal}
\title{The "mif" class}
\description{
  The \code{mif} class holds a fitted model and is created by a call to \code{\link{mif}}.
  See \code{\link{mif}} for usage.
}
\section{Objects from the Class}{
  Objects can be created by calls to the \code{\link{mif}} method on a \code{\link{pomp}} object.
  Such a call uses the MIF algorithm to fit the model parameters.
}
\section{Slots}{
  A \code{mif} object is derived from a \code{pfilterd.pomp} object and therefore has all the slots of such an object.
  A description of \code{mif}-specific slots follows.
  \describe{
    \item{ivps}{
      character; contains the names of initial-value parameters (IVPs).
      These are parameters which are to be estimated using fixed-lag smoothing.
    }
    \item{pars}{
      character; contains the names of parameters to be estimated using MIF.
    }
    \item{transform}{
      logical scalar; was the estimation done on the transformed scale?
    }
    \item{Nmif}{
      integer; number of MIF iterations that have been completed.
    }
    \item{particles}{
      function of prototype \code{particles(Np,center,sd,...)} that draws particles from a distribution centered on \code{center} and with width proportional to \code{sd}.
      This function can be optionally specified by the user.
    }
    \item{var.factor}{
      the scaling coefficient relating the width of the initial particle distribution to \code{rw.sd}
    }
    \item{ic.lag}{
      the fixed lag used in the estimation of initial-value parameters (IVPs)
    }
    \item{cooling.type, cooling.fraction}{
      the cooling schedule specifications
    }
    \item{method}{
      character; the mif update method.
    }
    \item{random.walk.sd}{
      named vector containing the random-walk variance to be used for ordinary parameters.
    }
    \item{conv.rec}{
      the \dQuote{convergence record}: a matrix containing a record of the parameter values, log likelihoods, and other pertinent information, with one row for each MIF iteration.
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
  See \code{\link{mif}}, \link{mif-methods}, \link{particles-mif}
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{mif}}, \link{mif-methods}, \code{\link{pfilter}}, \code{\link{pomp}}, \link{pomp-class}}
\keyword{models}
\keyword{ts}
