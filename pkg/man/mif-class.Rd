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
  Objects can be created by calls to the \code{\link{mif}} method on an \code{\link{pomp}} object.
  Such a call uses the MIF algorithm to fit the model parameters.
}
\section{Slots}{
  A \code{mif} object is derived from a \code{pomp} object and therefore has all the slots of such an object.
  See \code{\link{pomp-class}} for details.
  A full description of slots in a \code{mif} object follows.
  \describe{
    \item{ivps}{
      A character vector containing the names of initial-value parameters (IVPs).
      These are parameters which are to be estimated using fixed-lag smoothing.
    }
    \item{pars}{
      A character vector containing the names of parameters to be estimated using MIF.
    }
    \item{Nmif}{
      Number of MIF iterations that have been completed.
    }
    \item{particles}{
      A function of prototype \code{particles(Np,center,sd,...)} that draws particles from a distribution centered on \code{center} and with width proportional to \code{sd}.
      This function can be optionally specified by the user.
    }
    \item{alg.pars}{
      A named list of algorithm parameters.
      This consists of
      \code{Np}, the number of particles to use in filtering; 
      \code{var.factor}, the scaling coefficient relating the width of the initial particle distribution to \code{rw.sd};
      \code{ic.lag}, the fixed lag used in the estimation of initial-value parameters (IVPs);
      and \code{cooling.factor}, the exponential cooling factor, where \code{0<cooling.factor<1}.
    }
    \item{random.walk.sd}{
      A named vector containing the random-walk variance to be used for ordinary parameters.
    }
    \item{pred.mean}{
      Matrix of prediction means.
      See \code{\link{pfilter}}.
    }
    \item{pred.var}{
      Matrix of prediction variances.
      See \code{\link{pfilter}}.
    }
    \item{filter.mean}{
      Matrix of filtering means.
      See \code{\link{pfilter}}.
    }
    \item{eff.sample.size}{
      A vector containing the effective number of  particles at each time point.
      See \code{\link{pfilter}}.
    }
    \item{cond.loglik}{
      A vector containing the conditional log likelihoods at each time point.
      See \code{\link{pfilter}}.
    }
    \item{conv.rec}{
      The \dQuote{convergence record}: a matrix containing a record of the parameter values, log likelihoods, and other pertinent information, with one row for each MIF iteration.
    }
    \item{loglik}{
      A numeric value containing the value of the log likelihood, as evaluated for the random-parameter model.
      Note that this will not be equal to the log likelihood for the fixed-parameter model.
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
  See \code{\link{mif}}, \link{mif-methods}, \link{particles-mif}, \link{pfilter-mif}.
}
\references{
  E. L. Ionides, C. Bret\\'o, & A. A. King,
  Inference for nonlinear dynamical systems,
  Proc. Natl. Acad. Sci. U.S.A., 103:18438--18443, 2006.

  A. A. King, E. L. Ionides, M. Pascual, and M. J. Bouma,
  Inapparent infections and cholera dynamics,
  Nature, 454:877--880, 2008.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{mif}}, \link{mif-methods}, \code{\link{pomp}}, \link{pomp-class}}
\keyword{models}
\keyword{ts}
