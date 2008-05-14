\name{mif}
\docType{methods}
\alias{mif}
\alias{mif,mif-method}
\alias{mif-mif}
\alias{mif,pomp-method}
\alias{mif-pomp}
\alias{continue}
\alias{continue,mif-method}
\alias{continue-mif}
\title{The MIF algorithm}
\description{
  The MIF algorithm for estimating the parameters of a partially-observed Markov process.
}
\section{Running the MIF algorithm}{
  The call sequence for \code{mif} is

  \code{mif(object, Nmif, start, pars, ivps = character(0),
    particles, rw.sd, alg.pars, weighted = TRUE,
    tol = 1e-17, warn = TRUE, max.fail = 0)}

  Description of arguments:
  \describe{
    \item{object}{An object of class \code{pomp}.}
    \item{Nmif}{The number of MIF iterations to perform.}
    \item{start}{The initial guess of the parameters.  This must be a named vector.}
    \item{pars}{Character vector of names of ordinary parameters to be estimated.}
    \item{ivps}{Character vector of names of initial-value parameters to be estimated.}
    \item{particles}{
      Function of prototype \code{particles(Np,center,sd,...)} which sets up the initial particle matrix by drawing a sample of size \code{Np} from the initial particle distribution centered at \code{center} and of width \code{sd}.
      If \code{particles} is not supplied by the user, the default behavior is to draw the particles from a multivariate normal distribution with mean \code{center} and standard deviation \code{sd}.
    }
    \item{rw.sd}{The intensity of the random walk.}
    \item{alg.pars}{
      A named list of algorithm parameters.
      This consists of
      \item{Np}{the number of particles to use in filtering}
      \item{var.factor}{the scaling coefficient relating the width of the initial particle distribution to \code{rw.sd}}
      \item{ic.lag}{the timepoint for fixed-lag smoothing of initial-value parameters (IVPs)}
      \item{cooling.factor}{the exponential cooling factor, \eqn{a}, where \eqn{0<a<1}}
    }
    \item{weighted}{
      Should a weighted average be used?
      If \code{weighted=F}, the MIF update is not used;
      rather, an unweighed average of the filtering means is used for the update.
    }
    \item{tol}{
      Particles with log likelihood below \code{tol} are considered to be "lost".
      A filtering failure occurs when, at some time point, all particles are lost.
    }
    \item{warn}{
      Should a warning be generated when a filtering failure occurs?
    }
    \item{max.fail}{
      Maximum number of filtering failures permitted.
      If the number of failures exceeds this number, execution will terminate with an error.
    }
    }
}
\section{Re-running MIF Iterations}{
  To re-run a sequence of MIF iterations, one can use the \code{mif} method on a "mif" object.
  The call sequence is

  \code{mif(object, \dots)}

  Any additional arguments that are valid for the \code{mif} method of a "pomp" object (see \code{\link{pomp-class}}) can be given (with the exception of the \code{particles} function.
  Arguments not specified will take the values they have in the slots of \code{object}.
}
\section{Continuing MIF Iterations}{
  One can continue a series of MIF iterations from where one left off.
  The call sequence is

  \code{continue(object, Nmif, \dots)}

  This will perform \code{Nmif} additional MIF iterations on the "mif" object \code{object}.
  A call to \code{mif} to perform \code{Nmif=m} iterations followed by a call to \code{continue} to perform \code{Nmif=n} iterations will produce precisely the same effect as a single call to \code{mif} to perform \code{Nmif=m+n} iterations.

  All additional arguments are passed to \code{mif}.
  This feature can be used to change any of the parameters (except the \code{particles} function).
}
\section{Details}{
  \strong{It is the user's responsibility to ensure that the \code{particles} function satisfies the following conditions:}

  \code{particles} has at least the following arguments:
  \code{Np}, \code{center}, \code{sd}, and \code{\dots}.
  \code{Np} should be assumed to be an integer; \code{center} and \code{sd} will be named vectors of the same length.
  Additional arguments may be specified;
  these will be filled with the elements of the \code{userdata} slot of the underlying "pomp" object (see \code{\link{pomp-class}}).

  \code{particles} returns a \code{length(center)} x \code{Np} matrix with rownames.
  Each column represents a distinct particle.
  The rownames are used by the algorithms (see \code{mif}, \code{pfilter}).

  The center of the particle distribution returned by \code{particles} should be \code{center}.
  The width of the particle distribution should vary monotonically with \code{sd}.
  In particular, when \code{sd=0}, the \code{particles} should return matrices with \code{Np} identical columns, each corresponding to the parameters specified in \code{center}.
}
\references{
  E. L. Ionides, C. Bret{\'o}, & A. A. King,
  Inference for nonlinear dynamical systems,
  Proc. Natl. Acad. Sci. U.S.A., 103:18438--18443, 2006.
}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{
  \code{\link{mif-class}}, \code{\link{mif-methods}}, \code{\link{pomp}}, \code{\link{pomp-class}}, \code{\link{pfilter}}.
  See the "intro\_to\_pomp" vignette for an example.
}
\keyword{models}
\keyword{ts}
