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
\usage{
mif(object, \dots)
\S4method{mif}{pomp}(object, Nmif, start, pars, ivps = character(0),
    particles, rw.sd, alg.pars, weighted = TRUE,
    tol = 1e-17, warn = TRUE, max.fail = 0, verbose = FALSE)
\S4method{mif}{mif}(object, \dots)
\S4method{continue}{mif}(object, Nmif, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp}.
  }
  \item{Nmif}{
    The number of MIF iterations to perform.
  }
  \item{start}{
    The initial guess of the parameters.
    This must be a named vector.
  }
  \item{pars}{
    character vector;
    names of ordinary parameters to be estimated.
  }
  \item{ivps}{
    character vector;
    names of initial-value parameters to be estimated.
  }
  \item{particles}{
    Function of prototype \code{particles(Np,center,sd,...)} which sets up the initial particle matrix by drawing a sample of size \code{Np} from the initial particle distribution centered at \code{center} and of width \code{sd}.
    If \code{particles} is not supplied by the user, the default behavior is to draw the particles from a multivariate normal distribution with mean \code{center} and standard deviation \code{sd}.
  }
  \item{rw.sd}{
    numeric vector with names; the intensity of the random walk to be applied to parameters.
    The random walk is only applied to parameters named in \code{pars}.
    The algorithm requires that the random walk be nontrivial.
    Thus, each element in \code{rw.sd[pars]} must be positive.
    \code{rw.sd} is also used to scale the initial-value parameters (via the \code{particles} function).
    Therefore, each element of \code{rw.sd[ivps]} must be positive.
  }
  \item{alg.pars}{
    A named list of algorithm parameters.
    This consists of
    \describe{
      \item{Np}{
	a positive integer;
	the number of particles to use in filtering
      }
      \item{ic.lag}{
	a positive integer;
	the timepoint for fixed-lag smoothing of initial-value parameters
      }
      \item{var.factor}{
	a positive number;
	the scaling coefficient relating the width of the initial particle distribution to \code{rw.sd}
      }
      \item{cooling.factor}{
	a positive number not greater than 1;
	the exponential cooling factor, \code{alpha}.
      }
    }
  }
  \item{weighted}{
    Should a weighted average be used?
    If \code{weighted=F}, the MIF update is not used;
    instead, an unweighed average of the filtering means is used for the update.
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
  \item{verbose}{
    logical; if TRUE, print progress reports.
  }
  \item{\dots}{
    Additional arguments that can be used to override the defaults.
  }
}
\section{Re-running MIF Iterations}{
  To re-run a sequence of MIF iterations, one can use the \code{mif} method on a \code{mif} object.
  The call sequence is \code{mif(object)}.
  By default, the same parameters used for the original MIF run are re-used.
  If one does specify additional arguments, these will override the defaults.
  An exception is that one cannot override the \code{particles} function.
}
\section{Continuing MIF Iterations}{
  One can continue a series of MIF iterations from where one left off.
  The call sequence is \code{continue(object, Nmif)}.
  This will perform \code{Nmif} additional MIF iterations on the \code{mif} object \code{object}.
  A call to \code{mif} to perform \code{Nmif=m} iterations followed by a call to \code{continue} to perform \code{Nmif=n} iterations will produce precisely the same effect as a single call to \code{mif} to perform \code{Nmif=m+n} iterations.
  Additional arguments are passed to \code{mif}.
  This feature can be used to change any of the parameters (except the \code{particles} function).
}
\section{Details}{
  \strong{It is the user's responsibility to ensure that, if the optional \code{particles} argument is given, that the \code{particles} function satisfies the following conditions:}

  \code{particles} has at least the following arguments:
  \code{Np}, \code{center}, \code{sd}, and \code{\dots}.
  \code{Np} should be assumed to be an integer; \code{center} and \code{sd} will be named vectors of the same length.
  Additional arguments may be specified;
  these will be filled with the elements of the \code{userdata} slot of the underlying \code{pomp} object (see \code{\link{pomp-class}}).

  \code{particles} returns a \code{length(center)} x \code{Np} matrix with rownames.
  Each column represents a distinct particle.
  The rownames are used by the algorithms (see \code{mif}, \code{pfilter}).

  The center of the particle distribution returned by \code{particles} should be \code{center}.
  The width of the particle distribution should vary monotonically with \code{sd}.
  In particular, when \code{sd=0}, the \code{particles} should return matrices with \code{Np} identical columns, each corresponding to the parameters specified in \code{center}.
}
\references{
  E. L. Ionides, C. Bret\'o, & A. A. King,
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
