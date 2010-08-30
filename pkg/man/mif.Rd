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
\S4method{mif}{pomp}(object, Nmif = 1, start, pars, ivps = character(0),
    particles, rw.sd, Np, ic.lag, var.factor, cooling.factor,
    weighted = TRUE, tol = 1e-17, max.fail = 0,
    verbose = getOption("verbose"))
\S4method{mif}{mif}(object, Nmif, \dots)
\S4method{continue}{mif}(object, Nmif = 1, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp}.
  }
  \item{Nmif}{
    The number of MIF iterations to perform.
  }
  \item{start}{
    named numerical vector;
    the starting guess of the parameters.
  }
  \item{pars}{
    optional character vector naming the ordinary parameters to be estimated.
    Every parameter named in \code{pars} must have a positive random-walk standard deviation specified in \code{rw.sd}.
    Leaving \code{pars} unspecified is equivalent to setting it equal to the names of all parameters with a positive value of \code{rw.sd} that are not \code{ivps}.
  }
  \item{ivps}{
    optional character vector naming the initial-value parameters to be estimated.
    Every parameter named in \code{ivps} must have a positive random-walk standard deviation specified in \code{rw.sd}.
  }
  \item{particles}{
    Function of prototype \code{particles(Np,center,sd,...)} which sets up the starting particle matrix by drawing a sample of size \code{Np} from the starting particle distribution centered at \code{center} and of width \code{sd}.
    If \code{particles} is not supplied by the user, the default behavior is to draw the particles from a multivariate normal distribution with mean \code{center} and standard deviation \code{sd}.
  }
  \item{rw.sd}{
    numeric vector with names; the intensity of the random walk to be applied to parameters.
    The random walk is only applied to parameters named in \code{pars} (i.e., not to those named in \code{ivps}).
    The algorithm requires that the random walk be nontrivial, so each element in \code{rw.sd[pars]} must be positive.
    \code{rw.sd} is also used to scale the initial-value parameters (via the \code{particles} function).
    Therefore, each element of \code{rw.sd[ivps]} must be positive.
    The following must be satisfied:
    \code{names(rw.sd)} must be a subset of \code{names(start)},
    \code{rw.sd} must be non-negative (zeros are simply ignored),
    the name of every positive element of \code{rw.sd} must be in either \code{pars} or \code{ivps}.
  }
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
    the scaling coefficient relating the width of the starting particle distribution to \code{rw.sd}.
    In particular, the width of the distribution of particles at the start of the first MIF iteration will be \code{random.walk.sd*var.factor}.
  }
  \item{cooling.factor}{
    a positive number not greater than 1;
    the exponential cooling factor, \code{alpha}.
  }
  \item{weighted}{
    logical; if TRUE, the MIF update (a weighted average) is used.
    If FALSE, the MIF update is not used;
    instead, an unweighed average of the filtering means is used for the update.
  }
  \item{tol}{
    numeric scalar; particles with log likelihood below \code{tol} are considered to be \dQuote{lost}.
    A filtering failure occurs when, at some time point, all particles are lost.
  }
  \item{max.fail}{
    integer; maximum number of filtering failures permitted.
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
  By default, the same parameters used for the original MIF run are re-used (except for \code{weighted}, \code{tol}, \code{max.fail}, and \code{verbose}, the defaults of which are shown above).
  If one does specify additional arguments, these will override the defaults.
}
\section{Continuing MIF Iterations}{
  One can continue a series of MIF iterations from where one left off using the \code{continue} method.
  A call to \code{mif} to perform \code{Nmif=m} iterations followed by a call to \code{continue} to perform \code{Nmif=n} iterations will produce precisely the same effect as a single call to \code{mif} to perform \code{Nmif=m+n} iterations.
  By default, all the algorithmic parameters are the same as used in the original call to \code{mif}.
  Additional arguments will override the defaults.
}
\section{Details}{
  If \code{particles} is not specified, the default behavior is to draw the particles from a multivariate normal distribution.
  \strong{It is the user's responsibility to ensure that, if the optional \code{particles} argument is given, that the \code{particles} function satisfies the following conditions:}

  \code{particles} has at least the following arguments:
  \code{Np}, \code{center}, \code{sd}, and \code{\dots}.
  \code{Np} may be assumed to be a positive integer;
  \code{center} and \code{sd} will be named vectors of the same length.
  Additional arguments may be specified;
  these will be filled with the elements of the \code{userdata} slot of the underlying \code{pomp} object (see \code{\link{pomp-class}}).

  \code{particles} returns a \code{length(center)} x \code{Np} matrix with rownames matching the names of \code{center} and \code{sd}.
  Each column represents a distinct particle.

  The center of the particle distribution returned by \code{particles} should be \code{center}.
  The width of the particle distribution should vary monotonically with \code{sd}.
  In particular, when \code{sd=0}, the \code{particles} should return matrices with \code{Np} identical columns, each given by the parameters specified in \code{center}.
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
\seealso{
  \code{\link{mif-class}}, \code{\link{mif-methods}}, \code{\link{pomp}}, \code{\link{pomp-class}}, \code{\link{pfilter}}.
  See the \dQuote{intro_to_pomp} vignette for an example.
}
\keyword{models}
\keyword{ts}
