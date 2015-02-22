\name{abc}
\docType{methods}
\alias{abc}
\alias{abc-abc}
\alias{abc,probed.pomp-method}
\alias{abc-probed.pomp}
\alias{abc,pomp-method}
\alias{abc-pomp}
\alias{abc,abc-method}
\alias{abc-abc}
\alias{continue,abc-method}
\alias{continue-abc}
\alias{abc-class}
\title{The ABC algorithm}
\description{
  The approximate Bayesian computation (ABC) algorithm for estimating the parameters of a partially-observed Markov process.
}
\usage{
\S4method{abc}{pomp}(object, Nabc = 1, start, proposal,
    pars, rw.sd, probes, scale, epsilon,
    verbose = getOption("verbose"), \dots)
\S4method{abc}{probed.pomp}(object, probes,
    verbose = getOption("verbose"), \dots)
\S4method{abc}{abc}(object, Nabc, start, proposal,
    probes, scale, epsilon,
    verbose = getOption("verbose"), \dots)
\S4method{continue}{abc}(object, Nabc = 1, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp}.
  }
  \item{Nabc}{
    The number of ABC iterations to perform.
  }
  \item{start}{
    named numeric vector;
    the starting guess of the parameters.
  }
  \item{proposal}{
    optional function that draws from the proposal distribution.
    Currently, the proposal distribution must be symmetric for proper inference:
    it is the user's responsibility to ensure that it is.
    Several functions that construct appropriate proposal function are provided:
    see \link{MCMC proposal functions} for more information.
  }
  \item{pars}{Deprecated and ignored.  Will be removed in a future release.}
  \item{rw.sd}{
    Deprecated.  Will be removed in a future release.
    Specifying \code{rw.sd} is equivalent to setting \code{proposal=mvn.diag.rw(rw.sd)}.
  }
  \item{probes}{
    List of probes (AKA summary statistics).
    See \code{\link{probe}} for details.
  }
  \item{scale}{
    named numeric vector of scales.
  }
  \item{epsilon}{
    ABC tolerance.
  }
  \item{verbose}{
    logical; if TRUE, print progress reports.
  }
  \item{\dots}{
    Additional arguments.
    These are currently ignored.
  }
}
\value{
  An object of class \code{abc}.
  This class inherits from class \code{\link[=probed.pomp-class]{probed.pomp}} and contains the following additional slots:
  \describe{
    \item{pars, Nabc, dprior, hyperparams, scale, epsilon}{
      These slots hold the values of the corresponding arguments of the call to \code{abc}.
    }
  }
}
\section{Re-running ABC Iterations}{
  To re-run a sequence of ABC iterations, one can use the \code{abc} method on a \code{abc} object.
  By default, the same parameters used for the original ABC run are re-used (except for \code{tol}, \code{max.fail}, and \code{verbose}, the defaults of which are shown above).
  If one does specify additional arguments, these will override the defaults.
}
\section{Continuing ABC Iterations}{
  One can continue a series of ABC iterations from where one left off using the \code{continue} method.
  A call to \code{abc} to perform \code{Nabc=m} iterations followed by a call to \code{continue} to perform \code{Nabc=n} iterations will produce precisely the same effect as a single call to \code{abc} to perform \code{Nabc=m+n} iterations.
  By default, all the algorithmic parameters are the same as used in the original call to \code{abc}.
  Additional arguments will override the defaults.
}
\section{Details}{
  TO APPEAR.
}
\references{
  T. Toni and M. P. H. Stumpf,
  Simulation-based model selection for dynamical systems in systems and population biology,
  Bioinformatics 26:104--110, 2010.

  T. Toni, D. Welch, N. Strelkowa, A. Ipsen, and M. P. H. Stumpf,
  Approximate Bayesian computation scheme for parameter inference and model selection in dynamical systems
  Journal of the Royal Society, Interface 6:187--202, 2009.
}
\author{Edward L. Ionides \email{ionides at umich dot edu}, Aaron A. King \email{kingaa at umich dot edu}}
\seealso{
  \link{abc-methods}, \code{\link{pomp}}, \code{\link{probe}}.
  See the \dQuote{intro_to_pomp} vignette for an example
}
\keyword{ts}
