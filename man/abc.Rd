\title{Estimation by approximate Bayesian computation (ABC)}
\name{Approximate Bayesian computation}
\docType{methods}
\alias{abc}
\alias{ABC}
\alias{abc-abc}
\alias{abc,probed.pomp-method}
\alias{abc,ANY-method}
\alias{abc,missing-method}
\alias{abc-probed.pomp}
\alias{abc,pomp-method}
\alias{abc-pomp}
\alias{abc,abc-method}
\alias{abc-abc}
\alias{continue,abc-method}
\alias{continue-abc}
\alias{abc-class}
\alias{abc-methods}
\alias{conv.rec,abc-method}
\alias{conv.rec-abc}
\alias{conv.rec,abcList-method}
\alias{conv.rec-abcList}
\alias{plot-abc}
\alias{plot,abc-method}
\alias{plot-abcList}
\alias{plot,abcList-method}
\alias{abcList-class}
\alias{covmat}
\alias{covmat,abc-method}
\alias{covmat-abc}
\alias{covmat,abcList-method}
\alias{covmat-abcList}
\alias{c-abc}
\alias{c,abc-method}
\alias{c-abcList}
\alias{c,abcList-method}
\alias{[-abcList}
\alias{[,abcList-method}
\description{
  The approximate Bayesian computation (ABC) algorithm for estimating the parameters of a partially-observed Markov process.
}
\usage{
\S4method{abc}{pomp}(object, Nabc = 1, start,
    proposal, probes, scale, epsilon,
    verbose = getOption("verbose"), \dots)
\S4method{abc}{probed.pomp}(object, probes,
    verbose = getOption("verbose"), \dots)
\S4method{abc}{abc}(object, Nabc, start, proposal,
    probes, scale, epsilon,
    verbose = getOption("verbose"), \dots)
\S4method{continue}{abc}(object, Nabc = 1, \dots)
\S4method{conv.rec}{abc}(object, pars, \dots)
\S4method{conv.rec}{abcList}(object, \dots)
\S4method{plot}{abc}(x, y, pars, scatter = FALSE, \dots)
\S4method{plot}{abcList}(x, y, \dots)
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
  \item{pars}{Names of parameters.}
  \item{scatter}{
    optional logical;
    If \code{TRUE}, draw scatterplots.
    If \code{FALSE}, draw traceplots.
  }
  \item{x}{\code{abc} object.}
  \item{y}{Ignored.}
  \item{\dots}{
    Additional arguments.
    These are currently ignored.
  }
}
\section{Running ABC}{
  \code{abc} returns an object of class \code{abc}.
  One or more \code{abc} objects can be joined to form an \code{abcList} object.
}
\section{Re-running ABC iterations}{
  To re-run a sequence of ABC iterations, one can use the \code{abc} method on a \code{abc} object.
  By default, the same parameters used for the original ABC run are re-used (except for \code{tol}, \code{max.fail}, and \code{verbose}, the defaults of which are shown above).
  If one does specify additional arguments, these will override the defaults.
}
\section{Continuing ABC iterations}{
  One can continue a series of ABC iterations from where one left off using the \code{continue} method.
  A call to \code{abc} to perform \code{Nabc=m} iterations followed by a call to \code{continue} to perform \code{Nabc=n} iterations will produce precisely the same effect as a single call to \code{abc} to perform \code{Nabc=m+n} iterations.
  By default, all the algorithmic parameters are the same as used in the original call to \code{abc}.
  Additional arguments will override the defaults.
}
\section{Methods}{
  Methods that can be used to manipulate, display, or extract information from an \code{abc} object:
  \describe{
    \item{\code{conv.rec(object, pars)} }{
      returns the columns of the convergence-record matrix corresponding to the names in \code{pars}.
      By default, all rows are returned.
    }
    \item{c}{
      Concatenates \code{abc} objects into an \code{abcList}.
    }
    \item{plot}{
      Diagnostic plots.
    }
    \item{\code{covmat(object, start, thin, expand)}}{
      computes the empirical covariance matrix of the ABC samples beginning with iteration \code{start} and thinning by factor \code{thin}.
      It expands this by a factor \code{expand^2/n}, where \code{n} is the number of parameters estimated.
      The intention is that the resulting matrix is a suitable input to the proposal function \code{\link{mvn.rw}}.
    }
  }
}
\references{
  T. Toni and M. P. H. Stumpf,
  Simulation-based model selection for dynamical systems in systems and population biology,
  Bioinformatics 26:104--110, 2010.

  T. Toni, D. Welch, N. Strelkowa, A. Ipsen, and M. P. H. Stumpf,
  Approximate Bayesian computation scheme for parameter inference and model selection in dynamical systems
  Journal of the Royal Society, Interface 6:187--202, 2009.
}
\author{Edward L. Ionides, Aaron A. King}
\seealso{
  \code{\link{pomp}}, \code{\link{probe}}, \link{MCMC proposal distributions}, and the tutorials on the \href{https://kingaa.github.io/pomp/}{package website}.  
}
\keyword{ts}
