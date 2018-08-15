\title{Estimation by approximate Bayesian computation (ABC)}
\name{Approximate Bayesian computation}
\docType{methods}
\alias{abc}
\alias{ABC}
\alias{abc,probed_pomp-method}
\alias{abc,ANY-method}
\alias{abc,missing-method}
\alias{abc-probed_pomp}
\alias{abc,pomp-method}
\alias{abc-pomp}
\alias{abc,abcd_pomp-method}
\alias{abc-abcd_pomp}
\alias{continue,abcd_pomp-method}
\alias{continue-abcd_pomp}
\alias{abc-methods}
\alias{traces,abcd_pomp-method}
\alias{traces-abcd_pomp}
\alias{traces,abcList-method}
\alias{traces-abcList}
\alias{plot-abcd_pomp}
\alias{plot,abcd_pomp-method}
\alias{plot-Abc}
\alias{plot,Abc-method}
\alias{c-Abc}
\alias{c,Abc-method}
\alias{covmat}
\alias{covmat,abcd_pomp-method}
\alias{covmat-abcd_pomp}
\alias{covmat,abcList-method}
\alias{covmat-abcList}
\description{
  The approximate Bayesian computation (ABC) algorithm for estimating the parameters of a partially-observed Markov process.
}
\usage{
\S4method{abc}{pomp}(object, Nabc = 1, start,
    proposal, probes, scale, epsilon,
    verbose = getOption("verbose"), \dots)
\S4method{abc}{probed_pomp}(object, probes,
    verbose = getOption("verbose"), \dots)
\S4method{abc}{abcd_pomp}(object, Nabc, start, proposal,
    probes, scale, epsilon,
    verbose = getOption("verbose"), \dots)
\S4method{continue}{abcd_pomp}(object, Nabc = 1, \dots)
\S4method{traces}{abcd_pomp}(object, pars, \dots)
\S4method{traces}{abcList}(object, \dots)
}
\arguments{
  \item{object}{
    An object of class \sQuote{pomp}.
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
  \item{\dots}{
    Additional arguments are passed to \code{\link{pomp}}, allowing one to supply new or modify existing model characteristics or components.
  }
}
\section{Running ABC}{
  \code{abc} returns an object of class \sQuote{abcd_pomp}.
  One or more \sQuote{abcd_pomp} objects can be joined to form an \sQuote{abcList} object.
}
\section{Re-running ABC iterations}{
  To re-run a sequence of ABC iterations, one can use the \code{abc} method on a \sQuote{abcd_pomp} object.
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
  Methods that can be used to manipulate, display, or extract information from an \sQuote{abc} object:
  \describe{
    \item{\code{traces(object, pars)} }{
      returns the history of the ABC iterations as an object of class \sQuote{\link[coda]{mcmc}} or \sQuote{\link[coda]{mcmc.list}}.
      Thus, it can be used as input for the various \pkg{coda} convergence diagnostics.
    }
    \item{c}{
      Concatenates \sQuote{abcd_pomp} objects into an \code{abcList}.
    }
    \item{\code{plot(x, pars, scatter)}}{
       displays diagnostic plots of the traces of the parameters named in \code{pars} either against ABC iteration number (the default, \code{scatter=FALSE}) or as a scatterplot (\code{scatter = TRUE}).
    }
    \item{\code{covmat(object, start, thin, expand)}}{
      computes the empirical covariance matrix of the ABC samples beginning with iteration \code{start} and thinning by factor \code{thin}.
      It expands this by a factor \code{expand^2/n}, where \code{n} is the number of parameters estimated.
      The intention is that the resulting matrix is a suitable input to the proposal function \code{\link{mvn.rw}}.
    }
  }
}
\references{
  J.-M. Marin, P. Pudlo, C. P. Robert, and R. J. Ryder,
  Approximate Bayesian computational methods.
  Statistics and Compuing 22:1167--1180, 2012.

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
