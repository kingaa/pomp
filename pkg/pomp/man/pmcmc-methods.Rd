\name{pmcmc-methods}
\docType{methods}
\alias{pmcmc-methods}
\alias{logLik,pmcmc-method}
\alias{logLik-pmcmc}
\alias{conv.rec,pmcmc-method}
\alias{conv.rec-pmcmc}
\alias{conv.rec,pmcmcList-method}
\alias{conv.rec-pmcmcList}
\alias{plot-pmcmc}
\alias{plot,pmcmc-method}
\alias{plot-pmcmcList}
\alias{plot,pmcmcList-method}
\alias{pmcmcList-class}
\alias{c-pmcmc}
\alias{c,pmcmc-method}
\alias{c-pmcmcList}
\alias{c,pmcmcList-method}
\alias{[-pmcmcList}
\alias{[,pmcmcList-method}
\title{Methods of the "pmcmc" class}
\description{Methods of the "pmcmc" class.}
\usage{
\S4method{logLik}{pmcmc}(object, \dots)
\S4method{conv.rec}{pmcmc}(object, pars, \dots)
\S4method{conv.rec}{pmcmcList}(object, \dots)
\S4method{plot}{pmcmc}(x, y = NULL, \dots)
\S4method{plot}{pmcmcList}(x, y, \dots)
\S4method{c}{pmcmc}(x, \dots, recursive = FALSE)
\S4method{c}{pmcmcList}(x, \dots, recursive = FALSE)
}
\arguments{
  \item{object, x}{The \code{pmcmc} or \code{pmcmcList} object.}
  \item{pars}{Names of parameters.}
  \item{y, recursive}{Ignored.}
  \item{\dots}{
    Further arguments (either ignored or passed to underlying functions).
  }
}
\section{Methods}{
  \describe{
    \item{conv.rec}{
      \code{conv.rec(object, pars)} returns the columns of the convergence-record matrix corresponding to the names in \code{pars} as an object of class \code{\link[coda]{mcmc}} or \code{\link[coda]{mcmc.list}}.
    }
    \item{plot}{
      Diagnostic plots.
    }
    \item{logLik}{
      Returns the value in the \code{loglik} slot.
    }
    \item{c}{
      Concatenates \code{pmcmc} objects into a \code{pmcmcList}.
    }
  }
}
\references{
  C. Andrieu, A. Doucet and R. Holenstein,
  Particle Markov chain Monte Carlo methods, 
  J. Roy. Stat. Soc B, to appear, 2010.

  C. Andrieu and G.O. Roberts,
  The pseudo-marginal approach for efficient computation,
  Ann Stat 37:697-725, 2009.
}
\author{Edward L. Ionides \email{ionides at umich dot edu}, Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{pmcmc}}, \code{\link{pomp}}, \code{\link{pomp-class}}, \code{\link{pfilter}}}
\keyword{models}
\keyword{ts}
