\name{pfilter-methods}
\alias{logLik,pfilterd.pomp-method}
\alias{logLik-pfilterd.pomp}
\alias{$,pfilterd.pomp-method}
\alias{$-pfilterd.pomp}
\alias{pred.mean}
\alias{pred.mean,pfilterd.pomp-method}
\alias{pred.mean-pfilterd.pomp}
\alias{pred.var}
\alias{pred.var,pfilterd.pomp-method}
\alias{pred.var-pfilterd.pomp}
\alias{filter.mean}
\alias{filter.mean,pfilterd.pomp-method}
\alias{filter.mean-pfilterd.pomp}
\title{Methods of the "pfilterd.pomp" class}
\description{Methods of the "pfilterd.pomp" class.}
\usage{
\S4method{logLik}{pfilterd.pomp}(object, \dots)
\S4method{pred.mean}{pfilterd.pomp}(object, pars, \dots)
\S4method{pred.var}{pfilterd.pomp}(object, pars, \dots)
\S4method{filter.mean}{pfilterd.pomp}(object, pars, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pfilterd.pomp} or inheriting class \code{pfilterd.pomp}.
  }
  \item{pars}{Names of parameters.}
  \item{\dots}{
    Additional arguments unused at present.
  }
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{pfilter}}, \link{pomp-class}}
\keyword{models}
\keyword{ts}
