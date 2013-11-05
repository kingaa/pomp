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
\alias{eff.sample.size}
\alias{eff.sample.size,pfilterd.pomp-method}
\alias{eff.sample.size-pfilterd.pomp}
\alias{cond.logLik}
\alias{cond.logLik,pfilterd.pomp-method}
\alias{cond.logLik-pfilterd.pomp}
\alias{as,pfilterd.pomp-method}
\alias{as.data.frame.pfilterd.pomp}
\alias{coerce,pfilterd.pomp,data.frame-method}
\title{Methods of the "pfilterd.pomp" class}
\description{Methods of the "pfilterd.pomp" class.}
\usage{
\S4method{logLik}{pfilterd.pomp}(object, \dots)
\S4method{pred.mean}{pfilterd.pomp}(object, pars, \dots)
\S4method{pred.var}{pfilterd.pomp}(object, pars, \dots)
\S4method{filter.mean}{pfilterd.pomp}(object, pars, \dots)
\S4method{eff.sample.size}{pfilterd.pomp}(object, \dots)
\S4method{cond.logLik}{pfilterd.pomp}(object, \dots)
\S4method{as}{pfilterd.pomp}(object, class)
\S4method{coerce}{pfilterd.pomp,data.frame}(from, to = "data.frame", strict = TRUE)
\S3method{as.data.frame}{pfilterd.pomp}(x, row.names, optional, \dots)
}
\arguments{
  \item{object, x}{
    An object of class \code{pfilterd.pomp} or inheriting class \code{pfilterd.pomp}.
  }
  \item{pars}{Names of parameters.}
  \item{class}{
    character;
    name of the class to which \code{object} should be coerced.
  }
  \item{from, to}{
    the classes between which coercion should be performed.
  }
  \item{strict}{
    ignored.
  }
  \item{row.names, optional}{
    ignored.
  }
  \item{\dots}{
    Additional arguments unused at present.
  }
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{pfilter}}, \link{pomp-class}}
\keyword{models}
\keyword{ts}
