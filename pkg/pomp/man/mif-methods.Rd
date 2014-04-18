\name{mif-methods}
\docType{methods}
\alias{mif-methods}
\alias{logLik,mif-method}
\alias{logLik-mif}
\alias{conv.rec}
\alias{conv.rec,mif-method}
\alias{conv.rec-mif}
\alias{plot-mif}
\alias{plot,mif-method}
\alias{compare.mif}
\title{Methods of the "mif" class}
\description{Methods of the \code{mif} class.}
\usage{
\S4method{logLik}{mif}(object, \dots)
\S4method{conv.rec}{mif}(object, pars, transform = FALSE, \dots)
\S4method{plot}{mif}(x, y = NULL, \dots)
compare.mif(z)
}
\arguments{
  \item{object}{The \code{mif} object.}
  \item{pars}{Names of parameters.}
  \item{x}{The \code{mif} object.}
  \item{y}{Ignored.}
  \item{z}{A \code{mif} object or list of \code{mif} objects.}
  \item{transform}{
    optional logical;
    should the parameter transformations be applied?
    See \code{\link[=coef-pomp]{coef}} for details.
  }
  \item{\dots}{
    Further arguments (either ignored or passed to underlying functions).
  }
}
\section{Methods}{
  \describe{
    \item{conv.rec}{
      \code{conv.rec(object, pars = NULL)} returns the columns of the convergence-record matrix corresponding to the names in \code{pars}.
      By default, all rows are returned.
    }
    \item{logLik}{
      Returns the value in the \code{loglik} slot.
    }
    \item{compare.mif}{
      Given a \code{mif} object or a list of \code{mif} objects, \code{compare.mif} produces a set of diagnostic plots.
    }
    \item{plot}{
      Plots a series of diagnostic plots.
      When \code{x} is a \code{mif} object, \code{plot(x)} is equivalent to \code{compare.mif(list(x))}.
    }
  }
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{mif}}, \code{\link{pomp}}, \code{\link{pomp-class}}, \code{\link{pfilter}}}
\keyword{models}
\keyword{ts}
