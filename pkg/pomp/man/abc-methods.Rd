\name{abc-methods}
\docType{methods}
\alias{abc-methods}
\alias{conv.rec,abc-method}
\alias{conv.rec-abc}
\alias{plot-abc}
\alias{plot,abc-method}
\title{Methods of the "abc" class}
\description{Methods of the "abc" class.}
\usage{
\S4method{conv.rec}{abc}(object, pars, \dots)
\S4method{plot}{abc}(x, y, pars, scatter = FALSE, \dots)
}
\arguments{
  \item{object, x}{The \code{abc} object.}
  \item{pars}{Names of parameters.}
  \item{y}{Ignored.}
  \item{scatter}{
    optional logical;
    If \code{TRUE}, draw scatterplots.
    If \code{FALSE}, draw traceplots.
  }
  \item{\dots}{
    Further arguments (either ignored or passed to underlying functions).
  }
}
\section{Methods}{
  \describe{
    \item{conv.rec}{
      \code{conv.rec(object, pars)} returns the columns of the convergence-record matrix corresponding to the names in \code{pars}.
      By default, all rows are returned.
    }
    \item{plot}{
      Plots a series of diagnostic plots.
    }
  }
}
\author{Edward L. Ionides \email{ionides at umich dot edu}, Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{abc}}, \code{\link{pomp}}}
\keyword{ts}
