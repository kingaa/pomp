\name{mif2-methods}
\docType{methods}
\alias{mif2-methods}
\alias{logLik,mif2d.pomp-method}
\alias{logLik-mif2d.pomp}
\alias{conv.rec,mif2d.pomp-method}
\alias{conv.rec-mif2d.pomp}
\alias{plot-mif2d.pomp}
\alias{plot,mif2d.pomp-method}
\title{Methods of the "mif2d.class" class}
\description{Methods of the \code{mif2d.class} class.}
\usage{
\S4method{logLik}{mif2d.pomp}(object, \dots)
\S4method{conv.rec}{mif2d.pomp}(object, pars, transform = FALSE, \dots)
\S4method{plot}{mif2d.pomp}(x, y = NULL, \dots)
}
\arguments{
  \item{object, x}{The \code{mif2d.pomp} object.}
  \item{pars}{Names of parameters.}
  \item{y}{Ignored.}
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
    \item{plot}{
      Plots a series of diagnostic plots.
    }
  }
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
\seealso{\code{\link{mif}}, \code{\link{pomp}}, \code{\link{pomp-class}}, \code{\link{pfilter}}}
\keyword{models}
\keyword{ts}
