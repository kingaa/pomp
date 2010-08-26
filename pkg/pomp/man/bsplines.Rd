\name{B-splines}
\alias{bspline.basis}
\alias{periodic.bspline.basis}
\title{B-spline bases}
\description{
  These functions generate B-spline basis functions.
  \code{bspline.basis} gives a basis of spline functions.
  \code{periodic.bspline.basis} gives a basis of periodic spline functions.
}
\usage{
bspline.basis(x, nbasis, degree = 3)
periodic.bspline.basis(x, nbasis, degree = 3, period = 1)
}
\arguments{
  \item{x}{Vector at which the spline functions are to be evaluated.}
  \item{nbasis}{The number of basis functions to return.}
  \item{degree}{Degree of requested B-splines.}
  \item{period}{The period of the requested periodic B-splines.}
}
\value{
  \item{bspline.basis}{
    Returns a matrix with \code{length(x)} rows and \code{nbasis} columns.
    Each column contains the values one of the spline basis functions.
  }
  \item{periodic.bspline.basis}{
    Returns a matrix with \code{length(x)} rows and \code{nbasis} columns.
    The basis functions returned are periodic with period \code{period}.
  }
}
\details{
  Direct access to the underlying C routines is available.
  See the header file \dQuote{pomp.h} for details.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
x <- seq(0,2,by=0.01)
y <- bspline.basis(x,degree=3,nbasis=9)
matplot(x,y,type='l',ylim=c(0,1.1))
lines(x,apply(y,1,sum),lwd=2)

x <- seq(-1,2,by=0.01)
y <- periodic.bspline.basis(x,nbasis=5)
matplot(x,y,type='l')
}
\keyword{smooth}
