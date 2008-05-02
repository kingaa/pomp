\name{B-splines}
\alias{bspline.basis}
\alias{periodic.bspline.basis}
\title{B-spline bases}
\description{
  These functions generate B-spline basis functions as lookup tables.
  \code{bspline.basis} gives a set of basis functions.
  \code{periodic.bspline.basis} gives a basis of periodic spline functions.
}
\usage{
bspline.basis(x, degree = 3, knots)
periodic.bspline.basis(x, nbasis, degree = 3, period = 1)
}
\arguments{
  \item{x}{Vector at which the spline functions are to be evaluated.}
  \item{degree}{Degree of requested B-splines.}
  \item{knots}{Vector of positions of the knots.}
  \item{nbasis}{The number of basis functions to return.}
  \item{period}{The period of the requested periodic B-splines.}
}
\value{
  \item{bspline.basis}{
    Returns a matrix with \code{length(x)} rows and \code{nbasis=length(knots)-degree-1} columns.
    Each column contains the values one of the \code{nbasis} spline functions.
  }
  \item{periodic.bspline.basis}{
    Returns a matrix with \code{length(x)} rows and \code{nbasis} columns.
  }
}
\author{Aaron A. King (kingaa at umich dot edu)}
\examples{
x <- seq(-0.2,1.2,by=0.01)
y <- bspline.basis(x,degree=7,seq(0,1,length=14))
matplot(x,y)

x <- seq(-1,2,by=0.01)
y <- periodic.bspline.basis(x,nbasis=5)
matplot(x,y)
}
\keyword{smooth}
