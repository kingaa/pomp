\name{particles-mif}
\docType{methods}
\alias{particles}
\alias{particles,mif-method}
\alias{particles-mif}
\title{Generate particles from the user-specified distribution.}
\description{Generate particles from the user-specified distribution.}
\usage{
particles(object, \dots)
\S4method{particles}{mif}(object, Np = 1, center = coef(object), sd = 0, \dots)
}
\arguments{
  \item{object}{the "mif" object}
  \item{Np}{the number of particles, i.e., number of draws.}
  \item{center}{the central value of the distribution of particles}
  \item{sd}{the width of the distribution}
  \item{\dots}{additional arguments.  At present, these are ignored.}
}
\details{
  The \code{particles} method is used to set up the initial distribution
  of particles.  It is an interface to the user-specifed
  \code{particles} slot in the "mif" object.
}
\value{
  \code{particles} returns a list of two matrices.  \code{states}
  contains the state-variable portion of the particles; \code{params}
  contains the parameter portion.  Each has \code{Np} columns.
}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{mif}}, \link{mif-methods}, \code{\link{pomp}},
  \link{pomp-class}}
\keyword{models}
\keyword{ts}
