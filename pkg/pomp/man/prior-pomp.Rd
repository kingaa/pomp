\name{prior-pomp}
\docType{methods}
\alias{dprior}
\alias{dprior,pomp-method}
\alias{dprior-pomp}
\alias{rprior}
\alias{rprior,pomp-method}
\alias{rprior-pomp}
\keyword{internal}
\title{Evaluate or simulate from the prior probability density}
\description{
  The method \code{dprior} evaluates the prior probability density and \code{rprior} simulates from it.
  This function is part of the low-level interface to \code{pomp} objects.
  This help page does not give instructions on the implementation of models: see \code{\link{pomp}} for instructions.
}
\usage{
\S4method{dprior}{pomp}(object, params, log = FALSE, \dots)
\S4method{rprior}{pomp}(object, params, \dots)
}
\arguments{
  \item{object}{an object of class \code{pomp}.}
  \item{params}{
    a rank-2 array of parameters with columns corresponding to the columns of \code{x}.
  }
  \item{log}{if TRUE, log probabilities are returned.}
  \item{\dots}{at present, these are ignored.}
}
\details{
  These functions are essentially wrappers around the user-supplied \code{dprior} and \code{rprior} functions stored in the \code{pomp} object.
  For specifications on writing these functions, see \code{\link{pomp}}.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{pomp-class}}, \code{\link{pomp}}}
\keyword{models}
\keyword{ts}
