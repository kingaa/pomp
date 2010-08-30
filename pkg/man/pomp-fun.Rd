\name{pomp-fun}
\docType{methods}
\alias{pomp.fun-methods}
\alias{show,pomp.fun-method}
\alias{show-pomp.fun}
\alias{print,pomp.fun-method}
\alias{print-pomp.fun}
\alias{pomp.fun-class}
\alias{pomp.fun}
\keyword{internal}
\title{Definition and methods of the "pomp.fun" class}
\description{Definition and methods of the "pomp.fun" class}
\usage{
pomp.fun(f = NULL, PACKAGE, proto = NULL)
\S4method{show}{pomp.fun}(object)
\S4method{print}{pomp.fun}(x, \dots)
}
\arguments{
  \item{f}{A function or the name of a native routine.}
  \item{PACKAGE}{
    optional; the name of the dynamically-loadable library in which the native function \code{f} can be found.
  }
  \item{proto}{
    optional string; a prototype against which \code{f} will be checked.
  }
  \item{object, x}{The \code{pomp.fun} object.}
}
\details{
  The sQuote{pomp.fun} class helps to settle common issues associated with user-defined functions which can be defined either via R code or by a native, compiled routine.
  It is not exported to userland.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{
  \code{\link{pomp}}
}
