\name{pomp-fun}
\docType{methods}
\alias{pomp.fun-methods}
\alias{show,pomp.fun-method}
\alias{show-pomp.fun}
\alias{print,pomp.fun-method}
\alias{print-pomp.fun}
\alias{pomp.fun-class}
\alias{pomp.fun}
\alias{pomp.fun,function-method}
\alias{pomp.fun,character-method}
\alias{pomp.fun,missing-method}
\alias{pomp.fun,pomp.fun-method}
\alias{pomp.fun,Csnippet-method}
\keyword{internal}
\title{Definition and methods of the "pomp.fun" class}
\description{Definition and methods of the \code{pomp.fun} class}
\usage{
\S4method{pomp.fun}{function}(f, proto = NULL, \dots)
\S4method{pomp.fun}{character}(f, PACKAGE = NULL, 
      obsnames = character(0), statenames = character(0),
      paramnames = character(0), covarnames = character(0), \dots)
\S4method{pomp.fun}{Csnippet}(f, slotname = NULL, libname = NULL,
      obsnames = character(0), statenames = character(0),
      paramnames = character(0), covarnames = character(0), \dots)
\S4method{pomp.fun}{pomp.fun}(f, \dots)
\S4method{pomp.fun}{missing}(f, \dots)
\S4method{pomp.fun}{NULL}(f, \dots)
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
