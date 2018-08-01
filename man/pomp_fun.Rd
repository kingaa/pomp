\name{pomp-fun}
\docType{methods}
\alias{pomp_fun-methods}
\alias{show,pomp_fun-method}
\alias{show-pomp_fun}
\alias{print,pomp_fun-method}
\alias{print-pomp_fun}
\alias{pomp_fun-class}
\alias{pomp_fun}
\alias{pomp_fun,function-method}
\alias{pomp_fun,character-method}
\alias{pomp_fun,missing-method}
\alias{pomp_fun,pomp_fun-method}
\alias{pomp_fun,Csnippet-method}
\title{Definition and methods of the "pomp_fun" class}
\description{Definition and methods of the \sQuote{pomp_fun} class}
\usage{
\S4method{pomp_fun}{function}(f, proto = NULL, \dots)
\S4method{pomp_fun}{character}(f, PACKAGE = NULL,
      obsnames = character(0), statenames = character(0),
      paramnames = character(0), covarnames = character(0), \dots)
\S4method{pomp_fun}{Csnippet}(f, slotname = NULL, libname = NULL,
      obsnames = character(0), statenames = character(0),
      paramnames = character(0), covarnames = character(0), \dots)
\S4method{pomp_fun}{pomp_fun}(f, \dots)
\S4method{pomp_fun}{missing}(f, \dots)
\S4method{pomp_fun}{NULL}(f, \dots)
\S4method{show}{pomp_fun}(object)
\S4method{print}{pomp_fun}(x, \dots)
}
\arguments{
  \item{f}{A function or the name of a native routine.}
  \item{PACKAGE}{
    optional; the name of the dynamically-loadable library in which the native function \code{f} can be found.
  }
  \item{proto}{
    optional string; a prototype against which \code{f} will be checked.
  }
  \item{object, x}{The \sQuote{pomp_fun} object.}
}
\details{
  The \sQuote{pomp_fun} class helps to settle common issues associated with user-defined functions which can be defined either via R code or by a native, compiled routine.
  It is not exported to userland.
}
\author{Aaron A. King}
\seealso{
  \code{\link{pomp}}
}
\keyword{programming}
\keyword{internal}
