\name{init.state-pomp}
\docType{methods}
\alias{init.state}
\alias{init.state,pomp-method}
\alias{init.state-pomp}

\title{
  Return a matrix of initial conditions given a vector of parameters and an initial time.
}
\description{
  The method \code{init.state} returns a vector of initial conditions for the state process when given a vector of parameters \code{params} and an initial time \code{t0}.
}
\usage{
init.state(object, params, t0, \dots)
\S4method{init.state}{pomp}(object, params, t0, \dots)
}
\arguments{
  \item{object}{an object of class \code{pomp}.}
  \item{params}{
    a named vector of parameters.
  }
  \item{t0}{the initial time at which initial states are requested.}
  \item{\dots}{at present, these are ignored.}
}
\value{
  Returns a matrix of initial states (with rownames).
}
\references{}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{pomp-class}}}
\keyword{models}
\keyword{ts}
