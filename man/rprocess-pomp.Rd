\name{rprocess-pomp}
\docType{methods}
\alias{rprocess}
\alias{rprocess,pomp-method}
\alias{rprocess-pomp}

\title{Simulate the process model of a partially-observed Markov process}
\description{The method \code{rprocess} runs simulations of the the
  process-model portion of partially-observed Markov process.}
\usage{
rprocess(object, xstart, times, params, \dots)
\S4method{rprocess}{pomp}(object, xstart, times, params, \dots)
}
\arguments{
  \item{object}{an object of class \code{pomp}.}
  \item{xstart}{a rank-2 array containing the starting state of the
    system.  Columns of \code{xstart} correspond to states; rows to
    state variables.  If there is more than one column of \code{xstart},
    multiple independent simulations will be performed, one
    corresponding to each column.  Note that in this case, \code{params}
    must have the same number of columns as \code{xstart}.}
  \item{times}{a numerical vector containing times.  The first entry of
    \code{times} is the initial time (corresponding to \code{xstart}).
    Subsequent times are the values of time for which the state of the
    simulated processes are required.}
  \item{params}{a rank-2 array of parameters with the parameters
    corresponding to the columns of \code{xstart}.}
  \item{\dots}{at present, these are ignored.}
}
\references{}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{pomp-class}}}
\keyword{models}
\keyword{ts}
