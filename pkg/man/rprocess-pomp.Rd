\name{rprocess-pomp}
\docType{methods}
\alias{rprocess}
\alias{rprocess,pomp-method}
\alias{rprocess-pomp}
\keyword{internal}
\title{Simulate the process model of a partially-observed Markov process}
\description{
  The method \code{rprocess} runs simulations of the the process-model portion of partially-observed Markov process.
  This function is part of the low-level interface to \code{pomp} objects.
  This help page does not give instructions on the implementation of models: see \code{\link{pomp}} for instructions.
}
\usage{
rprocess(object, xstart, times, params, \dots)
\S4method{rprocess}{pomp}(object, xstart, times, params, offset, \dots)
}
\arguments{
  \item{object}{an object of class \code{pomp}.}
  \item{xstart}{
    a rank-2 array containing the starting state of the system.
    Columns of \code{xstart} correspond to states; rows to state variables.
    If there is more than one column of \code{xstart}, multiple independent simulations will be performed, one corresponding to each column.
    Note that in this case, \code{params} must have the same number of columns as \code{xstart}.
  }
  \item{times}{
    a numerical vector containing times.
    The first entry of \code{times} is the initial time (corresponding to \code{xstart}).
    Subsequent times are the additional times at which the state of the simulated processes are required.
  }
  \item{params}{
    a rank-2 array of parameters with the parameters corresponding to the columns of \code{xstart}.
  }
  \item{offset}{
    throw away the first \code{offset} times.
  }
  \item{\dots}{at present, these are ignored.}
}
\value{
  \code{rprocess} returns a rank-3 array with rownames.
  Suppose \code{x} is the array returned.
  Then \code{dim(x)=c(nvars,nreps,ntimes-offset)}, where \code{nvars} is the number of state variables (=\code{nrow(xstart)}), \code{nreps} is the number of independent realizations simulated (=\code{ncol(xstart)}), and \code{ntimes} is the length of the vector \code{times}.
  \code{x[,j,k]} is the value of the state process in the \code{j}-th realization at time \code{times[k+offset]}.
  The rownames of \code{x} must correspond to those of \code{xstart}. 
}
\details{
  This function is essentially a wrapper around the user-supplied \code{rprocess} slot of the \code{pomp} object.
  For specifications on writing such a function, see \code{\link{pomp}}.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{pomp-class}}, \code{\link{pomp}}}
\keyword{models}
\keyword{ts}
