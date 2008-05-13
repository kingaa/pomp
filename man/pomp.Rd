\name{pomp}
\alias{pomp}
\title{Partially-observed Markov process object.}
\description{
  Create a new \code{pomp} object.
}
\usage{
  pomp(data, times, t0, rprocess, dprocess, rmeasure, dmeasure,
       initializer, \dots)
}
\arguments{
  \item{data}{
    An array holding the data.
    This array is of dimensions \code{nobs} x \code{ntimes}, where \code{nobs} is the number of observed variables and \code{ntimes} is the number of times at which observations were made.
    It is also possible to specify \code{data} as a data-frame, in which case the argument \code{times} must be the name of the column of observation times.
    Note that if \code{data} is a data-frame, it will be coerced to an array with storage-mode `double'.
  }
  \item{times}{
    The times corresponding to the observations.
    \code{times} must be a strictly increasing numeric vector.
    If \code{data} is a data-frame, \code{times} should be the name of the column of observation times.
  }
  \item{t0}{
    The zero-time.
    This must be no later than the time of the first observation, \code{times[1]}.
  }
  \item{rprocess}{
    Function of prototype \code{rprocess(xstart,times,params,\dots)} which simulates from the unobserved process.
  }
  \item{dprocess}{
    Function of prototype \code{dprocess(x,times,params,log=FALSE,\dots)} which evaluates the likelihood of a sequence of consecutive state transitions.
  }
  \item{rmeasure}{
    Function of prototype \code{rmeasure(x,times,params,\dots)} which simulates from the observation process.
  }
  \item{dmeasure}{
    Function of prototype \code{dmeasure(y,x,times,params,log=FALSE,\dots)} which gives the likelihood of \code{y} given \code{x}.
  }
  \item{initializer}{
    Function of prototype \code{initializer(params,t0,\dots)} which yields initial conditions for the state process when given a vector, \code{params}, of parameters.
    By default, i.e., if it is unspecified when \code{pomp} is called, the initializer assumes any names in \code{params} that end in ".0" correspond to initial value parameters.
    These are simply copied over as initial conditions when \code{init.state} is called (see \code{\link{init.state-pomp}}.
    The names of the state variables are the same as the corresponding initial value parameters, but with the ".0" dropped.
  }
  \item{\dots}{
    Any additional arguments are stored in a slot `userdata' and are passed as arguments to each of the functions \code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, and \code{initializer} whenever they are evaluated.
    Using this mechanism, the user can store additional information necessary for the definition of the model.
  }
}
\details{
  \strong{It is the user's responsibility to ensure that the \code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, and \code{initializer} functions satisfy the following conditions:}
  
  \code{rprocess} must have at least the following arguments:
  \code{xstart}, \code{times}, \code{params}, and \code{\dots}.
  It can also take additional arguments.
  It is guaranteed that these will be filled with the corresponding elements the user has included as additional arguments in the construction of the \code{pomp} object.

  In calls to \code{rprocess}, \code{xstart} will be a rank-2 array (matrix) with rows corresponding to state variables and columns corresponding to independent realizations of the process.
  \code{params} will similarly be a rank-2 array with rows corresponding to parameters and columns corresponding to independent realizations.
  The columns of \code{params} are to be matched up with those of \code{xstart};
  in particular, they will agree in number.
  Both \code{xstart} and \code{params} must have rownames, which are available for use by the user.

  \code{rprocess} must return a rank-3 array with rownames.
  Suppose \code{x} is the array returned.
  Then \code{dim(x)=c(nvars,nreps,ntimes)}, where \code{nvars} is the number of state variables (=\code{nrow(xstart)}), \code{nreps} is the number of independent realizations simulated (=\code{ncol(xstart)}), and \code{ntimes} is the length of the vector \code{times}.
  \code{x[,j,k]} is the value of the state process in the \code{j}-th realization at time \code{times[k]}.
  In particular, \code{x[,,1]} must be identical to \code{xstart}.
  The rownames of \code{x} must correspond to those of \code{xstart}.

  \code{dprocess} must have at least the following arguments:
  \code{x}, \code{times}, \code{params}, \code{log}, and \code{\dots}.
  It may take additional arguments.
  It is guaranteed that these will be filled with the corresponding elements the user has included as additional arguments in the construction of the \code{pomp} object.

  In calls to \code{dprocess}, \code{x} will be an \code{nvars} x \code{nreps} x \code{ntimes} array, where these terms have the same meanings as above.
  \code{params} will be a rank-2 array with rows corresponding to individual parameters and columns corresponding to independent realizations.
  The columns of \code{params} are to be matched up with those of \code{x}; in particular, they will agree in number.
  Both \code{x} and \code{params} must have rownames, available for use by the user.

  \code{dprocess} must return a rank-2 array (matrix).
  Suppose \code{d} is the array returned.
  Then \code{dim(d)=c(nreps,ntimes-1)}.
  \code{d[j,k]} is the probability density of the transition from state \code{x[,j,k-1]} at time \code{times[k-1]} to state \code{x[,j,k]} at time \code{times[k]}.
  If \code{log=TRUE}, then the log of the p.d.f. is returned.
  \strong{In writing this function, you may assume that the transitions are consecutive.}

  \code{rmeasure} must have at least the arguments \code{x}, \code{times}, \code{params}, and \code{\dots}.
  It may take additional arguments, which will be filled with user-specified data as above.
  \code{x} must be a rank-3 array dimension \code{c(nvars,nreps,ntimes)}, where these variables have the same meanings as above.
  \code{times} is the corresponding set of times.
  \code{params} is a rank-2 array of dimension \code{c(npars,nreps)} as above.

  \code{rmeasure} must return a rank-3 array with rownames.
  If \code{y} is the returned array, then \code{dim(y)=c(nobs,nreps,ntimes)}, where \code{nobs} is the number of observable variables and \code{nreps}, \code{ntimes} agree with the corresponding dimensions of \code{x}.
  \code{y[,j,k]} must be the vector of observables in the \code{j}-th realization at time \code{times[k]}.

  \code{dmeasure} must have at least the arguments \code{y}, \code{x}, \code{times}, \code{params}, \code{log}, and \code{\dots}.
  \code{y} is a rank-2 array of observations (\code{nobs} x \code{ntimes});
  \code{x}, a rank-3 array of states;
  \code{params}, a rank-2 array containing parameters;
  and \code{times}, the corresponding observation times.
  It may take additional arguments which will be filled with user-specified data as above.
  \code{dmeasure} must return a rank-2 array of dimension \code{nreps} x \code{ntimes}.
  If \code{d} is the returned array, then \code{d[j,k]} is the p.d.f. of \code{y[,k]} given \code{x[,j,k]} at time \code{times[k]}.
  If \code{log=TRUE}, then the log of the p.d.f. is returned.

  \code{initializer} must have at least the arguments \code{params}, \code{t0}, and \code{\dots}.
  \code{params} is a named vector of parameters.
  \code{t0} will be the time at which initial conditions are desired.
  \code{initializer} must return a named vector of initial conditions.
}
\value{An object of class \code{pomp}.}
\examples{
## See the vignettes.
}
\section{Warning}{
  Some error checking is done, but complete error checking is impossible.  
  If the user-specified functions do not conform to the above specifications (see Details), then the results may be invalid.
}
\references{}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\link{pomp-class}}
\keyword{models}
\keyword{ts}
