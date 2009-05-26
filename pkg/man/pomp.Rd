\name{pomp}
\alias{pomp}
\title{Partially-observed Markov process object.}
\description{
  Create a new \code{pomp} object.
}
\usage{
  pomp(data, times, t0, \dots, rprocess, dprocess, rmeasure, dmeasure,
       measurement.model, skeleton.map, skeleton.vectorfield, initializer,
       covar, tcovar, statenames, paramnames, covarnames, PACKAGE)
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
    optional; a function of prototype \code{rprocess(xstart,times,params,\dots)} which simulates from the unobserved process.
  }
  \item{dprocess}{
    optional; a function of prototype \code{dprocess(x,times,params,log,\dots)} which evaluates the likelihood of a sequence of consecutive state transitions.
  }
  \item{rmeasure}{
    optional; the measurement model simulator.
    This can be specified in one of three ways:
    (1) as a function of prototype \code{rmeasure(x,t,params,\dots)} which makes a draw from the observation process given states \code{x}, time \code{t}, and parameters \code{params}.
    (2) as the name of a native (compiled) routine with prototype "pomp\_measure\_model\_simulator" as defined in the header file "pomp.h".
    In the above cases, if the measurement model depends on covariates, the optional argument \code{covars} will be filled with interpolated values at each call.
    (3) using the formula-based \code{measurement.model} facility (see below).
  }
  \item{dmeasure}{
    optional; the measurement model probability density function.
    This can be specified in one of three ways:
    (1) as a function of prototype \code{dmeasure(y,x,t,params,log,\dots)} which computes the p.d.f. of \code{y} given \code{x}, \code{t}, and \code{params}.
    (2) as the name of a native (compiled) routine with prototype "pomp\_measure\_model\_density" as defined in the header file "pomp.h".
    In the above cases, if the measurement model depends on covariates, the optional argument \code{covars} will be filled with interpolated values at each call.
    (3) using the formula-based \code{measurement.model} facility (see below).
    As might be expected, if \code{log=TRUE}, this function should return the log likelihood.
  }
  \item{measurement.model}{
    optional; a formula or list of formulae, specifying the measurement model.
    These formulae are parsed and used to generate \code{rmeasure} and \code{dmeasure} functions.
    If \code{measurement.model} is given it overrides any specification of \code{rmeasure} or \code{dmeasure}.
    See below for an example.
    NB: it will typically be possible to acclerate measurement model computations by writing \code{dmeasure} and/or \code{rmeasure} functions directly.
  }
  \item{skeleton.map}{
    optional.
    If we are dealing with a discrete-time Markov process, its deterministic skeleton is a map and can be specified using \code{skeleton.map} in one of two ways:
    (1) as a function of prototype \code{skeleton(x,t,params,\dots)} which evaluates the deterministic skeleton (whether vectorfield or map) of at state \code{x} and time \code{t} given the parameters \code{params}, or
    (2) as the name of a native (compiled) routine with prototype "pomp\_vectorfield\_map" as defined in the header file "pomp.h".
    If the deterministic skeleton depends on covariates, the optional argument \code{covars} will be filled with interpolated values of the covariates at the time \code{t}.
  }
  \item{skeleton.vectorfield}{
    optional.
    If we are dealing with a continuous-time Markov process, its deterministic skeleton is a vectorfield and can be specified using \code{skeleton.vectorfield} in either of the two ways described above, under \code{skeleton.map}.
    Note that it makes no sense to specify the deterministic skeleton both as a map and as a vectorfield: an attempt to do so will generate an error.
  }
  \item{initializer}{
    optional; a function of prototype \code{initializer(params,t0,\dots)} which yields initial conditions for the state process when given a vector, \code{params}, of parameters.
    By default, i.e., if it is unspecified when \code{pomp} is called, the initializer assumes any names in \code{params} that end in ".0" correspond to initial value parameters.
    These are simply copied over as initial conditions when \code{init.state} is called (see \code{\link{init.state-pomp}}).
    The names of the state variables are the same as the corresponding initial value parameters, but with the ".0" dropped.
  }
  \item{covar, tcovar}{
    An optional table of covariates: \code{covar} is the table (with one column per variable) and \code{tcovar} the corresponding times (one entry per row of \code{covar}).
    This can be in the form of a matrix or a data frame.
    In either case the columns are taken to be distinct covariates.
    If \code{covar} is a data frame, \code{tcovar} can be either the name or the index of the time variable.
    If a covariate table is supplied, then the value of each of the covariates is interpolated each time that \code{rmeasure}, \code{dmeasure}, or \code{skeleton} is evaluated.
    The resulting interpolated values are passed to the corresponding functions as a numeric vector named \code{covars}.
  }
  \item{statenames, paramnames, covarnames}{
    Optional character vectors specifying the names of state variables, parameters, or covariates, respectively.
    These are only used in the event that one or more of the basic functions (\code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, \code{skeleton}) are defined using native routines.
    In that case, these name vectors are matched against the corresponding names and the indices of the names are passed to the native routines.
  }
  \item{PACKAGE}{
    An optional string giving the name of the dynamically loaded library in which the native routines are to be found.
  }
  \item{\dots}{
    Any additional arguments are passed as arguments to each of the functions \code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, and \code{initializer} whenever they are evaluated.
  }
}
\details{
  \strong{
    It is not typically necessary (or easy) to define all of the functions \code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, and \code{initializer} in any given problem.
    Each algorithm makes use of a different subset of these functions.
    If one of these functions is defined, it is the user's responsibility to ensure that it satisfies the following conditions:
  }
  \describe{
    \item{\code{rprocess}}{
      if provided, must have at least the following arguments:
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
    }
    \item{\code{dprocess}}{
      if provided, must have at least the following arguments:
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
    }
    \item{\code{rmeasure}}{
      if provided, must take at least the arguments \code{x}, \code{t}, \code{params}, and \code{\dots}.
      It may take additional arguments, which will be filled with user-specified data as above.
      \code{x} may be assumed to be a named vector of length \code{nvars}, (which has the same meanings as above).
      \code{t} is a scalar quantity, the time at which the measurement is made.
      \code{params} may be assumed to be a named vector of length \code{npars}.
      
      \code{rmeasure} must return a named vector.
      If \code{y} is the returned vector, then \code{length(y)=nobs}, where \code{nobs} is the number of observable variables.
    }
    \item{\code{dmeasure}}{
      if provided, must take at least the arguments \code{y}, \code{x}, \code{times}, \code{params}, \code{log}, and \code{\dots}.
      \code{y} may be assumed to be a named vector of length \code{nobs} containing (actual or simulated) values of the observed variables;
      \code{x} will be a named vector of length \code{nvar} containing state variables
      \code{params}, a named vector containing parameters;
      and \code{t}, a scalar, the corresponding observation time.
      It may take additional arguments which will be filled with user-specified data as above.
      \code{dmeasure} must return a single numeric value, the p.d.f. of \code{y} given \code{x} at time \code{t}.
      If \code{log=TRUE}, then the log of the p.d.f. is returned.
    }
    \item{\code{skeleton}}{
      If \code{skeleton} is an R function, it must have at least the arguments \code{x}, \code{t}, \code{params}, and \code{\dots}.
      \code{x} is a numeric vector containing the coordinates of a point in state space at which evaluation of the skeleton is desired.
      \code{t} is a numeric value giving the time at which evaluation of the skeleton is desired.
      Of course, these will be irrelevant in the case of an autonomous skeleton.
      \code{params} is a numeric vector holding the parameters.
      The optional argument \code{covars} is a numeric vector containing the values of the covariates at the time \code{t}.
      \code{covars} will have one value for each column of the covariate table specified when the \code{pomp} object was created.
      \code{covars} is constructed from the covariate table (see \code{covar}, below) by interpolation.
      \code{skeleton} may take additional arguments, which will be filled, as above, with user-specified data.
      \code{skeleton} must return a numeric vector of the same length as \code{x}.
      The return value is interpreted as the vectorfield (if the dynamical system is continuous) or the value of the map (if the dynamical system is discrete), at the point \code{x} at time \code{t}.

      If \code{skeleton} is the name of a native routine, this routine must be of prototype "pomp\_vectorfield\_map" as defined in the header "pomp.h" (see the "include" directory in the installed package directory).
    }
    \item{\code{initializer}}{
      if provided, must have at least the arguments \code{params}, \code{t0}, and \code{\dots}.
      \code{params} is a named vector of parameters.
      \code{t0} will be the time at which initial conditions are desired.
      \code{initializer} must return a named vector of initial conditions.
    }
  }
}
\value{An object of class \code{pomp}.}
\examples{
## Simple example: a 2-D Brownian motion, observed with normal error
## first, set up the pomp object:

rw2 <- pomp(
            rprocess = function (xstart, times, params, ...) { 
              ## this function simulates two independent random walks with intensities s1, s2
              nreps <- ncol(params)
              ntimes <- length(times)
              dt <- diff(times)
              x <- array(0,dim=c(2,nreps,ntimes))
              rownames(x) <- rownames(xstart)
              noise.sds <- params[c('s1','s2'),]
              x[,,1] <- xstart
              for (j in 2:ntimes) {
                ## we are mimicking a continuous-time process, so the increments have SD ~ sqrt(dt)
                ## note that we do not have to assume that 'times' are equally spaced
                x[,,j] <- rnorm(n=2*nreps,mean=x[,,j-1],sd=noise.sds*dt[j-1])
              }
              x
            },
            dprocess = function (x, times, params, log, ...) { 
              ## given a sequence of consecutive states in 'x', this function computes the p.d.f.
              nreps <- ncol(params)
              ntimes <- length(times)
              dt <- diff(times)
              d <- array(0,dim=c(2,nreps,ntimes-1))
              noise.sds <- params[c('s1','s2'),]
              for (j in 2:ntimes)
                d[,,j-1] <- dnorm(x[,,j]-x[,,j-1],mean=0,sd=noise.sds*dt[j-1],log=TRUE)
              if (log) {
                apply(d,c(2,3),sum)
              } else {
                exp(apply(d,c(2,3),sum))
              }
            },
            measurement.model=list(
              y1 ~ norm(mean=x1,sd=tau),
              y2 ~ norm(mean=x2,sd=tau)
            ),
            skeleton = function (x, t, params, covars, ...) {
              ## since this is a continuous-time process, the skeleton is the vectorfield
              ## since the random walk is unbiased, the vectorfield is identically zero
              rep(0,length=length(x))
            },
            times=1:100,
            data=rbind(
              y1=rep(0,100),
              y2=rep(0,100)
              ),
            t0=0
            )

## specify some parameters
p <- rbind(s1=c(2,2,3),s2=c(0.1,1,2),tau=c(1,5,0),x1.0=c(0,0,5),x2.0=c(0,0,0))

## simulate
examples <- simulate(rw2,params=p)
rw2 <- examples[[1]]  ## by default, simulate produces a list of pomp objects
plot(rw2)

t <- seq(0,20)
X <- simulate(rw2,params=p[,1],nsim=10,states=TRUE,obs=TRUE,times=t)

## compute the process model likelihoods
f <- dprocess(
              rw2,
              x=X$states,
              times=t,
              params=matrix(
                p[,1],
                nrow=nrow(p),
                ncol=10,
                dimnames=list(rownames(p),NULL)
                ),
              log=TRUE
              )
apply(f,1,sum)

## compute the measurement likelihoods
g <- dmeasure(
              rw2,
              y=X$obs[,7,],
              x=X$states,
              times=t,
              params=matrix(
                p[,1],
                nrow=nrow(p),
                ncol=10,
                dimnames=list(rownames(p),NULL)
                ),
              log=TRUE
              )
apply(g,1,sum)

## For more examples, see the vignettes.
}
\section{Warning}{
  Some error checking is done, but complete error checking is impossible.  
  If the user-specified functions do not conform to the above specifications (see Details), then the results may be invalid.
  Each algorithm that uses a \code{pomp} object uses some subset of the four basic components (\code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}).
  It is unlikely that any algorithm will use all of these.
}
\references{}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{
  \link{pomp-class},
  \link{pomp-methods},
  \code{\link{rprocess}},
  \code{\link{dprocess}},
  \code{\link{rmeasure}},
  \code{\link{dmeasure}},
  \code{\link{skeleton}},
  \code{\link{init.state}},
  \link{euler}
}
\keyword{models}
\keyword{ts}
