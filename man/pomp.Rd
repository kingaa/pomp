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
    These are simply copied over as initial conditions when \code{init.state} is called (see \code{\link{init.state-pomp}}).
    The names of the state variables are the same as the corresponding initial value parameters, but with the ".0" dropped.
  }
  \item{\dots}{
    Any additional arguments are stored in a slot `userdata' and are passed as arguments to each of the functions \code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, and \code{initializer} whenever they are evaluated.
    Using this mechanism, the user can store additional information necessary for the definition of the model.
  }
}
\details{
  \strong{It is the user's responsibility to ensure that the \code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, and \code{initializer} functions satisfy the following conditions:}
  
  \describe{
    \item{\code{rprocess}}{
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
    }
    \item{\code{dprocess}}{
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
    }
    \item{\code{rmeasure}}{
      \code{rmeasure} must have at least the arguments \code{x}, \code{times}, \code{params}, and \code{\dots}.
      It may take additional arguments, which will be filled with user-specified data as above.
      \code{x} must be a rank-3 array dimension \code{c(nvars,nreps,ntimes)}, where these variables have the same meanings as above.
      \code{times} is the corresponding set of times.
      \code{params} is a rank-2 array of dimension \code{c(npars,nreps)} as above.
      
      \code{rmeasure} must return a rank-3 array with rownames.
      If \code{y} is the returned array, then \code{dim(y)=c(nobs,nreps,ntimes)}, where \code{nobs} is the number of observable variables and \code{nreps}, \code{ntimes} agree with the corresponding dimensions of \code{x}.
      \code{y[,j,k]} must be the vector of observables in the \code{j}-th realization at time \code{times[k]}.
    }
    \item{\code{dmeasure}}{
      \code{dmeasure} must have at least the arguments \code{y}, \code{x}, \code{times}, \code{params}, \code{log}, and \code{\dots}.
      \code{y} is a rank-2 array of observations (\code{nobs} x \code{ntimes});
      \code{x}, a rank-3 array of states;
      \code{params}, a rank-2 array containing parameters;
      and \code{times}, the corresponding observation times.
      It may take additional arguments which will be filled with user-specified data as above.
      \code{dmeasure} must return a rank-2 array of dimension \code{nreps} x \code{ntimes}.
      If \code{d} is the returned array, then \code{d[j,k]} is the p.d.f. of \code{y[,k]} given \code{x[,j,k]} at time \code{times[k]}.
      If \code{log=TRUE}, then the log of the p.d.f. is returned.
    }
    \item{\code{initializer}}{
      \code{initializer} must have at least the arguments \code{params}, \code{t0}, and \code{\dots}.
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
              nsims <- ncol(params)
              ntimes <- length(times)
              dt <- diff(times)
              x <- array(0,dim=c(2,nsims,ntimes))
              rownames(x) <- rownames(xstart)
              noise.sds <- params[c('s1','s2'),]
              x[,,1] <- xstart
              for (j in 2:ntimes) {
                ## we are mimicking a continuous-time process, so the increments have SD ~ sqrt(dt)
                ## note that we do not have to assume that 'times' are equally spaced
                x[,,j] <- rnorm(n=2*nsims,mean=x[,,j-1],sd=noise.sds*dt[j-1])
              }
              x
            },
            dprocess = function (x, times, params, log = FALSE, ...) { 
              ## given a sequence of consecutive states in 'x', this function computes the p.d.f.
              nsims <- ncol(params)
              ntimes <- length(times)
              dt <- diff(times)
              d <- array(0,dim=c(2,nsims,ntimes-1))
              noise.sds <- params[c('s1','s2'),]
              for (j in 2:ntimes)
                d[,,j-1] <- dnorm(x[,,j]-x[,,j-1],mean=0,sd=noise.sds*dt[j-1],log=TRUE)
              if (log) {
                apply(d,c(2,3),sum)
              } else {
                exp(apply(d,c(2,3),sum))
              }
            },
            rmeasure = function (x, times, params, ...) {
              ## noisy observations of the two walks with common noise SD 'tau'
              nsims <- dim(x)[2]
              ntimes <- dim(x)[3]
              y <- array(0,dim=c(2,nsims,ntimes))
              rownames(y) <- c('y1','y2')
              for (j in 1:nsims) {
                for (k in 1:ntimes) {
                  y[,j,k] <- rnorm(2,mean=x[,j,k],sd=params['tau',j])
                }
              }
              y
            },
            dmeasure = function (y, x, times, params, log = FALSE, ...) {
              ## noisy observations of the two walks with common noise SD 'tau'
              d1 <- dnorm(
                          x=y['y1',],
                          mean=x['x1',,],
                          sd=params['tau',],
                          log=TRUE
                          )
              d2 <- dnorm(
                          x=y['y2',],
                          mean=x['x2',,],
                          sd=params['tau',],
                          log=TRUE
                          )
              if (log) {
                d1+d2
              } else {
                exp(d1+d2)
              }
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
}
\references{}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\link{pomp-class}, \link{euler}, \code{vignette(intro_to_pomp)}}
\keyword{models}
\keyword{ts}
