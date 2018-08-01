\name{Low-level-interface}
\docType{methods}
\title{pomp low-level interface}
\alias{pomp low-level interface}
\alias{rprocess}
\alias{rprocess,pomp-method}
\alias{rprocess-pomp}
\alias{dprocess}
\alias{dprocess,pomp-method}
\alias{dprocess-pomp}
\alias{rmeasure}
\alias{rmeasure,pomp-method}
\alias{rmeasure-pomp}
\alias{dmeasure}
\alias{dmeasure,pomp-method}
\alias{dmeasure-pomp}
\alias{dprior}
\alias{dprior,pomp-method}
\alias{dprior-pomp}
\alias{rprior}
\alias{rprior,pomp-method}
\alias{rprior-pomp}
\alias{init.state}
\alias{init.state,pomp-method}
\alias{init.state-pomp}
\alias{skeleton}
\alias{skeleton,pomp-method}
\alias{skeleton-pomp}
\alias{trajectory}
\alias{trajectory,pomp-method}
\alias{trajectory-pomp}
\alias{pompLoad}
\alias{pompLoad-pomp}
\alias{pompLoad,pomp-method}
\alias{pompUnload}
\alias{pompUnload-pomp}
\alias{pompUnload,pomp-method}
\alias{show,unshowable-method}
\alias{show,listies-method}
\alias{print,unshowable-method}
\alias{print,listies-method}
\description{
  A \sQuote{pomp} object implements a partially observed Markov process (\acronym{POMP}) model.
  Basic operations on this model (with shorthand terms) include:
  \enumerate{
    \item simulation of the state process given parameters (\code{rprocess})
    \item evaluation of the likelihood of a given state trajectory given parameters (\code{dprocess})
    \item simulation of the observation process given the states and parameters (\code{rmeasure})
    \item evaluation of the likelihood of a set of observations given the states and parameters (\code{dmeasure})
    \item simulation from the prior probability distribution (\code{rprior})
    \item evaluation of the prior probability density (\code{dprior})
    \item simulation from the distribution of initial states, given parameters (\code{init.state})
    \item evaluation of the deterministic skeleton at a point in state space, given parameters (\code{skeleton})
    \item computation of a trajectory of the deterministic skeleton given parameters (\code{trajectory})
  }
  \pkg{pomp} provides S4 methods that implement each of these basic operations.
  These operations can be combined to implement statistical inference methods that depend only on a model's \acronym{POMP} structure.
  For convenience, parameter transformations may also be enclosed in a \sQuote{pomp} object.

  This page documents these elements.
}
\usage{
\S4method{rprocess}{pomp}(object, xstart, times, params, offset = 0, \dots)
\S4method{dprocess}{pomp}(object, x, times, params, log = FALSE, \dots)
\S4method{rmeasure}{pomp}(object, x, times, params, \dots)
\S4method{dmeasure}{pomp}(object, y, x, times, params, log = FALSE, \dots)
\S4method{dprior}{pomp}(object, params, log = FALSE, \dots)
\S4method{rprior}{pomp}(object, params, \dots)
\S4method{init.state}{pomp}(object, params, t0, nsim, \dots)
\S4method{skeleton}{pomp}(object, x, t, params, \dots)
\S4method{trajectory}{pomp}(object, params, times, t0, as.data.frame = FALSE, \dots,
          verbose = getOption("verbose", FALSE))
\S4method{pompLoad}{pomp}(object, \dots)
\S4method{pompUnload}{pomp}(object, \dots)
}
\arguments{
  \item{object}{an object of class \sQuote{pomp}.}
  \item{xstart}{
    an \code{nvar} x \code{nrep} matrix containing the starting state of the system.
    Columns of \code{xstart} correspond to states; rows to components of the state vector.
    One independent simulation will be performed for each column.
    Note that in this case, \code{params} must also have \code{nrep} columns.
  }
  \item{x}{
    a rank-3 array containing states of the unobserved process.
    The dimensions of \code{x} are \code{nvars} x \code{nrep} x \code{ntimes}, where \code{nvars} is the number of state variables, \code{nrep} is the number of replicates, and \code{ntimes} is the length of \code{times}.
  }
  \item{y}{
    a matrix containing observations.
    The dimensions of \code{y} are \code{nobs} x \code{ntimes}, where \code{nobs} is the number of observables and \code{ntimes} is the length of \code{times}.
  }
  \item{times, t}{
    a numeric vector (length \code{ntimes}) containing times.
    These must be in non-decreasing order.
  }
  \item{params}{
    a \code{npar} x \code{nrep} matrix of parameters.
    Each column is an independent parameter set and is paired with the corresponding column of \code{x} or \code{xstart}.

    In the case of \code{init.state}, \code{params} is a named vector of parameters.
  }
  \item{offset}{
    integer;
    the first \code{offset} times in \code{times} will not be returned.
  }
  \item{t0}{the initial time at which initial states are requested.}
  \item{nsim}{
    optional integer; the number of initial states to simulate.
    By default, this is equal to the number of columns of \code{params}.
  }
  \item{log}{if TRUE, log probabilities are returned.}
  \item{as.data.frame}{
    logical; if \code{TRUE}, return the result as a data-frame.
  }
  \item{\dots}{
    In \code{trajectory}, additional arguments are passed to the ODE integrator (if the skeleton is a vectorfield) and ignored if it is a map.
    See \code{\link[deSolve]{ode}} for a description of the additional arguments accepted.

    In all other cases, additional arguments are ignored.
  }
  \item{verbose}{
    logical; if \code{TRUE}, more information will be displayed.
  }
}
\section{rprocess}{
  \code{rprocess} simulates the process-model portion of partially-observed Markov process.

  When \code{rprocess} is called, the first entry of \code{times} is taken to be the initial time
  (i.e., that corresponding to \code{xstart}).
  Subsequent times are the additional times at which the state of the simulated processes are required.

  \code{rprocess} returns a rank-3 array with rownames.
  Suppose \code{x} is the array returned.
  Then \preformatted{dim(x)=c(nvars,nrep,ntimes-offset),} where \code{nvars} is the number of state variables (=\code{nrow(xstart)}), \code{nrep} is the number of independent realizations simulated (=\code{ncol(xstart)}), and \code{ntimes} is the length of the vector \code{times}.
  \code{x[,j,k]} is the value of the state process in the \code{j}-th realization at time \code{times[k+offset]}.
  The rownames of \code{x} must correspond to those of \code{xstart}.
}
\section{dprocess}{
  \code{dprocess} evaluates the probability density of a sequence of consecutive state transitions.

  \code{dprocess} returns a matrix of dimensions \code{nrep} x \code{ntimes-1}.
  If \code{d} is the returned matrix, \code{d[j,k]} is the likelihood of the transition from state \code{x[,j,k-1]} at time \code{times[k-1]} to state \code{x[,j,k]} at time \code{times[k]}.
}
\section{rmeasure}{
  \code{rmeasure} simulate the measurement model given states and parameters.

  \code{rmeasure} returns a rank-3 array of dimensions \code{nobs} x \code{nrep} x \code{ntimes}, where \code{nobs} is the number of observed variables.
}
\section{dmeasure}{
  \code{dmeasure} evaluates the probability density of observations given states.

  \code{dmeasure} returns a matrix of dimensions \code{nreps} x \code{ntimes}.
  If \code{d} is the returned matrix, \code{d[j,k]} is the likelihood of the observation \code{y[,k]} at time \code{times[k]} given the state \code{x[,j,k]}.
}
\section{dprior, rprior}{
  \code{dprior} evaluates the prior probability density and \code{rprior} simulates from the prior.
}
\section{init.state}{
  \code{init.state} returns an \code{nvar} x \code{nsim} matrix of state-process initial conditions when given an \code{npar} x \code{nsim} matrix of parameters, \code{params}, and an initial time \code{t0}.
  By default, \code{t0} is the initial time defined when the \sQuote{pomp} object ws constructed.
  If \code{nsim} is not specified, then \code{nsim=ncol(params)}.
}
\section{skeleton}{
  The method \code{skeleton} evaluates the deterministic skeleton at a point or points in state space, given parameters.
  In the case of a discrete-time system, the skeleton is a map.
  In the case of a continuous-time system, the skeleton is a vectorfield.
  NB: \code{skeleton} just evaluates the deterministic skeleton;
  it does not iterate or integrate.

  \code{skeleton} returns an array of dimensions \code{nvar} x \code{nrep} x \code{ntimes}.
  If \code{f} is the returned matrix, \code{f[i,j,k]} is the i-th component of the deterministic skeleton at time \code{times[k]} given the state \code{x[,j,k]} and parameters \code{params[,j]}.
}
\section{trajectory}{
  \code{trajectory} computes a trajectory of the deterministic skeleton of a Markov process.
  In the case of a discrete-time system, the deterministic skeleton is a map and a trajectory is obtained by iterating the map.
  In the case of a continuous-time system, the deterministic skeleton is a vector-field; \code{trajectory} uses the numerical solvers in \pkg{\link[deSolve]{deSolve}} to integrate the vectorfield.

  \code{trajectory} returns an array of dimensions \code{nvar} x \code{nrep} x \code{ntimes}.
  If \code{x} is the returned matrix, \code{x[i,j,k]} is the i-th component of the state vector at time \code{times[k]} given parameters \code{params[,j]}.

  When the skeleton is a vectorfield, \code{trajectory} integrates it using \code{\link[deSolve]{ode}}.
  When the skeleton is a map, \code{trajectory} iterates it.
  By default, time is advanced 1 unit per iteration.
  The user can change this behavior by specifying the desired timestep using the argument \code{skelmap.delta.t} in the construction of the \sQuote{pomp} object.
}
\section{Parameter transformations}{
  User-defined parameter transformations enclosed in the \sQuote{pomp} object can be accessed via \code{\link{partrans}}.
}
\section{\code{pompLoad}, \code{pompUnload}}{
  \code{pompLoad} and \code{pompUnload} cause compiled codes associated with \code{object} to be dynamically linked or unlinked, respectively.
  When C snippets are used in the construction of a \sQuote{pomp} object, the resulting shared-object library is dynamically loaded (linked) before each use, and unloaded afterward.
  These functions are provided because in some instances, greater control may be desired.
  These functions have no effect on shared-object libraries linked by the user.
}
\examples{
pompExample(ricker)

p <- parmat(c(r=42,c=1,phi=10,sigma=0.3,N_0=7,e.0=0),10)
t <- c(1:10,20,30)
t0 <- 0
x0 <- init.state(ricker,params=p,t0=t0)
x <- rprocess(ricker,xstart=x0,times=c(t0,t),params=p,offset=1)
y <- rmeasure(ricker,params=p,x=x,times=t)
ll <- dmeasure(ricker,y=y[,3,,drop=FALSE],x=x,times=t,params=p,log=TRUE)
apply(ll,1,sum)
f <- skeleton(ricker,x=x,t=t,params=p)
z <- trajectory(ricker,params=p,times=t,t0=t0)

## short arguments are recycled:
p <- c(r=42,phi=10,c=1,sigma=0.3,N_0=7,e.0=0)
t <- c(1:10,20,30)
t0 <- 0
x0 <- init.state(ricker,params=p,t0=t0)
x <- rprocess(ricker,xstart=x0,times=c(t0,t),params=p,offset=1)
y <- rmeasure(ricker,params=p,x=x,times=t)
ll <- dmeasure(ricker,y=y,x=x,times=t,params=p,log=TRUE)
f <- skeleton(ricker,x=x,t=t,params=p)
z <- trajectory(ricker,params=p,times=t,t0=t0)
}
\author{Aaron A. King}
\seealso{\code{\link{pomp}}, \link{pomp methods}}
\keyword{programming}
