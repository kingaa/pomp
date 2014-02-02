\name{mif}
\docType{methods}
\alias{mif}
\alias{mif,mif-method}
\alias{mif-mif}
\alias{mif,pfilterd.pomp-method}
\alias{mif-pfilterd.pomp}
\alias{mif,pomp-method}
\alias{mif-pomp}
\alias{continue}
\alias{continue,mif-method}
\alias{continue-mif}
\title{The MIF algorithm}
\description{
  The MIF algorithm for estimating the parameters of a partially-observed Markov process.
}
\usage{
mif(object, \dots)
\S4method{mif}{pomp}(object, Nmif = 1, start, pars, ivps = character(0),
    particles, rw.sd, Np, ic.lag, var.factor,
    cooling.type, cooling.fraction, cooling.factor,
    method = c("mif","unweighted","fp","mif2"),
    tol = 1e-17, max.fail = Inf,
    verbose = getOption("verbose"), transform = FALSE, \dots)
\S4method{mif}{pfilterd.pomp}(object, Nmif = 1, Np, tol, \dots)
\S4method{mif}{mif}(object, Nmif, start, pars, ivps,
    particles, rw.sd, Np, ic.lag, var.factor,
    cooling.type, cooling.fraction,
    method, tol, transform, \dots)
\S4method{continue}{mif}(object, Nmif = 1, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp}.
  }
  \item{Nmif}{
    The number of MIF iterations to perform.
  }
  \item{start}{
    named numerical vector;
    the starting guess of the parameters.
  }
  \item{pars}{
    optional character vector naming the ordinary parameters to be estimated.
    Every parameter named in \code{pars} must have a positive random-walk standard deviation specified in \code{rw.sd}.
    Leaving \code{pars} unspecified is equivalent to setting it equal to the names of all parameters with a positive value of \code{rw.sd} that are not \code{ivps}.
  }
  \item{ivps}{
    optional character vector naming the initial-value parameters (IVPs) to be estimated.
    Every parameter named in \code{ivps} must have a positive random-walk standard deviation specified in \code{rw.sd}.
    If \code{pars} is empty, i.e., only IVPs are to be estimated, see below \dQuote{"Using MIF to estimate initial-value parameters only"}.
  }
  \item{particles}{
    Function of prototype \code{particles(Np,center,sd,...)} which sets up the starting particle matrix by drawing a sample of size \code{Np} from the starting particle distribution centered at \code{center} and of width \code{sd}.
    If \code{particles} is not supplied by the user, the default behavior is to draw the particles from a multivariate normal distribution with mean \code{center} and standard deviation \code{sd}.
  }
  \item{rw.sd}{
    numeric vector with names; the intensity of the random walk to be applied to parameters.
    The random walk is only applied to parameters named in \code{pars} (i.e., not to those named in \code{ivps}).
    The algorithm requires that the random walk be nontrivial, so each element in \code{rw.sd[pars]} must be positive.
    \code{rw.sd} is also used to scale the initial-value parameters (via the \code{particles} function).
    Therefore, each element of \code{rw.sd[ivps]} must be positive.
    The following must be satisfied:
    \code{names(rw.sd)} must be a subset of \code{names(start)},
    \code{rw.sd} must be non-negative (zeros are simply ignored),
    the name of every positive element of \code{rw.sd} must be in either \code{pars} or \code{ivps}.
  }
  \item{Np}{
    the number of particles to use in filtering.
    This may be specified as a single positive integer, in which case the same number of particles will be used at each timestep.
    Alternatively, if one wishes the number of particles to vary across timestep, one may specify \code{Np} either as a vector of positive integers (of length \code{length(time(object,t0=TRUE))}) or as a function taking a positive integer argument.
    In the latter case, \code{Np(k)} must be a single positive integer, representing the number of particles to be used at the \code{k}-th timestep:
    \code{Np(0)} is the number of particles to use going from \code{timezero(object)} to \code{time(object)[1]},
    \code{Np(1)}, from \code{timezero(object)} to \code{time(object)[1]},
    and so on, while when \code{T=length(time(object,t0=TRUE))},
    \code{Np(T)} is the number of particles to sample at the end of the time-series.
  }
  \item{ic.lag}{
    a positive integer;
    the timepoint for fixed-lag smoothing of initial-value parameters.
    The \code{mif} update for initial-value parameters consists of replacing them by their filtering mean at time \code{times[ic.lag]}, where \code{times=time(object)}.
    It makes no sense to set \code{ic.lag>length(times)};
    if it is so set, \code{ic.lag} is set to \code{length(times)} with a warning.
    For \code{method='mif2'}, the default is \code{ic.lag=length(times)}.
  }
  \item{var.factor}{
    a positive number;
    the scaling coefficient relating the width of the starting particle distribution to \code{rw.sd}.
    In particular, the width of the distribution of particles at the start of the first MIF iteration will be \code{random.walk.sd*var.factor}.
  }
  \item{cooling.type, cooling.fraction, cooling.factor}{
    specifications for the cooling schedule, i.e., the manner in which the intensity of the parameter perturbations is reduced with successive filtering iterations.
    \code{cooling.type} specifies the nature of the cooling schedule.
    When \code{cooling.type="geometric"}, on the n-th MIF iteration, the relative perturbation intensity is \code{cooling.fraction^(n/50)}.
    When \code{cooling.type="hyperbolic"}, on the n-th MIF iteration, the relative perturbation intensity is \code{(s+1)(s+n)}, where \code{(s+1)/(s+50)=cooling.fraction}.
    \code{cooling.fraction} is the relative magnitude of the parameter perturbations after 50 MIF iterations.
    \code{cooling.factor} is now deprecated:
    to achieve the old behavior, use \code{cooling.type="geometric"} and \code{cooling.fraction=(cooling.factor)^50}.
  }
  \item{method}{
    \code{method} sets the update rule used in the algorithm.
    \code{method="mif"} uses the iterated filtering update rule (Ionides 2006, 2011);
    \code{method="unweighted"} updates the parameter to the unweighted average of the filtering means of the parameters at each time;
    \code{method="fp"} updates the parameter to the filtering mean at the end of the time series.
    \code{method="mif2"} [EXPLAIN]
  }
  \item{tol}{
    See the description under \code{\link{pfilter}}.
  }
  \item{max.fail}{
    See the description under \code{\link{pfilter}}.
  }
  \item{verbose}{
    logical; if TRUE, print progress reports.
  }
  \item{transform}{
    logical;
    if \code{TRUE}, optimization is performed on the transformed scale.
  }
  \item{\dots}{
    additional arguments that override the defaults.
  }
}
\section{Re-running MIF Iterations}{
  To re-run a sequence of MIF iterations, one can use the \code{mif} method on a \code{mif} object.
  By default, the same parameters used for the original MIF run are re-used (except for \code{weighted}, \code{tol}, \code{max.fail}, and \code{verbose}, the defaults of which are shown above).
  If one does specify additional arguments, these will override the defaults.
}
\section{Continuing MIF Iterations}{
  One can resume a series of MIF iterations from where one left off using the \code{continue} method.
  A call to \code{mif} to perform \code{Nmif=m} iterations followed by a call to \code{continue} to perform \code{Nmif=n} iterations will produce precisely the same effect as a single call to \code{mif} to perform \code{Nmif=m+n} iterations.
  By default, all the algorithmic parameters are the same as used in the original call to \code{mif}.
  Additional arguments will override the defaults.
}
\section{Using MIF to estimate initial-value parameters only}{
  One can use MIF's fixed-lag smoothing to estimate only initial value parameters (IVPs).
  In this case, \code{pars} is left empty and the IVPs to be estimated are named in \code{ivps}.
  If \code{theta} is the current parameter vector, then at each MIF iteration, \code{Np} particles are drawn from a distribution centered at \code{theta} and with width proportional to \code{var.factor*rw.sd}, a particle filtering operation is performed, and \code{theta} is replaced by the filtering mean at \code{time(object)[ic.lag]}.
  Note the implication that, when \code{mif} is used in this way on a time series any longer than \code{ic.lag}, unnecessary work is done.
  If the time series in \code{object} is longer than \code{ic.lag}, consider replacing \code{object} with \code{window(object,end=ic.lag)}. 
}
\section{Details}{
  If \code{particles} is not specified, the default behavior is to draw the particles from a multivariate normal distribution.
  \strong{It is the user's responsibility to ensure that, if the optional \code{particles} argument is given, that the \code{particles} function satisfies the following conditions:}

  \code{particles} has at least the following arguments:
  \code{Np}, \code{center}, \code{sd}, and \code{\dots}.
  \code{Np} may be assumed to be a positive integer;
  \code{center} and \code{sd} will be named vectors of the same length.
  Additional arguments may be specified;
  these will be filled with the elements of the \code{userdata} slot of the underlying \code{pomp} object (see \code{\link{pomp-class}}).

  \code{particles} returns a \code{length(center)} x \code{Np} matrix with rownames matching the names of \code{center} and \code{sd}.
  Each column represents a distinct particle.

  The center of the particle distribution returned by \code{particles} should be \code{center}.
  The width of the particle distribution should vary monotonically with \code{sd}.
  In particular, when \code{sd=0}, the \code{particles} should return matrices with \code{Np} identical columns, each given by the parameters specified in \code{center}.
}
\references{
  E. L. Ionides, C. Bret\\'o, & A. A. King,
  Inference for nonlinear dynamical systems,
  Proc. Natl. Acad. Sci. U.S.A., 103:18438--18443, 2006.

  E. L. Ionides, A. Bhadra, Y. Atchad{\\'e}, & A. A. King,
  Iterated filtering,
  Annals of Statistics, 39:1776--1802, 2011.

  A. A. King, E. L. Ionides, M. Pascual, and M. J. Bouma,
  Inapparent infections and cholera dynamics,
  Nature, 454:877--880, 2008.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{
  \code{\link{mif-methods}}, \code{\link{pomp}}, \code{\link{pomp-class}}, \code{\link{pfilter}}.
  See the \dQuote{intro_to_pomp} vignette for examples.
}
\keyword{ts}
