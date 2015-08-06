\name{Iterated filtering}
\title{Maximum likelihood by iterated filtering}
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
\alias{mif-class}
\alias{mif-methods}
\alias{conv.rec}
\alias{conv.rec,mif-method}
\alias{conv.rec-mif}
\alias{conv.rec,mifList-method}
\alias{conv.rec-mifList}
\alias{coef,mifList-method}
\alias{coef.rec-mifList}
\alias{plot-mif}
\alias{plot,mif-method}
\alias{plot-mifList}
\alias{plot,mifList-method}
\alias{mifList-class}
\alias{c-mif}
\alias{c,mif-method}
\alias{c-mifList}
\alias{c,mifList-method}
\alias{[-mifList}
\alias{[,mifList-method}
\description{
  Iterated filtering algorithms for estimating the parameters of a partially-observed Markov process.
  Running \code{mif} causes the iterated filtering algorithm to run for a specified number of iterations.
  At each iteration, the particle filter is performed on a perturbed version of the model.
  Specifically, parameters to be estimated are subjected to random perturbations at each observation.
  This extra variability effectively smooths the likelihood surface and combats particle depletion by introducing diversity into the population of particles.
  At the iterations progress, the magnitude of the perturbations is diminished according to a user-specified cooling schedule.
  For most purposes, \code{mif} has been superseded by \code{\link{mif2}}.
}
\usage{
\S4method{mif}{pomp}(object, Nmif = 1, start, ivps = character(0),
    particles, rw.sd, Np, ic.lag, var.factor = 1,
    cooling.type, cooling.fraction.50,
    method = c("mif","unweighted","fp","mif2"),
    tol = 1e-17, max.fail = Inf,
    verbose = getOption("verbose"), transform = FALSE, \dots)
\S4method{mif}{pfilterd.pomp}(object, Nmif = 1, Np, tol, \dots)
\S4method{mif}{mif}(object, Nmif, start, ivps,
    particles, rw.sd, Np, ic.lag, var.factor,
    cooling.type, cooling.fraction.50,
    method, tol, transform, \dots)
\S4method{continue}{mif}(object, Nmif = 1, \dots)
\S4method{conv.rec}{mif}(object, pars, transform = FALSE, \dots)
\S4method{conv.rec}{mifList}(object, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp}.
  }
  \item{Nmif}{
    The number of filtering iterations to perform.
  }
  \item{start}{
    named numerical vector;
    the starting guess of the parameters.
  }
  \item{ivps}{
    optional character vector naming the initial-value parameters (IVPs) to be estimated.
    Every parameter named in \code{ivps} must have a positive random-walk standard deviation specified in \code{rw.sd}.
    If there are no regular parameters with positive \code{rw.sd}, i.e., only IVPs are to be estimated, see below \dQuote{"Using \code{mif} to estimate initial-value parameters only"}.
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

    For \code{method="mif2"}, the default is \code{ic.lag=length(times)}.
  }
  \item{var.factor}{
    optional positive scalar;
    the scaling coefficient relating the width of the starting particle distribution to \code{rw.sd}.
    In particular, the width of the distribution of particles at the start of the first \code{mif} iteration will be \code{random.walk.sd*var.factor}.
    By default, \code{var.factor=1}.
  }
  \item{cooling.type, cooling.fraction.50}{
    specifications for the cooling schedule, i.e., the manner in which the intensity of the parameter perturbations is reduced with successive filtering iterations.
    \code{cooling.type} specifies the nature of the cooling schedule.

    When \code{cooling.type="geometric"}, on the n-th \code{mif} iteration, the relative perturbation intensity is \code{cooling.fraction.50^(n/50)}.

    When \code{cooling.type="hyperbolic"}, on the n-th \code{mif} iteration, the relative perturbation intensity is \code{(s+1)/(s+n)}, where \code{(s+1)/(s+50)=cooling.fraction.50}.
    \code{cooling.fraction.50} is the relative magnitude of the parameter perturbations after 50 \code{mif} iterations.
  }
  \item{method}{
    \code{method} sets the update rule used in the algorithm.
    \code{method="mif"} uses the iterated filtering update rule (Ionides 2006, 2011);
    \code{method="unweighted"} updates the parameter to the unweighted average of the filtering means of the parameters at each time;
    \code{method="fp"} updates the parameter to the filtering mean at the end of the time series.
    \code{method="mif2"} implements an incomplete version of the iterated Bayes map method of Ionides (2015).
    The latter method is, by every indication, both more efficient and more stable.
  }
  \item{tol, max.fail}{
    See the description under \code{\link{pfilter}}.
  }
  \item{verbose}{
    logical; if TRUE, print progress reports.
  }
  \item{transform}{
    logical;
    if \code{TRUE}, optimization is performed on the transformed scale, as defined by the user-supplied parameter transformations (see \code{\link{pomp}}).
  }
  \item{\dots}{
    additional arguments that override the defaults.
  }
  \item{pars}{names of parameters.}
}
\value{
  Upon successful completion, \code{mif} returns an object of class \code{mif}.
  The latter inherits from the \code{\link{pfilterd.pomp}} and \code{\link{pomp}} classes.
}
\section{IF2}{
  A more full-featured version of the improved iterated filtering algorithm (IF2) is implemented as \code{\link{mif2}}.
}
\section{Regular parameters vs initial-value parameters}{
  Initial-value parameters (IVPs) differ from regular parameters in that the majority of the information about these parameters is restricted to the early part of the time series.
  That is, increasing the length of the time series provides progressively less additional information about IVPs than it does about regular parameters.
  In \code{mif}, while regular parameters are perturbed at the initial time and after every observation, IVPs are perturbed only at the initial time.
}
\section{Re-running \code{mif} Iterations}{
  To re-run a sequence of \code{mif} iterations, one can use the \code{mif} method on a \code{mif} object.
  By default, the same parameters used for the original \code{mif} run are re-used (except for \code{tol}, \code{max.fail}, and \code{verbose}, the defaults of which are shown above).
  If one does specify additional arguments, these will override the defaults.
}
\section{Continuing \code{mif} Iterations}{
  One can resume a series of \code{mif} iterations from where one left off using the \code{continue} method.
  A call to \code{mif} to perform \code{Nmif=m} iterations followed by a call to \code{continue} to perform \code{Nmif=n} iterations will produce precisely the same effect as a single call to \code{mif} to perform \code{Nmif=m+n} iterations.
  By default, all the algorithmic parameters are the same as used in the original call to \code{mif}.
  Additional arguments will override the defaults.
}
\section{Using \code{mif} to estimate initial-value parameters only}{
  One can use \code{mif}'s fixed-lag smoothing to estimate only initial value parameters (IVPs).
  In this case, \code{pars} is left empty and the IVPs to be estimated are named in \code{ivps}.
  If \code{theta} is the current parameter vector, then at each \code{mif} iteration, \code{Np} particles are drawn from a distribution centered at \code{theta} and with width proportional to \code{var.factor*rw.sd}, a particle filtering operation is performed, and \code{theta} is replaced by the filtering mean at \code{time(object)[ic.lag]}.
  Note the implication that, when \code{mif} is used in this way on a time series any longer than \code{ic.lag}, unnecessary work is done.
  If the time series in \code{object} is longer than \code{ic.lag}, consider replacing \code{object} with \code{window(object,end=ic.lag)}. 
}
\section{Methods}{
  Methods that can be used to manipulate, display, or extract information from a \code{mif} object:
  \describe{
    \item{conv.rec}{
      \code{conv.rec(object, pars = NULL)} returns the columns of the convergence-record matrix corresponding to the names in \code{pars}.
      By default, all rows are returned.
    }
    \item{logLik}{
      Returns the value in the \code{loglik} slot.
      NB: this is \emph{not} the same as the likelihood of the model at the MLE!
    }
    \item{c}{
      Concatenates \code{mif} objects into a \code{mifList}.
    }
    \item{plot}{
      Plots a series of diagnostic plots when applied to a \code{mif} or \code{mifList} object.
    }
  }
}
\section{Details}{
  If \code{particles} is not specified, the default behavior is to draw the particles from a multivariate normal distribution.
  \strong{It is the user's responsibility to ensure that, if the optional \code{particles} argument is given, that the \code{particles} function satisfies the following conditions:}

  \code{particles} has at least the following arguments:
  \code{Np}, \code{center}, \code{sd}, and \code{\dots}.
  \code{Np} may be assumed to be a positive integer;
  \code{center} and \code{sd} will be named vectors of the same length.
  Additional arguments may be specified;
  these will be filled with the elements of the \code{userdata} slot of the underlying \code{pomp} object (see \code{\link{pomp}}).

  \code{particles} returns a \code{length(center)} x \code{Np} matrix with rownames matching the names of \code{center} and \code{sd}.
  Each column represents a distinct particle.

  The center of the particle distribution returned by \code{particles} should be \code{center}.
  The width of the particle distribution should vary monotonically with \code{sd}.
  In particular, when \code{sd=0}, the \code{particles} should return matrices with \code{Np} identical columns, each given by the parameters specified in \code{center}.
}
\references{
  E. L. Ionides, C. Breto, & A. A. King,
  Inference for nonlinear dynamical systems,
  Proc. Natl. Acad. Sci. U.S.A., 103:18438--18443, 2006.

  E. L. Ionides, A. Bhadra, Y. Atchad{\\'e}, & A. A. King,
  Iterated filtering,
  Annals of Statistics, 39:1776--1802, 2011.

  E. L. Ionides, D. Nguyen, Y. Atchad{\\'e}, S. Stoev, and A. A. King. 
  Inference for dynamic and latent variable models via iterated, perturbed Bayes maps.
  Proc. Natl. Acad. Sci. U.S.A., 112:719--724, 2015.

  A. A. King, E. L. Ionides, M. Pascual, and M. J. Bouma,
  Inapparent infections and cholera dynamics,
  Nature, 454:877--880, 2008.
}
\author{Aaron A. King}
\seealso{
  \code{\link{pomp}}, \code{\link{pfilter}}, \code{\link{mif2}}
}
\keyword{optimize}
\keyword{ts}
