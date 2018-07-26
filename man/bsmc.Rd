\name{Bayesian sequential Monte Carlo}
\title{The Liu and West Bayesian particle filter}
\alias{Bayesian sequential Monte Carlo}
\alias{bsmc2}
\alias{bsmc2-pomp}
\alias{bsmc2,ANY-method}
\alias{bsmc2,missing-method}
\alias{bsmc2,pomp-method}
\alias{prior_samples}
\alias{prior_samples,bsmcd.pomp-method}
\alias{prior_samples-bsmcd.pomp}
\alias{posterior_samples}
\alias{posterior_samples,bsmcd.pomp-method}
\alias{posterior_samples-bsmcd.pomp}
\alias{plot,bsmcd.pomp-method}
\alias{plot-bsmcd.pomp}
\alias{logEvidence}
\alias{logEvidence,bsmcd.pomp-method}
\alias{logEvidence-bsmcd.pomp}
\alias{cond.logEvidence}
\alias{cond.logEvidence,bsmcd.pomp-method}
\alias{cond.logEvidence-bsmcd.pomp}
\alias{eff.sample.size,bsmcd.pomp-method}
\alias{eff.sample.size-bsmcd.pomp}
\description{
  Modified version of the Liu and West (2001) algorithm.
}
\usage{
\S4method{bsmc2}{pomp}(object, params, Np, est, smooth = 0.1,
     tol = 1e-17, max.fail = 0, transform = FALSE, \dots,
     verbose = getOption("verbose"))
}
\arguments{
  \item{object}{
    An object of class \code{pomp} or inheriting class \code{pomp}.
  }
  \item{params, Np}{
    Specifications for the prior distribution of particles.
    See details below.
  }
  \item{est}{
    Names of the parameters that are to be estimated.
    No updates will be made to the other parameters.
    If \code{est} is not specified, all parameters for which there is variation in \code{params} will be estimated.
  }
  \item{smooth}{
    Kernel density smoothing parameter.
    The compensating shrinkage factor will be \code{sqrt(1-smooth^2)}.
    Thus, \code{smooth=0} means that no noise will be added to parameters.
    The general recommendation is that the value of \code{smooth} should be chosen close to 0 (e.g., \code{shrink} ~ 0.1).
  }
  \item{tol, max.fail}{
    Particles with log likelihood below \code{tol} are considered to be \dQuote{lost}.
    A filtering failure occurs when, at some time point, all particles are lost.
    When all particles are lost, the conditional log likelihood at that time point is set to be \code{log(tol)}.
    \code{max.fail} gives the maximum number of filtering failures tolerated.
    If the number of filtering failures exceeds this number, execution will terminate with an error.
  }
  \item{transform}{
    logical;
    if \code{TRUE}, the algorithm operates on the transformed scale.
  }
  \item{\dots}{
    passed to \code{pomp} to supply or modify model components or characteristics.
  }
  \item{verbose}{
    logical; if \code{TRUE}, print diagnostic messages.
  }
}
\value{
  An object of class \dQuote{bsmcd.pomp}.
  One can plot the results using the \code{plot} method.
  \code{plot(x, \dots, pars, thin)}, when \code{x} is a \code{bsmcd.pomp} object, produces density and scatter plots of the prior and posterior distributions.
  One can restrict to particular parameters by naming them in the \code{pars} argument.
  By default, all the samples are plotted.
  When the number of samples is very large, it can be helpful to plot a subsample:
  \code{thin} specifies the size of this random subsample.

  One can also extract various items from such an object:
  \describe{
    \item{\code{prior_samples()}}{
      extracts the samples from the prior distribution that were used in the computation.
    }
    \item{\code{posterior_samples()}}{
      extracts the samples from the (approximate) posterior distribution.
    }
    \item{\code{eff.sample.size()}}{
      returns the effective sample size at each observation time.
    }
    \item{\code{logEvidence()} and \code{cond.logEvidence()}}{
      returns the log evidence of the data and the log evidence of each observation conditional on the preceding ones, respectively.
    }
    \item{\code{coef()}}{
      returns a vector containing the posterior medians.
    }
 }
}
\details{
  There are two ways to specify the prior distribution of particles.
  If \code{params} is unspecified or is a named vector, \code{Np} draws are made from the prior distribution, as specified by \code{\link{rprior}}.
  Alternatively, \code{params} can be specified as an \code{npars} x \code{Np} matrix (with rownames).

  \code{bsmc2} uses a version of the original algorithm, but discards the auxiliary particle filter.
  The modification appears to give superior performance for the same amount of effort.
}
\author{
  Michael Lavine (lavine at math dot umass dot edu),
  Matthew Ferrari (mferrari at psu dot edu),
  Aaron A. King (kingaa at umich dot edu),
  Edward L. Ionides (ionides at umich dot edu)
}
\references{
  Liu, J. and M. West.
  Combining Parameter and State Estimation in Simulation-Based Filtering.
  In A. Doucet, N. de Freitas, and N. J. Gordon, editors,
  Sequential Monte Carlo Methods in Practice, pages 197-224.
  Springer, New York, 2001.
}
\seealso{
  \code{\link{pomp}}, \code{\link{pfilter}}
}
\keyword{ts}
