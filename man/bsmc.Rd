\name{Bayesian sequential Monte Carlo}
\title{The Liu and West Bayesian particle filter}
\alias{Bayesian sequential Monte Carlo}
\alias{bsmc}
\alias{bsmc-pomp}
\alias{bsmc,pomp-method}
\alias{bsmc2}
\alias{bsmc2-pomp}
\alias{bsmc2,pomp-method}
\alias{$,bsmcd.pomp-method}
\alias{$-bsmcd.pomp}
\alias{plot,bsmcd.pomp-method}
\alias{plot-bsmcd.pomp}
\description{
  Modified versions of the Liu and West (2001) algorithm.
}
\usage{
\S4method{bsmc2}{pomp}(object, params, Np, est, smooth = 0.1,
     tol = 1e-17, verbose = getOption("verbose"), max.fail = 0,
     transform = FALSE, \dots)
\S4method{bsmc}{pomp}(object, params, Np, est, smooth = 0.1,
     ntries = 1, tol = 1e-17, lower = -Inf, upper = Inf,
     verbose = getOption("verbose"), max.fail = 0,
     transform = FALSE, \dots)
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
    Names of the rows of \code{params} that are to be estimated.
    No updates will be made to the other parameters.
    If \code{est} is not specified, all parameters for which there is variation in \code{params} will be estimated.
  }
  \item{smooth}{
    Kernel density smoothing parameters.
    The compensating shrinkage factor will be \code{sqrt(1-smooth^2)}.
    Thus, \code{smooth=0} means that no noise will be added to parameters.
    Generally, the value of \code{smooth} should be chosen close to 0 (i.e., \code{shrink~0.1}).
  }
  \item{ntries}{
    Number of draws from \code{rprocess} per particle used to estimate the expected value of the state process at time \code{t+1} given the state and parameters at time \code{t}.
  }
  \item{tol}{
    Particles with log likelihood below \code{tol} are considered to be \dQuote{lost}.
    A filtering failure occurs when, at some time point, all particles are lost.
    When all particles are lost, the conditional log likelihood at that time point is set to be \code{log(tol)}.
  }
  \item{lower, upper}{
    optional; lower and upper bounds on the priors.
    This is useful in case there are box constraints satisfied by the priors.
    The posterior is guaranteed to lie within these bounds.
  }
  \item{verbose}{
    logical; if \code{TRUE}, print diagnostic messages.
  }
  \item{max.fail}{
    The maximum number of filtering failures allowed.
    If the number of filtering failures exceeds this number, execution will terminate with an error.
  }
  \item{transform}{
    logical;
    if \code{TRUE}, the algorithm operates on the transformed scale.
  }
  \item{\dots}{
    currently ignored.
  }
}
\value{
  An object of class \dQuote{bsmcd.pomp}.
  The \dQuote{params} slot of this object will hold the parameter posterior medians.
  The slots of this class include:
  \item{post}{
    A matrix containing draws from the approximate posterior distribution.
  }
  \item{prior}{
    A matrix containing draws from the prior distribution (identical to \code{params} on call).
  }
  \item{eff.sample.size}{
    A vector containing the effective number of particles at each time point.
  }
  \item{smooth}{
    The smoothing parameter used (see above).
  }
  \item{nfail}{
    The number of filtering failures encountered.
  }
  \item{cond.log.evidence}{
    A vector containing the conditional log evidence scores at each time point.
  }
  \item{log.evidence}{
    The estimated log evidence.
  }
  \item{weights}{
    The resampling weights for each particle.
  }
}
\details{
  There are two ways to specify the prior distribution of particles.
  If \code{params} is unspecified or is a named vector, \code{Np} draws are made from the prior distribution, as specified by \code{\link{rprior}}.
  Alternatively, \code{params} can be specified as an \code{npars} x \code{Np} matrix (with rownames).

  \code{bsmc} uses version of the original algorithm that includes a plug-and-play auxiliary particle filter.
  \code{bsmc2} discards this auxiliary particle filter and appears to give superior performance for the same amount of effort.
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
