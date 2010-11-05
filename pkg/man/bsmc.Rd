\name{bsmc}
\alias{bsmc}
\alias{bsmc-pomp}
\alias{bsmc,pomp-method}
\title{Liu and West Bayesian Particle Filter}
\description{
  Generates draws from the posterior distribution for the parameters using the Liu and West algorithm.
  \code{bsmc} gives draws from the posterior.
}
\usage{
\S4method{bsmc}{pomp}(object, params, est, smooth = 0.1, ntries = 1,
     tol = 1e-17, lower = -Inf, upper = Inf, seed = NULL,
     verbose = getOption("verbose"), max.fail = 0, \dots)
}
\arguments{
\item{object}{
    An object of class \code{pomp} or inheriting class \code{pomp}.
  }
  \item{params}{
    A \code{npars} x \code{np} matrix containing the parameters corresponding to the initial state values in \code{xstart}.
    The matrix should be Np columns long, where Np is the number of particles.  The values for each row should be Np draws from the prior distribution for the parameter. 
    It is permissible to supply \code{params} as a named numeric vector, i.e., without a \code{dim} attribute.
    In this case, all particles will inherit the same parameter values, which is equivalent to a degenerate prior.
    If some parameter values are to be held fixed in the fitting (i.e. params.fixed = TRUE), then those elements of params (rows if params is given as a matrix), must be set to NA.
  }
  \item{est}{
    Names of the rows of \code{params} that are to be estimated.
    No updates will be made to the other parameters.
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
    Particles with log likelihood below \code{tol} are considered to be "lost".
    A filtering failure occurs when, at some time point, all particles are lost.
    When all particles are lost, the conditional log likelihood at that time point is set to be \code{log(tol)}.
  }
  \item{lower, upper}{
    optional; lower and upper bounds on the priors.
    This is useful in case there are box constraints satisfied by the priors.
    The posterior is guaranteed to lie within these bounds.
  }
  \item{seed}{
    optional; an object specifying if and how the random number generator should be initialized (\sQuote{seeded}).
    If \code{seed} is an integer, it is passed to \code{set.seed} prior to any simulation and is returned as the \dQuote{seed} element of the return list.
    By default, the state of the random number generator is not changed and the value of \code{.Random.seed} on the call is stored in the \dQuote{seed} element of the return list.
  }
  \item{verbose}{
    logical; if \code{TRUE}, print diagnostic messages.
  }
  \item{max.fail}{
    The maximum number of filtering failures allowed.
    If the number of filtering failures exceeds this number, execution will terminate with an error.
  }
  \item{\dots}{
    currently ignored.
  }
}
\value{
  A list with the following elements:
  \item{post}{
    A matrix containing draws from the approximate posterior distribution.
  }
  \item{prior}{
    A matrix containing draws from the prior distribution (identical to \code{params} on call).
  }
  \item{eff.sample.size}{
    A vector containing the effective number of particles at each time point.
  }
  \item{cond.loglik}{
    A vector containing the conditional log likelihoods at each time point.
  }
  \item{smooth}{
    The smoothing parameter used (see above).
  }
  \item{seed}{
    The state of the random number generator at the time \code{bsmc} was called.
    If the argument \code{seed} was specified, this is a copy;
    if not, this is the internal state of the random number generator at the time of call.
  }
  \item{nfail}{
    The number of filtering failures encountered.
  }
  \item{loglik}{
    The estimated log-likelihood.
  }
  \item{weights}{
    The resampling weights for each particle.
  }
}
\author{
  Michael Lavine (lavine at math dot umass dot edu),
  Matthew Ferrari (mferrari at psu dot edu),
  Aaron A. King
}
\examples{
## See the vignettes for examples.
}
\seealso{\link{pomp-class}}
\keyword{ts}
