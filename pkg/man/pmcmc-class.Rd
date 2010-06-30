\name{pmcmc-class}
\docType{class}
\alias{pmcmc-class}
\title{The "pmcmc" class}
\description{
  The \code{pmcmc} class holds a fitted model and is created by a call to \code{\link{pmcmc}}.
  See \code{\link{pmcmc}} for usage.
}
\section{Objects from the Class}{
  Objects can be created by calls to the \code{\link{pmcmc}} method on an \code{\link{pomp}} object.
  Such a call uses the Particle MCMC algorithm to fit the model parameters.
}
\section{Slots}{
  A \code{pmcmc} object is derived from a \code{pomp} object and therefore has all the slots of such an object.
  See \code{\link{pomp-class}} for details.
  A full description of slots in a \code{pmcmc} object follows.
  \describe{
    \item{pars}{
      A character vector containing the names of parameters to be estimated using PMCMC.
    }
    \item{Np}{
      Number of particles to use in each filtering operation.
    }
    \item{Nmcmc}{
      Number of PMCMC iterations that have been completed.
    }
    \item{dprior}{
      A function of prototype \code{dprior(params,hyperparams,log)} that evaluates the prior density. This defaults to an improper uniform prior.
    }
    \item{hyperparams}{
      A list passed as an argument to \code{dprior}.
    }
    \item{random.walk.sd}{
      A named vector containing the random-walk variances used to parameterize a Gaussian random walk MCMC proposal.
    }
    \item{filter.mean}{
      Matrix of filtering means.
      See \code{\link{pfilter}}.
    }
    \item{eff.sample.size}{
      A vector containing the effective number of  particles at each time point.
      See \code{\link{pfilter}}.
    }
    \item{cond.loglik}{
      A vector containing the conditional log likelihoods at each time point.
      See \code{\link{pfilter}}.
    }
    \item{conv.rec}{
      The \dQuote{convergence record}: a matrix containing a record of the parameter values, log likelihoods, log prior probabilities, and other pertinent information, with one row for each PMCMC iteration.
    }
    \item{loglik}{
      A numeric value containing the value of the log likelihood, as evaluated for the random-parameter model.
      Note that this will not be equal to the log likelihood for the fixed-parameter model.
    }
    \item{log.prior}{
      A numeric value containing the log of the prior density evaluated at the parameter vector in the params slot.
    }
    \item{data, times, t0, rprocess, dprocess, dmeasure, rmeasure,
      skeleton.type, skeleton, initializer, states, params,
      statenames, paramnames, covarnames, tcovar, covar,
      PACKAGE, userdata}{
      Inherited from the \code{pomp} class.
    }
  }
}
\section{Extends}{
  Class \code{pomp}, directly.
  See \code{\link{pomp-class}}.
}
\section{Methods}{
  See \code{\link{pmcmc}}, \link{pmcmc-methods}.
}
\references{
  C. Andrieu, A. Doucet and R. Holenstein,
  Particle Markov chain Monte Carlo methods, 
  J. Roy. Stat. Soc B, to appear, 2010.

  C. Andrieu and G.O. Roberts,
  The pseudo-marginal approach for efficient computation,
  Ann Stat 37:697-725, 2009.
}
\author{Edward L. Ionides \email{ionides at umich dot edu}, Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{pmcmc}}, \link{pmcmc-methods}, \code{\link{pomp}}, \link{pomp-class}}
\keyword{models}
\keyword{ts}
