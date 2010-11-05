\name{pmcmc-class}
\docType{class}
\keyword{internal}
\alias{pmcmc-class}
\title{The "pmcmc" class}
\description{
  The \code{pmcmc} class holds a fitted model and is created by a call to \code{\link{pmcmc}}.
  See \code{\link{pmcmc}} for usage.
}
\section{Objects from the Class}{
  Objects can be created by calls to the \code{\link{pmcmc}} method on a \code{\link{pomp}} or \code{pfilterd.pomp} object.
  Such a call uses the Particle MCMC algorithm to fit the model parameters.
}
\section{Slots}{
  A \code{pmcmc} object is derived from a \code{pfilterd.pomp} object and therefore has all the slots of such an object.
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
