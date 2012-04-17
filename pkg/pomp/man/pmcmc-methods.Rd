\name{pmcmc-methods}
\docType{methods}
\alias{pmcmc-methods}
\alias{dprior}
\alias{dprior,pmcmc-method}
\alias{dprior-pmcmc}
\alias{logLik,pmcmc-method}
\alias{logLik-pmcmc}
\alias{conv.rec,pmcmc-method}
\alias{conv.rec-pmcmc}
\alias{filter.mean,pmcmc-method}
\alias{filter.mean-pmcmc}
\alias{plot-pmcmc}
\alias{plot,pmcmc-method}
\alias{compare.pmcmc}
\title{Methods of the "pmcmc" class}
\description{Methods of the "pmcmc" class.}
\usage{
\S4method{logLik}{pmcmc}(object, \dots)
\S4method{conv.rec}{pmcmc}(object, pars, \dots)
\S4method{filter.mean}{pmcmc}(object, pars, \dots)
\S4method{plot}{pmcmc}(x, y = NULL, \dots)
\S4method{dprior}{pmcmc}(object, params, log = FALSE, \dots)
compare.pmcmc(z)
}
\arguments{
  \item{object, x}{The \code{pmcmc} object.}
  \item{pars}{Names of parameters.}
  \item{y}{Ignored.}
  \item{z}{A \code{pmcmc} object or list of \code{pmcmc} objects.}
  \item{params}{
    Named vector of parameters.
  }
  \item{log}{if TRUE, log probabilities are returned.}
  \item{\dots}{
    Further arguments (either ignored or passed to underlying functions).
  }
}
\section{Methods}{
  \describe{
    \item{conv.rec}{
      \code{conv.rec(object, pars = NULL)} returns the columns of the convergence-record matrix corresponding to the names in \code{pars}.
      By default, all rows are returned.
    }
    \item{logLik}{
      Returns the value in the \code{loglik} slot.
    }
    \item{dprior}{
       \code{dprior(object,params,log)} evaluates the prior density at \code{params} with values of the hyperparameters given by \code{object@hyperparams}.
    } 
    \item{pmcmc}{
      Re-runs the PMCMC iterations.
      See the documentation for \code{\link{pmcmc}}.
    }
    \item{compare.pmcmc}{
      Given a \code{pmcmc} object or a list of \code{pmcmc} objects, \code{compare.pmcmc} produces a set of diagnostic plots.
    }
    \item{plot}{
      Plots a series of diagnostic plots.
      When \code{x} is a \code{pmcmc} object, \code{plot(x)} is equivalent to \code{compare.pmcmc(list(x))}.
    }
    \item{filter.mean}{
      \code{filter.mean(object, pars = NULL)} returns the rows of the filtering-mean matrix corresponding to the names in \code{pars}.
      By default, all rows are returned.
    }
    \item{print}{
      Prints a summary of the \code{pmcmc} object.
    }
    \item{show}{
      Displays the \code{pmcmc} object.
    }
    \item{pfilter}{
      See \code{\link{pfilter}}.
    }
  }
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
\seealso{\code{\link{pmcmc}}, \code{\link{pomp}}, \code{\link{pomp-class}}, \code{\link{pfilter}}}
\keyword{models}
\keyword{ts}
