\name{pmcmc}
\docType{methods}
\alias{pmcmc}
\alias{pmcmc,pmcmc-method}
\alias{pmcmc-pmcmc}
\alias{pmcmc,pfilterd.pomp-method}
\alias{pmcmc-pfilterd.pomp}
\alias{pmcmc,pomp-method}
\alias{pmcmc-pomp}
\alias{continue,pmcmc-method}
\alias{continue-pmcmc}
\alias{pmcmc-class}
\title{The PMCMC algorithm}
\description{
  The Particle MCMC algorithm for estimating the parameters of a partially-observed Markov process.
}
\usage{
\S4method{pmcmc}{pomp}(object, Nmcmc = 1, start, pars, rw.sd, Np,
    tol = 1e-17, max.fail = 0, verbose = getOption("verbose"), \dots)
\S4method{pmcmc}{pfilterd.pomp}(object, Nmcmc = 1, Np, tol, \dots)
\S4method{pmcmc}{pmcmc}(object, Nmcmc, start, pars, rw.sd, Np,
    tol, max.fail = 0, verbose = getOption("verbose"), \dots)
\S4method{continue}{pmcmc}(object, Nmcmc = 1, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp}.
  }
  \item{Nmcmc}{
    The number of PMCMC iterations to perform.
  }
  \item{start}{
    named numeric vector;
    the starting guess of the parameters.
  }
  \item{pars}{
    optional character vector naming the ordinary parameters to be estimated.
    Every parameter named in \code{pars} must have a positive random-walk standard deviation specified in \code{rw.sd}.
    Leaving \code{pars} unspecified is equivalent to setting it equal to the names of all parameters with a positive value of \code{rw.sd}.
  }
  \item{rw.sd}{
    numeric vector with names; used to parameterize a Gaussian random walk MCMC proposal. The random walk is only applied to parameters named in \code{pars}. The algorithm requires that the random walk be nontrivial, so each element in \code{rw.sd[pars]} must be positive.
    The following must be satisfied:
    \code{names(rw.sd)} must be a subset of \code{names(start)},
    \code{rw.sd} must be non-negative (zeros are simply ignored),
    the name of every positive element of \code{rw.sd} must be in \code{pars}.
  }
  \item{Np}{
    a positive integer;
    the number of particles to use in each filtering operation.
  }
  \item{tol}{
    numeric scalar; particles with log likelihood below \code{tol} are considered to be \dQuote{lost}.
    A filtering failure occurs when, at some time point, all particles are lost.
  }
  \item{max.fail}{
    integer; maximum number of filtering failures permitted.
    If the number of failures exceeds this number, execution will terminate with an error.
  }
  \item{verbose}{
    logical; if TRUE, print progress reports.
  }
  \item{\dots}{
    Additional arguments.
    These are currently ignored.
  }
}
\value{
  An object of class \code{pmcmc}.
}
\section{Re-running PMCMC Iterations}{
  To re-run a sequence of PMCMC iterations, one can use the \code{pmcmc} method on a \code{pmcmc} object.
  By default, the same parameters used for the original PMCMC run are re-used (except for \code{tol}, \code{max.fail}, and \code{verbose}, the defaults of which are shown above).
  If one does specify additional arguments, these will override the defaults.
}
\section{Continuing PMCMC Iterations}{
  One can continue a series of PMCMC iterations from where one left off using the \code{continue} method.
  A call to \code{pmcmc} to perform \code{Nmcmc=m} iterations followed by a call to \code{continue} to perform \code{Nmcmc=n} iterations will produce precisely the same effect as a single call to \code{pmcmc} to perform \code{Nmcmc=m+n} iterations.
  By default, all the algorithmic parameters are the same as used in the original call to \code{pmcmc}.
  Additional arguments will override the defaults.
}
\section{Details}{
  \code{pmcmc} implements an MCMC algorithm in which the true likelihood of the data is replaced by an unbiased estimate computed by a particle filter.
  This gives an asymptotically correct Bayesian procedure for parameter estimation (Andrieu and Roberts, 2009).
  An extension to give a correct Bayesian posterior distribution of unobserved state variables (as in Andrieu et al, 2010) has not yet been implemented.
}
\references{
  C. Andrieu, A. Doucet and R. Holenstein, 
  Particle Markov chain Monte Carlo methods, 
  J. R. Stat. Soc. B, to appear, 2010.

  C. Andrieu and G.O. Roberts,
  The pseudo-marginal approach for efficient computation,
  Ann. Stat. 37:697-725, 2009.
}
\author{Edward L. Ionides \email{ionides at umich dot edu}, Aaron A. King \email{kingaa at umich dot edu}}
\seealso{
  \code{\link{pmcmc-class}}, \code{\link{pmcmc-methods}}, \code{\link{pomp}}, \code{\link{pomp-class}}, \code{\link{pfilter}}.
}
\keyword{ts}
