\name{pfilter}
\alias{pfilter}
\alias{pfilter,pomp-method}
\alias{pfilter-pomp}
\alias{pfilter,mif-method}
\alias{pfilter-mif}
\title{Particle filter}
\description{
  Run a plain vanilla particle filter.
  Resampling is performed after each observation.
}
\usage{
pfilter(object, \dots)
\S4method{pfilter}{pomp}(object, params, Np, tol = 1e-17,
    max.fail = 0, pred.mean = FALSE, pred.var = FALSE,
    filter.mean = FALSE, save.states = FALSE,
    verbose = getOption("verbose"), \dots)
\S4method{pfilter}{mif}(object, params, Np, tol = 1e-17,
    max.fail = 0, pred.mean = FALSE, pred.var = FALSE,
    filter.mean = FALSE, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp} or inheriting class \code{pomp}.
  }
  \item{params}{
    A \code{npars} x \code{Np} matrix containing the parameters corresponding to the initial state values in \code{xstart}.
    This must have a 'rownames' attribute.
    It is permissible to supply \code{params} as a named numeric vector, i.e., without a \code{dim} attribute.
    In this case, all particles will inherit the same parameter values.
  }
  \item{Np}{
    integer; number of particles to use.
    When \code{object} is of class \code{mif}, this is by default the same number of particles used in the \code{mif} iterations.
  }
  \item{tol}{
    positive numeric scalar; particles with log likelihood below \code{tol} are considered to be "lost".
    A filtering failure occurs when, at some time point, all particles are lost.
    When all particles are lost, the conditional log likelihood at that time point is set to be \code{log(tol)}.
  }
  \item{max.fail}{
    integer; the maximum number of filtering failures allowed.
    If the number of filtering failures exceeds this number, execution will terminate with an error.
  }
  \item{pred.mean}{
    logical; if \code{TRUE}, the prediction means are calculated for the state variables and parameters.
  }
  \item{pred.var}{
    logical; if \code{TRUE}, the prediction variances are calculated for the state variables and parameters.
  }
  \item{filter.mean}{
    logical; if \code{TRUE}, the filtering means are calculated for the state variables and parameters.
  }
  \item{save.states}{
    logical; if \code{TRUE}, the state-vector for each particle is saved and returned.
  }
  \item{verbose}{
    logical; if \code{TRUE}, progress information is reported as \code{pfilter} works.
  }
  \item{\dots}{
    Additional arguments unused at present.
  }
}
\value{
  A list with the following elements:
  \item{pred.mean}{
    The matrix of prediction means.
    The rows correspond to states and parameters (if appropriate), in that order, the columns to successive observations in the time series contained in \code{object}.
  }
  \item{pred.var}{
    The matrix of prediction variances, in the same format as \code{pred.mean}.
  }
  \item{filter.mean}{
    The matrix of filtering means, in the same format as \code{pred.mean}.
  }
  \item{eff.sample.size}{
    A vector containing the effective number of particles at each time point.
  }
  \item{cond.loglik}{
    A vector containing the conditional log likelihoods at each time point.
  }
  \item{states}{
    If \code{saves.states=TRUE}, the array of state-vectors at each time point, for each particle.
    An array with dimensions \code{nvars}-by-\code{Np}-by-\code{ntimes}.
    In particular, \code{states[,i,t]} can be considered a sample from \eqn{f[X|y_{1:t}]}.
  }
  \item{nfail}{
    The number of filtering failures encountered.
  }
  \item{loglik}{
    The estimated log-likelihood.
  }
}
\examples{
## See the vignettes for examples.
}
\references{
  M. S. Arulampalam, S. Maskell, N. Gordon, & T. Clapp.
  A Tutorial on Particle Filters for Online Nonlinear, Non-Gaussian Bayesian Tracking.
  IEEE Trans. Sig. Proc. 50:174--188, 2002.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\link{pomp-class}}
\keyword{models}
\keyword{ts}
