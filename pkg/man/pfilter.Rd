\name{pfilter}
\alias{pfilter}
\alias{pfilter,pomp-method}
\alias{pfilter-pomp}
\alias{pfilter,mif-method}
\alias{pfilter-mif}
\title{Particle filter}
\description{
  Run a particle filter.
}
\usage{
pfilter(object, \dots)
\S4method{pfilter}{pomp}(object, params, Np, tol = 1e-17,
    warn = TRUE, max.fail = 0, pred.mean = FALSE, pred.var = FALSE,
    filter.mean = FALSE, .rw.sd, \dots)
\S4method{pfilter}{mif}(object, params, Np, tol = 1e-17, warn = TRUE,
    max.fail = 0, pred.mean = FALSE, pred.var = FALSE,
    filter.mean = FALSE, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp} or inheriting class \code{pomp}.
  }
  \item{params}{
    A \code{npars} x \code{np} matrix containing the parameters corresponding to the initial state values in \code{xstart}.
    This must have a 'rownames' attribute.
    It is permissible to supply \code{params} as a named numeric vector, i.e., without a \code{dim} attribute.
    In this case, all particles will inherit the same parameter values.
  }
  \item{Np}{
    Number of particles to use.
    When \code{object} is of class \code{mif}, this is by default the same number of particles used in the \code{mif} iterations.
  }
  \item{tol}{
    Particles with log likelihood below \code{tol} are considered to be "lost".
    A filtering failure occurs when, at some time point, all particles are lost.
    When all particles are lost, the conditional log likelihood at that time point is set to be \code{log(tol)}.
  }
  \item{warn}{
    Should filtering failures generate warnings?
  }
  \item{max.fail}{
    The maximum number of filtering failures allowed.
    If the number of filtering failures exceeds this number, execution will terminate with an error.
  }
  \item{pred.mean}{
    If \code{TRUE}, the prediction means are calculated for the state variables and parameters.
  }
  \item{pred.var}{
    If \code{TRUE}, the prediction variances are calculated for the state variables and parameters.
  }
  \item{filter.mean}{
    If \code{TRUE}, the filtering means are calculated for the state variables and parameters.
  }
  \item{.rw.sd}{
    For internal use with the MIF algorithm.
    If \code{TRUE}, the specified random walk SD is used.
  }
  \item{\dots}{
    Additional arguments unused at present.
  }
}
\value{
  A list with the following elements:
  \item{pred.mean}{
    The \code{nvars+npars} x \code{ntimes} matrix of prediction means, where \code{ntimes} is the length of the time series contained in \code{object}.
    The rows correspond to states and parameters, in that order.
  }
  \item{pred.variance}{
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
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\link{pomp-class}}
\keyword{models}
\keyword{ts}
