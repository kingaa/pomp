\name{Ensemble Kalman filters}
\title{Ensemble Kalman filters}
\alias{ensembe Kalman filter}
\alias{ensemble adjustment Kalman filter}
\alias{enkf}
\alias{enkf,pomp-method}
\alias{enkf-pomp}
\alias{eakf}
\alias{eakf,pomp-method}
\alias{eakf-pomp}
\alias{kalmand.pomp-class}
\alias{kalmand.pomp}
\alias{logLik,kalmand.pomp-method}
\alias{logLik-kalmand.pomp}
\alias{$,kalmand.pomp-method}
\alias{$-kalmand.pomp}
\alias{pred.mean,kalmand.pomp-method}
\alias{pred.mean-kalmand.pomp}
\alias{filter.mean,kalmand.pomp-method}
\alias{filter.mean-kalmand.pomp}
\alias{cond.logLik,kalmand.pomp-method}
\alias{cond.logLik-kalmand.pomp}
\alias{as,kalmand.pomp-method}
\alias{as.data.frame.kalmand.pomp}
\alias{coerce,kalmand.pomp,data.frame-method}
\description{
  The ensemble Kalman filter and ensemble adjustment Kalman filter.
}
\usage{
\S4method{enkf}{pomp}(object, params, Np, h, R,
              verbose = getOption("verbose"), \dots)
\S4method{eakf}{pomp}(object, params, Np, C, R,
              verbose = getOption("verbose"), \dots)
\S4method{logLik}{kalmand.pomp}(object, \dots)
\S4method{cond.logLik}{kalmand.pomp}(object, \dots)
\S4method{pred.mean}{kalmand.pomp}(object, pars, \dots)
\S4method{filter.mean}{kalmand.pomp}(object, pars, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp} or inheriting class \code{pomp}.
  }
  \item{params}{
    optional named numeric vector containing the parameters at which the filtering should be performed.
    By default, \code{params = coef(object)}.
  }
  \item{Np}{
    the number of particles to use.
  }
  \item{verbose}{
    logical; if \code{TRUE}, progress information is reported.
  }
  \item{h}{
    function returning the expected value of the observation given the state.
  }
  \item{C}{
    matrix converting state vector into expected value of the observation.
  }
  \item{R}{
    matrix; variance of the measurement noise.
  }
  \item{pars}{Names of variables.}
  \item{\dots}{
    additional arguments (currently ignored).
  }
}
\value{
  An object of class \code{kalmand.pomp}.
  This class inherits from class \code{\link[=pomp-class]{pomp}}.
}
\section{Methods}{
  \describe{
    \item{logLik}{
      Extracts the estimated log likelihood.
    }
    \item{cond.logLik}{
      Extracts the estimated conditional log likelihood
      \deqn{\ell_t(\theta) = \mathrm{Prob}[y_t \vert y_1, \dots, y_{t-1}],}{ell_t(theta)=Prob[y_t | y_1, \dots, y_(t-1)],}
      where \eqn{y_t} are the data, at time \eqn{t}.
    }
    \item{pred.mean}{
      Extract the mean of the approximate prediction distribution.
      This prediction distribution is that of \deqn{X_t \vert y_1,\dots,y_{t-1},}{X_t | y_1,\dots,y_(t-1),} where \eqn{X_t}, \eqn{y_t} are the state vector and data, respectively, at time \eqn{t}.
    }
    \item{filter.mean}{
      Extract the mean of the filtering distribution, which is that of \deqn{X_t \vert y_1,\dots,y_t,}{X_t | y_1,\dots,y_t,} where \eqn{X_t}, \eqn{y_t} are the state vector and data, respectively, at time \eqn{t}.
    }
  }
}
\author{Aaron A. King}
\seealso{
  \code{\link{pomp}}, \code{\link{pfilter}},
  and the tutorials on the \href{http://kingaa.github.io/pomp}{package website}.
}
\references{
  Evensen, G. (1994)
  Sequential data assimilation with a nonlinear quasi-geostrophic model using Monte Carlo methods to forecast error statistics
  Journal of Geophysical Research: Oceans 99:10143--10162

  Evensen, G. (2009)
  Data assimilation: the ensemble Kalman filter
  Springer-Verlag.

  Anderson, J. L. (2001)
  An Ensemble Adjustment Kalman Filter for Data Assimilation
  Monthly Weather Review 129:2884--2903
}
\keyword{ts}
\keyword{internals}
