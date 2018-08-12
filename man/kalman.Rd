\name{Ensemble Kalman filters}
\title{Ensemble Kalman filters}
\alias{ensembe Kalman filter}
\alias{ensemble adjustment Kalman filter}
\alias{enkf}
\alias{enkf,ANY-method}
\alias{eakf,ANY-method}
\alias{enkf,missing-method}
\alias{eakf,missing-method}
\alias{enkf,pomp-method}
\alias{enkf-pomp}
\alias{eakf}
\alias{eakf,pomp-method}
\alias{eakf-pomp}
\alias{kalmand_pomp}
\alias{logLik,kalmand_pomp-method}
\alias{logLik-kalmand_pomp}
\alias{pred.mean,kalmand_pomp-method}
\alias{pred.mean-kalmand_pomp}
\alias{filter.mean,kalmand_pomp-method}
\alias{filter.mean-kalmand_pomp}
\alias{cond.logLik,kalmand_pomp-method}
\alias{cond.logLik-kalmand_pomp}
\alias{forecast}
\alias{forecast,kalmand_pomp-method}
\alias{forecast-kalmand_pomp}
\alias{as,kalmand_pomp-method}
\alias{as.data.frame.kalmand_pomp}
\alias{coerce,kalmand_pomp,data.frame-method}
\description{
  The ensemble Kalman filter and ensemble adjustment Kalman filter.
}
\usage{
\S4method{enkf}{pomp}(object, params, Np, h, R, \dots,
              verbose = getOption("verbose", FALSE))
\S4method{eakf}{pomp}(object, params, Np, C, R, \dots,
              verbose = getOption("verbose", FALSE))
\S4method{logLik}{kalmand_pomp}(object, \dots)
\S4method{cond.logLik}{kalmand_pomp}(object, \dots)
\S4method{pred.mean}{kalmand_pomp}(object, vars, \dots)
\S4method{filter.mean}{kalmand_pomp}(object, vars, \dots)
\S4method{forecast}{kalmand_pomp}(object, vars, \dots)
}
\arguments{
  \item{object}{
    An object of class \sQuote{pomp} or inheriting class \sQuote{pomp}.
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
  \item{vars}{Names of variables.}
  \item{\dots}{
    additional arguments (currently ignored).
  }
}
\value{
  An object of class \sQuote{kalmand_pomp}.
  This class inherits from class \sQuote{pomp}.
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
    \item{forecast}{
      Extract the mean of the forecasting distribution, which is that of \deqn{Y_t \vert y_1,\dots,y_{t-1},}{Y_t | y_1,\dots,y_{t-1},}.
    }
  }
}
\author{Aaron A. King}
\seealso{
  \code{\link{pomp}}, \code{\link{pfilter}},
  and the tutorials on the \href{https://kingaa.github.io/pomp}{package website}.
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
