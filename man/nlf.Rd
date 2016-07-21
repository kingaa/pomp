\name{Nonlinear forecasting}
\title{Parameter estimation my maximum simulated quasi-likelihood (nonlinear forecasting)}
\alias{Nonlinear forecasting}
\alias{nlf}
\alias{nlf,ANY-method}
\alias{nlf,missing-method}
\alias{nlf,pomp-method}
\alias{nlf-pomp}
\alias{nlf,nlfd.pomp-method}
\alias{nlf-nlfd.pomp}
\alias{nlfd.pomp-class}
\alias{logLik,nlfd.pomp-method}
\alias{logLik-nlfd.pomp}
\alias{$,nlfd.pomp-method}
\alias{$-nlfd.pomp}
\description{
  \code{nlf} calls an optimizer to maximize the nonlinear forecasting (NLF) goodness of fit.
  The latter is computed by simulating data from a model, fitting a nonlinear autoregressive model to the simulated time series, and quantifying the ability of the resulting fitted model to predict the data time series.
  NLF is an \sQuote{indirect inference} method using a quasi-likelihood as the objective function.
}
\usage{
\S4method{nlf}{pomp}(object, start, est, lags, period = NA, tensor = FALSE,
     nconverge=1000, nasymp=1000, seed = 1066,
     transform.data, nrbf = 4,
     method = c("subplex", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
     skip.se = FALSE, verbose = getOption("verbose"), bootsamp = NULL,
     lql.frac = 0.1, se.par.frac = 0.1, eval.only = FALSE,
     transform = FALSE, \dots)
\S4method{nlf}{nlfd.pomp}(object, start, est, lags, period, tensor,
     nconverge, nasymp, seed, transform.data, nrbf, method,
     lql.frac, se.par.frac, transform, \dots)
}
\arguments{
  \item{object}{
    A \code{pomp} object, with the data and model to fit to it.
  }
  \item{start}{
    Named numeric vector with guessed parameters.
  }
  \item{est}{
    Vector containing the names or indices of parameters to be estimated.
  }
  \item{lags}{
    A vector specifying the lags to use when constructing the nonlinear autoregressive prediction model.
    The first lag is the prediction interval.
  }
  \item{period}{
    numeric; \code{period=NA} means the model is nonseasonal.
    period>0 is the period of seasonal forcing in 'real time'.
  }
  \item{tensor}{
    logical;
    if FALSE, the fitted model is a generalized additive model with time mod period as one of the predictors, i.e., a gam with time-varying intercept. 
    If TRUE, the fitted model is a gam with lagged state variables as predictors and time-periodic coefficients, constructed using tensor products of basis functions of state variables with basis functions of time.
  }
  \item{nconverge}{
    number of convergence timesteps to be discarded from the model simulation.
  }
  \item{nasymp}{
    number of asymptotic timesteps to be recorded from the model simulation.
  }
  \item{seed}{
    integer specifying the random number seed to use.
    When fitting, it is usually best to always run the simulations with the same sequence of random numbers, which is accomplished by setting \code{seed} to an integer.
    If you want a truly random simulation, set \code{seed=NULL}.
  }
  \item{transform}{
    logical;
    if \code{TRUE}, parameters are optimized on the transformed scale.
  }
  \item{transform.data}{
    optional function.
    If specified, forecasting is performed using data and model simulations transformed by this function.
    By default, \code{transform.data} is the identity function, i.e., no transformation is performed.
    The main purpose of \code{transform.data} is to achieve approximately multivariate normal forecasting errors.
    If data are univariate, \code{transform.data} should take a scalar and return a scalar.
    If data are multivariate, \code{transform.data} should assume a vector input and return a vector of the same length.
  }
  \item{nrbf}{
    integer scalar;
    the number of radial basis functions to be used at each lag.
  }
  \item{method}{
    Optimization method.
    Choices are \code{\link[subplex]{subplex}} and any of the methods used by \code{\link{optim}}.
  }
  \item{skip.se}{
    logical; if \code{TRUE}, skip the computation of standard errors.
  }
  \item{verbose}{
    logical; if \code{TRUE}, the negative log quasilikelihood and parameter values are printed at each iteration of the optimizer.
  }
  \item{bootsamp}{
    vector of integers;
    used to have the quasi-loglikelihood evaluated using a bootstrap re-sampling of the data set.
  } 
  \item{lql.frac}{
    target fractional change in log quasi-likelihood for quadratic standard error estimate
  } 
  \item{se.par.frac}{
    initial parameter-change fraction for quadratic standard error estimate
  } 
  \item{eval.only}{
    logical; if \code{TRUE}, no optimization is attempted and the quasi-loglikelihood value is evaluated at the \code{start} parameters.
  } 
  \item{\dots}{
    Arguments that will be passed to \code{optim} or \code{subplex} in the \code{control} list.
  }
}
\details{
  This runs an optimizer to maximize \code{nlf.objfun}.
}
\value{
  An object of class \code{nlfd.pomp}.
  \code{logLik} applied to such an object returns the log quasi likelihood.
  The \code{$} method allows extraction of arbitrary slots from the \code{nlfd.pomp} object.
}
\references{
  The following papers describe and motivate the NLF approach to model fitting:

  Ellner, S. P., Bailey, B. A., Bobashev, G. V., Gallant, A. R., Grenfell, B. T. and Nychka D. W. (1998)
  Noise and nonlinearity in measles epidemics: combining mechanistic and statistical approaches to population modeling.
  \emph{American Naturalist} \bold{151}, 425--440.

  Kendall, B. E., Briggs, C. J., Murdoch, W. W., Turchin, P., Ellner, S. P., McCauley, E., Nisbet, R. M. and Wood S. N. (1999)
  Why do populations cycle? A synthesis of statistical and mechanistic modeling approaches.
  \emph{Ecology} \bold{80}, 1789--1805.

  Kendall, B. E., Ellner, S. P., McCauley, E., Wood, S. N., Briggs, C. J., Murdoch, W. W. and Turchin, P. (2005)
  Population cycles in the pine looper moth (\emph{Bupalus piniarius}): dynamical tests of mechanistic hypotheses.
  \emph{Ecological Monographs} \bold{75}, 259--276.
}
\author{Stephen P. Ellner, Bruce E. Kendall, Aaron A. King}
\keyword{optimize}
\keyword{ts}
