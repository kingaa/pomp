\name{spect}
\alias{spect}
\alias{spect,pomp-method}
\alias{spect-pomp}
\alias{spect,spect.pomp-method}
\alias{spect-spect.pomp}
\alias{spect.match}
\alias{spect.match,pomp-method}
\alias{spect.match-pomp}
\alias{spect.match,spect.pomp-method}
\alias{spect.match-spect.pomp}
\title{Power spectrum computation for partially-observed Markov processes.}
\description{
  \code{spect} estimates the power spectrum of time series data and model simulations and compares the results.
  It can be used to diagnose goodness of fit and/or as the basis for frequency-domain parameter estimation (\code{spect.match}).

  \code{spect.match} tries to match the power spectrum of the model to that of the data.
  It calls an optimizer to adjust model parameters to minimize the discrepancy between simulated and actual data.
}
\usage{
  \S4method{spect}{pomp}(object, params, vars, kernel.width, nsim, seed = NULL,
      transform = identity,
      detrend = c("none","mean","linear","quadratic"),
      \dots)
  \S4method{spect}{spect.pomp}(object, params, vars, kernel.width, nsim, seed = NULL, transform,
      detrend, \dots)
  spect.match(object, start, est = character(0),
              vars, nsim, seed = NULL,
              kernel.width, transform = identity, 
              detrend = c("none","mean","linear","quadratic"),
              weights, method = c("subplex","Nelder-Mead","SANN"),
              verbose = getOption("verbose"),
              eval.only = FALSE, fail.value = NA, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp}.
  }
  \item{params}{
    optional named numeric vector of model parameters.
    By default, \code{params=coef(object)}.
  }
  \item{vars}{
    optional; names of observed variables for which the power spectrum will be computed.
    This must be a subset of \code{rownames(data.array(object))}.
    By default, the spectrum will be computed for all observables.
  }
  \item{kernel.width}{
    width parameter for the smoothing kernel used for calculating the estimate of the spectrum.
  }
  \item{nsim}{
    number of model simulations to be computed.
  }
  \item{seed}{
    optional; if non-\code{NULL}, the random number generator will be initialized with this seed for simulations.
    See \link{simulate-pomp}.
  }
  \item{transform}{
    function; this transformation will be applied to the observables prior to estimation of the spectrum, and prior to any detrending.
  }
  \item{detrend}{
    de-trending operation to perform.
    Options include no detrending, and subtraction of constant, linear, and quadratic trends from the data.
    Detrending is applied to each data series and to each model simulation independently.
  }
  \item{weights}{
    optional.
    The mismatch between model and data is measured by a weighted average of mismatch at each frequency.
    By default, all frequencies are weighted equally.
    \code{weights} can be specified either as a vector (which must have length equal to the number of frequencies) or as a function of frequency.
    If the latter, \code{weights(freq)} must return a nonnegative weight for each frequency.
  }
  \item{start}{
    named numeric vector; the initial guess of parameters.
  }
  \item{est}{
    character vector; the names of parameters to be estimated.
  }
  \item{method}{
    Optimization method.
    Choices are \code{\link[subplex]{subplex}} and any of the methods used by \code{\link{optim}}.
  }
  \item{verbose}{
    logical; print diagnostic messages?
  }
  \item{eval.only}{
    logical; if \code{TRUE}, no optimization is attempted.
    Instead, the probe-mismatch value is simply evaluated at the \code{start} parameters.
  }
  \item{fail.value}{
    optional scalar; if non-\code{NA}, this value is substituted for non-finite values of the objective function.
  }
  \item{\dots}{
    Additional arguments.
    In the case of \code{spec}, these are currently ignored.
    In the case of \code{spec.match}, these are passed to \code{optim} or \code{subplex} in the \code{control} list.
  }
}
\details{
  A call to \code{spect} results in the estimation of the power spectrum for the (transformed, detrended) data and \code{nsim} model simulations.
  The results of these computations are stored in an object of class \code{\link[=spect.pomp-class]{spect.pomp}}.


  A call to \code{spect.match} results in an attempt to optimize the agreement between model and data spectrum over the parameters named in \code{est}.
  The results, including coefficients of the fitted model and power spectra of fitted model and data, are stored in an object of class \code{\link[=spect.matched.pomp-class]{spect.matched.pomp}}.
}
\value{
  \code{spect} returns an object of class \code{\link[=spect.pomp-class]{spect.pomp}}.
  \code{spect.match} returns an object of class \code{\link[=spect.matched.pomp-class]{spect.matched.pomp}}.
}
\references{
  D.C. Reuman, R.A. Desharnais, R.F. Costantino, O. Ahmad, J.E. Cohen (2006)
  Power spectra reveal the influence of stochasticity on nonlinear population dynamics.
  \emph{Proceedings of the National Academy of Sciences} \bold{103}, 18860-18865.
  
  D.C. Reuman, R.F. Costantino, R.A. Desharnais, J.E. Cohen (2008)
  Color of environmental noise affects the nonlinear dynamics of cycling, stage-structured populations.
  \emph{Ecology Letters}, \bold{11}, 820-830.
}
\author{
  Daniel C. Reuman (d.reuman at imperial dot ac dot uk)

  Cai GoGwilt
  
  Aaron A. King (kingaa at umich dot edu)
}
\examples{
data(ou2)
good <- spect(
              ou2,
              vars=c("y1","y2"),
              kernel.width=3,
              detrend="mean",
              nsim=500
              )
summary(good)
plot(good)

ou2.bad <- ou2
coef(ou2.bad,c("x1.0","x2.0","alpha.1","alpha.4")) <- c(0,0,0.1,0.2)
bad <- spect(
             ou2.bad,
             vars=c("y1","y2"),
             kernel.width=3,
             detrend="mean",
             nsim=500
             )
summary(bad)
plot(bad)
}
\seealso{
  \link{pomp-class},
  \link{pomp-methods},
  \code{\link{probe}},
  \code{\link{probe.match}}
}
\keyword{models}
\keyword{ts}




