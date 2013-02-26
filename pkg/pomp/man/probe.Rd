\name{probe}
\alias{probe}
\alias{probe,pomp-method}
\alias{probe,probed.pomp-method}
\alias{probe-pomp}
\alias{probe-probed.pomp}
\alias{probe.match}
\alias{probe.match,pomp-method}
\alias{probe.match,probed.pomp-method}
\alias{probe.match,probe.matched.pomp-method}
\alias{probe.match-pomp}
\alias{probe.match-probed.pomp}
\alias{probe.match-probe.matched.pomp}
\alias{probed.pomp-class}
\alias{probe.matched.pomp-class}
\alias{probe.match.objfun}
\title{Probe a partially-observed Markov process.}
\description{
  \code{probe} applies one or more \dQuote{probes} to time series data and model simulations and compares the results.
  It can be used to diagnose goodness of fit and/or as the basis for \dQuote{probe-matching}, a generalized method-of-moments approach to parameter estimation.
  \code{probe.match} calls an optimizer to adjust model parameters to do probe-matching, i.e., to minimize the discrepancy between simulated and actual data.
  This discrepancy is measured using the \dQuote{synthetic likelihood} as defined by Wood (2010).
  \code{probe.match.objfun} constructs an objective function for probe-matching suitable for use in \code{optim}-like optimizers.
}
\usage{
\S4method{probe}{pomp}(object, probes, params, nsim, seed = NULL, \dots)
\S4method{probe}{probed.pomp}(object, probes, params, nsim, \dots)
\S4method{probe.match}{pomp}(object, start, est = character(0),
            probes, weights,
            nsim, seed = NULL,
            method = c("subplex","Nelder-Mead","SANN","BFGS","sannbox"),
            verbose = getOption("verbose"), 
            eval.only = FALSE, fail.value = NA, transform = FALSE, \dots)
\S4method{probe.match}{probed.pomp}(object, start, est = character(0),
            probes, weights,
            nsim, seed = NULL,
            method = c("subplex","Nelder-Mead","SANN","BFGS","sannbox"),
            verbose = getOption("verbose"), 
            eval.only = FALSE, fail.value = NA, transform = FALSE, \dots)
\S4method{probe.match}{probe.matched.pomp}(object, start, est,
            probes, weights,
            nsim, seed = NULL,
            method = c("subplex","Nelder-Mead","SANN","BFGS","sannbox"),
            verbose = getOption("verbose"), 
            eval.only = FALSE, fail.value, transform, \dots)
probe.match.objfun(object, params, est, probes, nsim = 1,
                   seed = NULL, fail.value = NA, transform = FALSE, \dots)

}
\arguments{
  \item{object}{
    An object of class \code{pomp}.
  }
  \item{probes}{
    A single probe or a list of one or more probes.
    A probe is simply a scalar- or vector-valued function of one argument that can be applied to the data array of a \code{pomp}.
    A vector-valued probe must always return a vector of the same size.
    A number of basic examples are provided with the package (see \link{basic.probes}).
  }
  \item{params}{
    optional named numeric vector of model parameters.
    By default, \code{params=coef(object)}.
  }
  \item{nsim}{
    The number of model simulations to be computed.
  }
  \item{seed}{
    optional; if non-\code{NULL}, the random number generator will be initialized with this seed for simulations.
    See \link{simulate-pomp}.
  }
  \item{start}{
    named numeric vector; the initial guess of parameters.
  }
  \item{est}{
    character vector; the names of parameters to be estimated.
  }
  \item{weights}{
    optional numeric vector of relative weights.
    Must be of the same length as \code{probes}.
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
    optional scalar;
    if non-\code{NA}, this value is substituted for non-finite values of the objective function.
    It should be a large number (i.e., bigger than any legitimate values the objective function is likely to take).
  }
  \item{transform}{
    logical;
    if \code{TRUE}, optimization is performed on the transformed scale.
  }
  \item{\dots}{
    Additional arguments.
    In the case of \code{probe}, these are currently ignored.
    In the case of \code{probe.match}, these are passed to \code{optim} or \code{subplex} in the \code{control} list.
  }
}
\details{
  A call to \code{probe} results in the evaluation of the probe(s) in \code{probes} on the data.
  Additionally, \code{nsim} simulated data sets are generated (via a call to \code{\link[=simulate-pomp]{simulate}}) and the probe(s) are applied to each of these.
  The results of the probe computations on real and simulated data are stored in an object of class \code{probed.pomp}.

  A call to \code{probe.match} results in an attempt to optimize the agreement between model and data, as measured by the specified probes, over the parameters named in \code{est}.
  The results, including coefficients of the fitted model and values of the probes for data and fitted-model simulations, are stored in an object of class \code{\link[=probe.matched.pomp-class]{probe.matched.pomp}}.

  The objective function minimized by \code{probe.match} --- in a form suitable for use with \code{\link{optim}}-like optimizers --- is created by a call to \code{probe.match.objfun}.
  Specifically, \code{probe.match.objfun} will return a function that takes a single numeric-vector argument that is assumed to cotain the parameters named in \code{est}, in that order.
  This function will return the negative synthetic log likelihood for the probes specified.
}
\value{
  \code{probe} returns an object of class \code{probed.pomp}.
  \code{probed.pomp} is derived from the \code{\link[=pomp-class]{pomp}} class and therefore have all the slots of \code{pomp}.
  In addition, a \code{probed.pomp} class has the following slots:
  \describe{
    \item{probes}{list of the probes applied.}
    \item{datvals, simvals}{values of each of the probes applied to the real and simulated data, respectively.}
    \item{quantiles}{fraction of simulations with probe values less than the value of the probe of the data.}
    \item{pvals}{
      two-sided p-values:
      fraction of the \code{simvals} that deviate more extremely from the mean of the \code{simvals} than does \code{datavals}.
    }
    \item{synth.loglik}{
      the log synthetic likelihood (Wood 2010).
      This is the likelihood assuming that the probes are multivariate-normally distributed.
    }
  }

  \code{probe.match} returns an object of class \code{probe.matched.pomp}, which is derived from class \code{probed.pomp}.
  \code{probe.matched.pomp} objects therefore have all the slots above plus the following:
  \describe{
    \item{est, weights, fail.value}{values of the corresponding arguments in the call to \code{spect.match}.}
    \item{value}{value of the objective function.}
    \item{evals}{
      number of function and gradient evaluations by the optimizer.
      See \code{\link{optim}}.
    }
    \item{convergence, msg}{
      Convergence code and message from the optimizer.
      See \code{\link{optim}}.
    }
  }

  \code{probe.match.objfun} returns a function suitable for use as an objective function in an \code{\link{optim}}-like optimizer.
}
\references{
  B. E. Kendall, C. J. Briggs, W. M. Murdoch, P. Turchin, S. P. Ellner, E. McCauley, R. M. Nisbet, S. N. Wood
  Why do populations cycle? A synthesis of statistical and mechanistic modeling approaches,
  Ecology, 80:1789--1805, 1999.

  S. N. Wood
  Statistical inference for noisy nonlinear ecological dynamic systems,
  Nature, 466: 1102--1104, 2010.
}
\author{
  Daniel C. Reuman (d.reuman at imperial dot ac dot uk)
  
  Aaron A. King (kingaa at umich dot edu)
}
\examples{
pompExample(ou2)
good <- probe(
              ou2,
              probes=list(
                y1.mean=probe.mean(var="y1"),
                y2.mean=probe.mean(var="y2"),
                y1.sd=probe.sd(var="y1"),
                y2.sd=probe.sd(var="y2"),
                extra=function(x)max(x["y1",])
                ),
              nsim=500
              )
summary(good)
plot(good)

bad <- probe(
             ou2,
             params=c(alpha.1=0.1,alpha.4=0.2,x1.0=0,x2.0=0,
                      alpha.2=-0.5,alpha.3=0.3,
                      sigma.1=3,sigma.2=-0.5,sigma.3=2,
                      tau=1),
             probes=list(
               y1.mean=probe.mean(var="y1"),
               y2.mean=probe.mean(var="y2"),
               y1.sd=probe.sd(var="y1"),
               y2.sd=probe.sd(var="y2"),
               extra=function(x)range(x["y1",])
               ),
             nsim=500
             )
summary(bad)
plot(bad)
}
\seealso{
  \link{pomp-class},
  \link{pomp-methods},
  \link{basic.probes},
  \code{\link{probe.match}}
}
\keyword{models}
\keyword{ts}




