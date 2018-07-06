\name{Probes and synthetic likelihood}
\title{Probe a partially-observed Markov process by computing summary statistics and the synthetic likelihood.}
\alias{Probes and synthetic likelihood}
\alias{probe}
\alias{probe,missing-method}
\alias{probe,ANY-method}
\alias{probe,pomp-method}
\alias{probe,probed.pomp-method}
\alias{probe-pomp}
\alias{probe-probed.pomp}
\alias{probe.match}
\alias{probe.match,ANY-method}
\alias{probe.match,missing-method}
\alias{probe.match,pomp-method}
\alias{probe.match,probed.pomp-method}
\alias{probe.match,probe.matched.pomp-method}
\alias{probe.match-pomp}
\alias{probe.match-probed.pomp}
\alias{probe.match-probe.matched.pomp}
\alias{probed.pomp-class}
\alias{probe.matched.pomp-class}
\alias{probe.match.objfun}
\alias{probe.match.objfun,ANY-method}
\alias{probe.match.objfun,missing-method}
\alias{probe.match.objfun,pomp-method}
\alias{probe.match.objfun,probed.pomp-method}
\alias{probe.match.objfun-pomp}
\alias{probe.match.objfun-probed.pomp}
\alias{probed.pomp-methods}
\alias{probe.matched.pomp-methods}
\alias{spect.pomp-methods}
\alias{spect.matched.pomp-methods}
\alias{summary,probed.pomp-method}
\alias{summary,probe.matched.pomp-method}
\alias{summary-probed.pomp}
\alias{summary-probe.matched.pomp}
\alias{values}
\alias{values,probed.pomp-method}
\alias{values,probe.matched.pomp-method}
\alias{values-probed.pomp}
\alias{values-probe.matched.pomp}
\alias{plot,probed.pomp-method}
\alias{plot,probe.matched.pomp-method}
\alias{plot-probed.pomp}
\alias{plot-probe.matched.pomp}
\alias{summary,spect.pomp-method}
\alias{summary,spect.matched.pomp-method}
\alias{summary-spect.pomp}
\alias{summary-spect.matched.pomp}
\alias{plot,spect.pomp-method}
\alias{plot,spect.matched.pomp-method}
\alias{plot,spect.pomp-method}
\alias{plot-spect.pomp}
\alias{as,probed.pomp-method}
\alias{logLik,probed.pomp-method}
\alias{logLik-probed.pomp}
\alias{$,probed.pomp-method}
\alias{$-probed.pomp}
\alias{$,probe.matched.pomp-method}
\alias{$-probe.matched.pomp}
\alias{coerce,probed.pomp,data.frame-method}
\description{
  \code{probe} applies one or more \dQuote{probes} to time series data and model simulations and compares the results.
  It can be used to diagnose goodness of fit and/or as the basis for \dQuote{probe-matching}, a generalized method-of-moments approach to parameter estimation.
  \code{probe.match} calls an optimizer to adjust model parameters to do probe-matching, i.e., to minimize the discrepancy between simulated and actual data.
  This discrepancy is measured using the \dQuote{synthetic likelihood} as defined by Wood (2010).
  \code{probe.match.objfun} constructs an objective function for probe-matching suitable for use in \code{optim}-like optimizers.
}
\usage{
\S4method{probe}{pomp}(object, probes, params, nsim, seed = NULL, \dots)
\S4method{probe}{probed.pomp}(object, probes, nsim, seed, \dots)
\S4method{probe.match.objfun}{pomp}(object, params, est, probes, nsim,
            seed = NULL, fail.value = NA, transform = FALSE, \dots)
\S4method{probe.match.objfun}{probed.pomp}(object, probes, nsim, seed, \dots)
\S4method{probe.match}{pomp}(object, start, est = character(0),
            probes, nsim, seed = NULL,
            method = c("subplex","Nelder-Mead","SANN","BFGS",
                       "sannbox","nloptr"),
            verbose = getOption("verbose"),
            fail.value = NA, transform = FALSE, \dots)
\S4method{probe.match}{probed.pomp}(object, probes, nsim, seed,
            \dots, verbose = getOption("verbose"))
\S4method{probe.match}{probe.matched.pomp}(object, est, probes,
            nsim, seed, transform, fail.value, \dots,
            verbose = getOption("verbose"))
\S4method{logLik}{probed.pomp}(object, \dots)
\S4method{values}{probed.pomp}(object, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp}.
  }
  \item{probes}{
    A single probe or a list of one or more probes.
    A probe is simply a scalar- or vector-valued function of one argument that can be applied to the data array of a \code{pomp}.
    A vector-valued probe must always return a vector of the same size.
    A number of useful examples are provided with the package: see \link{probe functions}).
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
  \item{method}{
    Optimization method.
    Choices refer to algorithms used in \code{\link{optim}}, \code{\link[subplex]{subplex}}, and \code{\link[nloptr]{nloptr}}.
  }
  \item{verbose}{
    logical; print diagnostic messages?
  }
  \item{fail.value}{
    optional numeric scalar;
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
    In the case of \code{probe.match}, these are passed to the optimizer (one of \code{\link{optim}}, \code{\link[subplex]{subplex}}, \code{\link[nloptr]{nloptr}}, or \code{\link{sannbox}}).
    These are passed via the optimizer's \code{control} list (in the case of \code{optim}, \code{subplex}, and \code{sannbox}) or the \code{opts} list (in the case of \code{nloptr}).
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
    \item{est, transform, fail.value}{values of the corresponding arguments in the call to \code{probe.match}.}
    \item{value}{value of the objective function at the optimum.}
    \item{evals}{
      number of function and gradient evaluations by the optimizer.
      See \code{\link{optim}}.
    }
    \item{convergence, msg}{
      Convergence code and message from the optimizer.
      See \code{\link{optim}} and \code{\link{nloptr}}.
    }
  }

  \code{probe.match.objfun} returns a function suitable for use as an objective function in an \code{\link{optim}}-like optimizer.
}
\section{Methods}{
  \describe{
    \item{plot}{
      displays diagnostic plots.
    }
    \item{summary}{
      displays summary information.
    }
    \item{values}{
      extracts the realized values of the probes on the data and on the simulations as a data frame in long format.
      The variable \code{.id} indicates whether the probes are from the data or simulations.
    }
    \item{logLik}{
      returns the synthetic likelihood for the probes.
      NB: in general, this is not the same as the likelihood.
    }
    \item{as, as.data.frame}{
      when a \sQuote{probed.pomp} is coerced to a \sQuote{data.frame}, the first row gives the probes applied to the data; the rest of the rows give the probes evaluated on simulated data.
      The rownames of the result can be used to distinguish these.
    }
  }
  In addition, slots of this object can be accessed via the \code{$} operator.
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
  Daniel C. Reuman, Aaron A. King
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
              nsim=300
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
             nsim=300
             )
summary(bad)
plot(bad)
}
\seealso{
  \link{pomp}, \link{probe functions}, \link{spect},
  and the tutorials on the \href{https://kingaa.github.io/pomp}{package website}.
}
\keyword{optimize}
\keyword{ts}
