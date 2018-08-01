\name{Probes and synthetic likelihood}
\title{Probe a partially-observed Markov process by computing summary statistics and the synthetic likelihood.}
\alias{Probes and synthetic likelihood}
\alias{probe}
\alias{probevals}
\alias{probe,missing-method}
\alias{probe,ANY-method}
\alias{probe,pomp-method}
\alias{probe,probed_pomp-method}
\alias{probe-pomp}
\alias{probe-probed_pomp}
\alias{probe.match}
\alias{probe.match,ANY-method}
\alias{probe.match,missing-method}
\alias{probe.match,pomp-method}
\alias{probe.match,probed_pomp-method}
\alias{probe.match,probe_matched_pomp-method}
\alias{probe.match-pomp}
\alias{probe.match-probed_pomp}
\alias{probe.match-probe_matched_pomp}
\alias{probe.match.objfun}
\alias{probe.match.objfun,ANY-method}
\alias{probe.match.objfun,missing-method}
\alias{probe.match.objfun,pomp-method}
\alias{probe.match.objfun,probed_pomp-method}
\alias{probe.match.objfun-pomp}
\alias{probe.match.objfun-probed_pomp}
\alias{probed_pomp-methods}
\alias{probe_matched_pomp-methods}
\alias{spectd_pomp-methods}
\alias{spect_matched_pomp-methods}
\alias{summary,probed_pomp-method}
\alias{summary,probe_matched_pomp-method}
\alias{summary-probed_pomp}
\alias{summary-probe_matched_pomp}
\alias{probevals,probed_pomp-method}
\alias{probevals-probed_pomp}
\alias{covmat-probed_pomp}
\alias{covmat,probed_pomp-method}
\alias{plot,probed_pomp-method}
\alias{plot-probed_pomp}
\alias{summary,spectd_pomp-method}
\alias{summary,spect_matched_pomp-method}
\alias{summary-spectd_pomp}
\alias{summary-spect_matched_pomp}
\alias{plot,spectd_pomp-method}
\alias{plot,spect_matched_pomp-method}
\alias{plot,spectd_pomp-method}
\alias{plot-spectd_pomp}
\alias{as,probed_pomp-method}
\alias{logLik,probed_pomp-method}
\alias{logLik-probed_pomp}
\alias{coerce,probed_pomp,data.frame-method}
\description{
  \code{probe} applies one or more \dQuote{probes} to time series data and model simulations and compares the results.
  It can be used to diagnose goodness of fit and/or as the basis for \dQuote{probe-matching}, a generalized method-of-moments approach to parameter estimation.
  \code{probe.match} calls an optimizer to adjust model parameters to do probe-matching, i.e., to minimize the discrepancy between simulated and actual data.
  This discrepancy is measured using the \dQuote{synthetic likelihood} as defined by Wood (2010).
  \code{probe.match.objfun} constructs an objective function for probe-matching suitable for use in \code{optim}-like optimizers.
}
\usage{
\S4method{probe}{pomp}(object, probes, params, nsim, seed = NULL, \dots,
            verbose = getOption("verbose", FALSE))
\S4method{probe}{probed_pomp}(object, probes, nsim, seed = NULL, \dots,
            verbose = getOption("verbose", FALSE))
\S4method{probe.match.objfun}{pomp}(object, params, est, probes, nsim,
            seed = NULL, fail.value = NA, transform = FALSE, \dots)
\S4method{probe.match.objfun}{probed_pomp}(object, probes, nsim, seed, \dots)
\S4method{probe.match}{pomp}(object, start, est = character(0),
            probes, nsim, seed = NULL,
            method = c("subplex","Nelder-Mead","SANN","BFGS",
                       "sannbox","nloptr"),
            verbose = getOption("verbose"),
            fail.value = NA, transform = FALSE, \dots)
\S4method{probe.match}{probed_pomp}(object, probes, nsim, seed,
            \dots, verbose = getOption("verbose"))
\S4method{probe.match}{probe_matched_pomp}(object, est, probes,
            nsim, seed, transform, fail.value, \dots,
            verbose = getOption("verbose"))
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
    In the case of \code{probe} and \code{probe.match.objfun}, additional arguments are passed to \code{\link{pomp}}, allowing one to supply new or modify existing model characteristics or components.

    In the case of \code{probe.match}, these are passed to the optimizer algorithm.
    These are passed via the optimizer's \code{control} list (in the case of \code{optim}, \code{subplex}, and \code{sannbox}) or the \code{opts} list (in the case of \code{nloptr}).
  }
}
\details{
  A call to \code{probe} results in the evaluation of the probe(s) in \code{probes} on the data.
  Additionally, \code{nsim} simulated data sets are generated (via a call to \code{\link[=simulate-pomp]{simulate}}) and the probe(s) are applied to each of these.
  The results of the probe computations on real and simulated data are stored in an object of class \code{probed_pomp}.

  A call to \code{probe.match} results in an attempt to optimize the agreement between model and data, as measured by the specified probes, over the parameters named in \code{est}.
  The results, including coefficients of the fitted model and values of the probes for data and fitted-model simulations, are stored in an object of class \code{probe_matched_pomp}.

  The objective function minimized by \code{probe.match}---in a form suitable for use with \code{\link{optim}}-like optimizers---is created by a call to \code{probe.match.objfun}.
  Specifically, \code{probe.match.objfun} will return a function that takes a single numeric-vector argument that is assumed to cotain the parameters named in \code{est}, in that order.
  This function will return the negative synthetic log likelihood for the probes specified.
}
\value{
  \code{probe} returns an object of class \code{probed_pomp}, which is derived from the \code{pomp} class and contains additional information about the \code{probe} calculation.

  This information can be summarized via a call to \code{summary} and displayed graphically via a call to \code{plot}.

  \code{probe.match} returns an object of class \code{probe_matched_pomp}, which is derived from class \code{probed_pomp}.
  It therefore contains the information described above as well as information on the convergence of the optimization algorithm.
  This can be displayed via a call to \code{summary}.

  \code{probe.match.objfun} returns a function suitable for use as an objective function in an \code{\link{optim}}-like optimizer.
}
\section{Methods}{
  \describe{
    \item{plot}{
      displays diagnostic plots.
    }
    \item{summary}{
      displays summary information.
      The summary includes quantiles (fractions of simulations with probe values less than those realized on the data) and the corresponding two-sided p-values.
      In addition, the \dQuote{synthetic likelihood} (Wood 2010) is computed, under the assumption that the probe values are multivariate-normally distributed.
    }
    \item{probevals}{
      extracts the realized values of the probes on the data and on the simulations.
      These are returned in a list of two elements, \code{datvals} and \code{simvals}.
    }
    \item{logLik}{
      returns the synthetic likelihood for the probes.
      NB: in general, this is not the same as the likelihood.
    }
    \item{as.data.frame, as(object,"data.frame")}{
      coercing a \sQuote{probed_pomp} to a \sQuote{data.frame}, gives the realized values of the probes on the data and on the simulations.
      The variable \code{.id} indicates whether the probes are from the data or simulations.
    }
  }
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
