\name{probe}
\alias{probe}
\alias{probe,pomp-method}
\alias{probe,probed.pomp-method}
\alias{probe-pomp}
\alias{probe-probed.pomp}
\alias{probe.match}
\alias{probe.match,pomp-method}
\alias{probe.match,probed.pomp-method}
\alias{probe.match-pomp}
\alias{probe.match-probed.pomp}
\title{Probe a partially-observed Markov process.}
\description{
  \code{probe} applies one or more \dQuote{probes} to time series data and model simulations and compares the results.
  It can be used to diagnose goodness of fit and/or as the basis for \dQuote{probe-matching}, a generalized method-of-moments approach to parameter estimation.
  \code{probe.match} calls an optimizer to adjust model parameters to do probe-matching, i.e., to minimize the discrepancy between simulated and actual data.
}
\usage{
  probe(object, probes, \dots)
  \S4method{probe}{pomp}(object, probes, nsim, seed = NULL, \dots)
  \S4method{probe}{probed.pomp}(object, probes, nsim, seed = NULL, \dots)
  probe.match(object, start, est = character(0),
              probes, weights,
              nsim, seed = NULL,
              method = c("subplex","Nelder-Mead","SANN"),
              verbose = getOption("verbose"), 
              eval.only = FALSE, fail.value = NA, \dots)
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
    optional scalar; if non-\code{NA}, this value is substituted for non-finite values of the objective function.
  }
  \item{\dots}{
    Additional arguments.
    In the case of \code{spec}, these are currently ignored.
    In the case of \code{spec.match}, these are passed to \code{optim} or \code{subplex} in the \code{control} list.
  }
}
\details{
  A call to \code{probe} results in the evaluation of the probe(s) in \code{probes} on the data.
  Additionally, \code{nsim} simulated data sets are generated (via a call to \code{\link[=simulate-pomp]{simulate}}) and the probe(s) are applied to each of these.
  The results of the probe computations on real and simulated data are stored in an object of class \code{\link[=probed.pomp-class]{probed.pomp}}.

  A call to \code{probe.match} results in an attempt to optimize the agreement between model and data, as measured by the specified probes, over the parameters named in \code{est}.
  The results, including coefficients of the fitted model and values of the probes for data and fitted-model simulations, are stored in an object of class \code{\link[=probe.matched.pomp-class]{probe.matched.pomp}}.
}
\value{
  \code{probe} returns an object of class \code{\link[=probed.pomp-class]{probed.pomp}}.
  \code{probe.match} returns an object of class \code{\link[=probe.matched.pomp-class]{probe.matched.pomp}}.
}
\references{
  B. E. Kendall, C. J. Briggs, W. M. Murdoch, P. Turchin, S. P. Ellner, E. McCauley, R. M. Nisbet, S. N. Wood
  Why do populations cycle? A synthesis of statistical and mechanistic modeling approaches,
  Ecology, 80:1789--1805, 1999.
}
\author{
  Daniel C. Reuman (d.reuman at imperial dot ac dot uk)
  
  Aaron A. King (kingaa at umich dot edu)
}
\examples{
data(ou2)
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

ou2.bad <- ou2
coef(ou2.bad,c("x1.0","x2.0","alpha.1","alpha.4")) <- c(0,0,0.1,0.2)
bad <- probe(
             ou2.bad,
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




