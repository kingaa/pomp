\name{probed.pomp-methods}
\docType{methods}
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
\title{Methods of the "probed.pomp", "probe.matched.pomp", "spect.pomp", and "spect.matched.pomp" classes}
\description{Methods of the \code{probed.pomp}, \code{probe.matched.pomp}, \code{spect.pomp}, and \code{spect.matched.pomp} classes}
\usage{
\S4method{summary}{probed.pomp}(object, \dots)
\S4method{plot}{probed.pomp}(x, y, \dots)
\S4method{summary}{probe.matched.pomp}(object, \dots)
\S4method{plot}{probe.matched.pomp}(x, y, \dots)
\S4method{summary}{spect.pomp}(object, \dots)
\S4method{logLik}{probed.pomp}(object, \dots)
\S4method{values}{probed.pomp}(object, \dots)
\S4method{plot}{spect.pomp}(x, y, max.plots.per.page = 4,
          plot.data = TRUE,
          quantiles = c(.025, .25, .5, .75, .975),
          quantile.styles = list(lwd=1, lty=1, col="gray70"),
          data.styles = list(lwd=2, lty=2, col="black"))
\S4method{summary}{spect.matched.pomp}(object, \dots)
\S4method{plot}{spect.matched.pomp}(x, y, \dots)
\S4method{as}{probed.pomp}(object, class)
}
\arguments{
  \item{object, x}{the object to be summarized or plotted.}
  \item{y}{ignored.}
  \item{max.plots.per.page}{maximum number of plots per page}
  \item{plot.data}{plot the data spectrum?}
  \item{quantiles}{quantiles to plot}
  \item{quantile.styles}{plot style parameters for the quantiles}
  \item{data.styles}{plot style parameters for the data spectrum}
  \item{class}{
    character;
    name of the class to which \code{object} should be coerced.
  }
  \item{\dots}{
    Further arguments (either ignored or passed to underlying functions).
  }
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
    \item{as}{
      when a \sQuote{probed.pomp} is coerced to a \sQuote{data.frame}, the first row gives the probes applied to the data; the rest of the rows give the probes evaluated on simulated data.
      The rownames of the result can be used to distinguish these.
    }
  }

  In addition, slots of this object can be accessed via the \code{$} operator.

}
\author{
  Daniel C. Reuman, Aaron A. King
}
\seealso{
  \code{\link{probe}},
  \code{\link[=probed.pomp-class]{probed.pomp}},
  \code{\link[=probe.matched.pomp-class]{probe.matched.pomp}},
  \code{\link{probe.match}}
}
\keyword{models}
\keyword{ts}
