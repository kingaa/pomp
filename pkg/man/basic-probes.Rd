\name{basic.probes}
\alias{basic.probes}
\alias{probe.mean}
\alias{probe.var}
\alias{probe.sd}
\alias{probe.period}
\alias{probe.quantile}
\alias{probe.acf}
\alias{probe.cov}
\alias{probe.cor}
\title{Some probes for partially-observed Markov processes}
\description{
  Several simple and configurable probes are provided in the package.
  These can be used directly and as examples for building custom probes.
}
\usage{
probe.mean(var, trim = 0, transform = identity, na.rm = TRUE)
probe.var(var, transform = identity, na.rm = TRUE)
probe.sd(var, transform = identity, na.rm = TRUE)
probe.period(var, kernel.width, transform = identity)
probe.quantile(var, prob, transform = identity)
probe.acf(var, lag = 0, type = c("correlation", "covariance", "partial"),
          transform = identity, \dots)
probe.cov(vars, lag, method = c("pearson", "kendall", "spearman"),
                 transform = identity)
probe.cor(vars, lag, method = c("pearson", "kendall", "spearman"),
                 transform = identity)
}
\arguments{
  \item{var}{
    character; the name or index of the observed variable.
  }
  \item{trim}{
    the fraction of observations to be trimmed (see \code{\link{mean}}).
  }
  \item{transform}{
    transformation to be applied to the data before the probe is computed.
  }
  \item{na.rm}{
    if \code{TRUE}, remove all NA observations prior to computing the probe.
  }
  \item{kernel.width}{
    width of modified Daniell smoothing kernel to be used in power-spectrum computation:
    see \code{\link{kernel}}.
  }
  \item{prob}{
    a single probability; the quantile to compute: see \code{\link{quantile}}.
  }
  \item{lag}{
    In \code{probe.acf}, the lag at which the ACF is computed.
    In \code{probe.cov} and \code{probe.cor}, the lag between time series.
  }
  \item{type}{
    The type of ACF to compute: see \code{\link{acf}}.
  }
  \item{vars}{
    length-2 character vector; the names or indices of the observed variables between which covariance or correlation is to be computed.
  }
  \item{method}{
    The type of covariance or correlation to compute: see \code{\link{cov}} and \code{\link{cor}}.
  }
  \item{\dots}{
    Additional arguments to be passed through to the probe computation.
  }
}
\value{
  A call to any one of these functions returns a probe function, suitable for use in \code{\link{probe}} or \code{\link{probe.match}}.
  That is, the function returned by each of these takes a data array (such as comes from a call to \code{\link{data.array}}) as input and returns a single numerical value.
}
\details{
  Each of these functions is relatively simple.
  See the source code for a complete understanding of what each does.
  \describe{
    \item{\code{probe.mean}, \code{probe.var}, \code{probe.sd}}{
      compute the mean, variance, and standard deviation of variable \code{var}, respectively.
    }
    \item{\code{probe.period}}{
      estimates the period of the Fourier component of the \code{var} series with largest power.
    }
    \item{\code{probe.quantile}}{
      estimates the \code{prob}-th quantile of variable \code{var}.
    }
    \item{\code{probe.acf}}{
      computes the autocorrelation or autocovariance function of variable \code{var} at lag \code{lag}.
    }
    \item{\code{probe.cov}, \code{probe.cor}}{
      compute the covariance and correlation, respectively, between the two variables named in \code{vars} at relative lag \code{lag}.}
  }
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
\seealso{
  \link{pomp-class},
  \link{pomp-methods}
}
\keyword{models}
\keyword{ts}
