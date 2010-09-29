\name{basic.probes}
\alias{basic.probes}
\alias{probe.mean}
\alias{probe.median}
\alias{probe.var}
\alias{probe.sd}
\alias{probe.period}
\alias{probe.quantile}
\alias{probe.cov}
\alias{probe.cor}
\alias{probe.acf}
\alias{probe.marginal}
\alias{probe.nlar}
\title{Some probes for partially-observed Markov processes}
\description{
  Several simple and configurable probes are provided in the package.
  These can be used directly and as examples for building custom probes.
}
\usage{
probe.mean(var, trim = 0, transform = identity, na.rm = TRUE)
probe.median(var, na.rm = TRUE)
probe.var(var, transform = identity, na.rm = TRUE)
probe.sd(var, transform = identity, na.rm = TRUE)
probe.marginal(var, ref, order = 3, diff = 1, transform = identity)
probe.nlar(var, lags, powers, transform = identity)
probe.acf(var, lag.max, type = c("covariance", "correlation"),
          transform = identity)
probe.cov(vars, lag, method = c("pearson", "kendall", "spearman"),
                 transform = identity)
probe.cor(vars, lag, method = c("pearson", "kendall", "spearman"),
                 transform = identity)
probe.period(var, kernel.width, transform = identity)
probe.quantile(var, prob, transform = identity)
}
\arguments{
  \item{var}{
    character; the name(s) of the observed variable(s).
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
  \item{lag.max}{
    The maximum lag at which the ACF is computed.
  }
  \item{lag}{
    In \code{probe.cov} and \code{probe.cor}, the lag between time series.
  }
  \item{type}{
    Compute autocorrelation or autocovariance?
  }
  \item{vars}{
    length-2 character vector;
    the names or indices of the observed variables between which covariance or correlation is to be computed.
  }
  \item{method}{
    The type of covariance or correlation to compute: see \code{\link{cov}} and \code{\link{cor}}.
  }
  \item{ref}{
    empirical reference distribution.
    Simulated data will be regressed against the values of \code{ref}, sorted and, optionally, differenced.
    The resulting regression coefficients capture information about the shape of the marginal distribution.
    A good choice for \code{ref} is the data itself.
  }
  \item{lags, powers}{
    the lags and powers, respectively, specifying the nonlinear autoregressive model that will be fit to the actual and simulated data.
    See Details, below, for a precise description.
  }
  \item{order}{
    order of polynomial regression.
  }
  \item{diff}{
    order of differencing to perform.
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
    \item{\code{probe.mean}, \code{probe.median}, \code{probe.var}, \code{probe.sd}}{
      return functions that compute the mean, median, variance, and standard deviation of variable \code{var}, respectively.
    }
    \item{\code{probe.period}}{
      returns a function that estimates the period of the Fourier component of the \code{var} series with largest power.
    }
    \item{\code{probe.marginal}}{
      returns a function that
      regresses the marginal distribution of variable \code{var} against the reference distribution \code{ref}.
      If \code{diff>0}, the data and the reference distribution are first differenced \code{diff} times and centered.
      Polynomial regression of order \code{order} is used.
      This probe returns \code{order} regression coefficients (the intercept is zero).
    }
    \item{\code{probe.nlar}}{
      returns a function that
      fit a nonlinear (polynomial) autoregressive model to the univariate series (variable \code{var}).
      Specifically, a model of the form \eqn{y_t = \sum \beta_k y_{t-\tau_k}^{p_k}+\epsilon_t}{y[t] = \sum beta[k] y[t-tau[k]]^p[k]+e[t]} will be fit, where \eqn{\tau_k}{tau[k]} are the \code{lags} and \eqn{p_k}{p[k]} are the \code{powers}.
      The data are first centered.
      This function returns the regression coefficients, \eqn{\beta_k}{beta[k]}.
    }
    \item{\code{probe.acf}}{
      returns a function that,
      if \code{type=="covariance"}, computes the autocovariance function of variable \code{var} at lags 0 through \code{lag.max};
      if \code{type=="correlation"}, computes the autocorrelation function of variable \code{var} at lags 1 through \code{lag.max}.
    }
    \item{\code{probe.cov}, \code{probe.cor}}{
      return functions that compute the covariance and correlation, respectively, between the two variables named in \code{vars} at relative lag \code{lag}.}
    \item{\code{probe.quantile}}{
      returns a function that estimates the \code{prob}-th quantile of variable \code{var}.
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
  Daniel C. Reuman (d.reuman at imperial dot ac dot uk)
  
  Aaron A. King (kingaa at umich dot edu)
}
\seealso{
  \link{pomp-class},
  \link{pomp-methods}
}
\keyword{models}
\keyword{ts}
