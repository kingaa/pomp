\name{Particle filter}
\title{Particle filter}
\alias{particle filter}
\alias{sequential Monte Carlo}
\alias{SMC}
\alias{pfilter}
\alias{pfilter,ANY-method}
\alias{pfilter,missing-method}
\alias{pfilter,pomp-method}
\alias{pfilter-pomp}
\alias{pfilter,pfilterd_pomp-method}
\alias{pfilter-pfilterd_pomp}
\alias{pfilterd_pomp-class}
\alias{pfilterd_pomp}
\alias{logLik,pfilterd_pomp-method}
\alias{logLik-pfilterd_pomp}
\alias{pred.mean}
\alias{pred.mean,pfilterd_pomp-method}
\alias{pred.mean-pfilterd_pomp}
\alias{pred.var}
\alias{pred.var,pfilterd_pomp-method}
\alias{pred.var-pfilterd_pomp}
\alias{filter.mean}
\alias{filter.mean,pfilterd_pomp-method}
\alias{filter.mean-pfilterd_pomp}
\alias{filter.traj}
\alias{filter.traj,pfilterd_pomp-method}
\alias{filter.traj-pfilterd_pomp}
\alias{eff.sample.size,pfilterd_pomp-method}
\alias{eff.sample.size-pfilterd_pomp}
\alias{cond.logLik,pfilterd_pomp-method}
\alias{cond.logLik-pfilterd_pomp}
\alias{as,pfilterd_pomp-method}
\alias{as.data.frame.pfilterd_pomp}
\alias{coerce,pfilterd_pomp,data.frame-method}
\description{
  A plain vanilla sequential Monte Carlo (particle filter) algorithm.
  Resampling is performed at each observation.
}
\usage{
\S4method{pfilter}{pomp}(object, params, Np, tol = 1e-17,
    max.fail = Inf, pred.mean = FALSE, pred.var = FALSE,
    filter.mean = FALSE, filter.traj = FALSE, save.states = FALSE,
    save.params = FALSE, verbose = getOption("verbose"), \dots)
\S4method{pfilter}{pfilterd_pomp}(object, params, Np, tol, \dots)
\S4method{logLik}{pfilterd_pomp}(object, \dots)
\S4method{cond.logLik}{pfilterd_pomp}(object, \dots)
\S4method{eff.sample.size}{pfilterd_pomp}(object, \dots)
\S4method{pred.mean}{pfilterd_pomp}(object, pars, \dots)
\S4method{pred.var}{pfilterd_pomp}(object, pars, \dots)
\S4method{filter.mean}{pfilterd_pomp}(object, pars, \dots)
\S4method{filter.traj}{pfilterd_pomp}(object, vars, \dots)
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
    This may be specified as a single positive integer, in which case the same number of particles will be used at each timestep.
    Alternatively, if one wishes the number of particles to vary across timesteps, one may specify \code{Np} either as a vector of positive integers of length \preformatted{length(time(object,t0=TRUE))} or as a function taking a positive integer argument.
    In the latter case, \code{Np(k)} must be a single positive integer, representing the number of particles to be used at the \code{k}-th timestep:
    \code{Np(0)} is the number of particles to use going from \code{timezero(object)} to \code{time(object)[1]},
    \code{Np(1)}, from \code{timezero(object)} to \code{time(object)[1]},
    and so on, while when \code{T=length(time(object,t0=TRUE))},
    \code{Np(T)} is the number of particles to sample at the end of the time-series.
    When \code{object} is of class \sQuote{mif}, this is by default the same number of particles used in the \code{mif} iterations.

    One should omit \code{Np} if \code{params} is a matrix of parameters, with one column for each particle.
    In this case, obviously, the number of particles is \code{ncol(params)}.
  }
  \item{tol}{
    positive numeric scalar; particles with likelihood less than \code{tol} are considered to be incompatible with the data.
    See the section on \emph{Filtering Failures} below for more information.
  }
  \item{max.fail}{
    integer; the maximum number of filtering failures allowed (see below).
    If the number of filtering failures exceeds this number, execution will terminate with an error.
    By default, \code{max.fail} is set to infinity, so no error can be triggered.
  }
  \item{pred.mean}{
    logical; if \code{TRUE}, the prediction means are calculated for the state variables and parameters.
  }
  \item{pred.var}{
    logical; if \code{TRUE}, the prediction variances are calculated for the state variables and parameters.
  }
  \item{filter.mean}{
    logical; if \code{TRUE}, the filtering means are calculated for the state variables and parameters.
  }
  \item{filter.traj}{
    logical; if \code{TRUE}, a filtered trajectory is returned for the state variables and parameters.
  }
  \item{save.states, save.params}{
    logical.
    If \code{save.states=TRUE}, the state-vector for each particle at each time is saved in the \code{saved.states} slot of the returned \sQuote{pfilterd_pomp} object.
    If \code{save.params=TRUE}, the parameter-vector for each particle at each time is saved in the \code{saved.params} slot of the returned \sQuote{pfilterd_pomp} object.
  }
  \item{verbose}{
    logical; if \code{TRUE}, progress information is reported as \code{pfilter} works.
  }
  \item{pars}{Names of parameters.}
  \item{vars}{Names of state variables.}
  \item{\dots}{
    Additional arguments are passed to \code{\link{pomp}}, allowing one to supply new or modify existing model characteristics or components.
  }
}
\value{
  An object of class \sQuote{pfilterd_pomp}.
  This class inherits from class \sQuote{pomp}.
  The following additional slots can be accessed via the \code{$} operator:
  \describe{
    \item{saved.states}{
      If \code{pfilter} was called with \code{save.states=TRUE}, this is the list of state-vectors at each time point, for each particle.
      It is a length-\code{ntimes} list of \code{nvars}-by-\code{Np} arrays.
      In particular, \code{saved.states[[t]][,i]} can be considered a sample from \eqn{f[X_t|y_{1:t}]}.
    }
    \item{saved.params}{
      If \code{pfilter} was called with \code{save.params=TRUE}, this is the list of parameter-vectors at each time point, for each particle.
      It is a length-\code{ntimes} list of \code{npars}-by-\code{Np} arrays.
      In particular, \code{saved.params[[t]][,i]} is the parameter portion of the i-th particle at time \eqn{t}.
    }
    \item{Np, tol, nfail}{
      the number of particles used, failure tolerance, and number of filtering failures (see below), respectively.
    }
  }
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
    \item{eff.sample.size}{
      Extracts the (time-dependent) estimated effective sample size, computed as
      \deqn{\left(\sum_i\!w_{it}^2\right)^{-1},}{1/(sum(w_it^2)),}
      where \eqn{w_{it}}{w_it} is the normalized weight of particle \eqn{i} at time \eqn{t}.
    }
    \item{pred.mean, pred.var}{
      Extract the mean and variance of the approximate prediction distribution.
      This prediction distribution is that of \deqn{X_t \vert y_1,\dots,y_{t-1},}{X_t | y_1,\dots,y_(t-1),} where \eqn{X_t}, \eqn{y_t} are the state vector and data, respectively, at time \eqn{t}.
    }
    \item{filter.mean}{
      Extract the mean of the filtering distribution, which is that of \deqn{X_t \vert y_1,\dots,y_t,}{X_t | y_1,\dots,y_t,} where \eqn{X_t}, \eqn{y_t} are the state vector and data, respectively, at time \eqn{t}.
    }
  }
}
\section{Filtering failures}{
  If the degree of disagreement between model and data becomes sufficiently large, a \dQuote{filtering failure} results.
  A filtering failure occurs when, at some time point, none of the \code{Np} particles is compatible with the data.
  In particular, if the conditional likelihood of a particle at any time is below the tolerance value \code{tol}, then that particle is considered to be uninformative and its likelihood is taken to be zero.
  A filtering failure occurs when this is the case for all particles.
  A warning is generated when this occurs unless the cumulative number of failures exceeds \code{max.fail}, in which case an error is generated.
}
\examples{
pompExample(gompertz)
pf <- pfilter(gompertz,Np=1000)	## use 1000 particles
plot(pf)
logLik(pf)
cond.logLik(pf)			## conditional log-likelihoods
eff.sample.size(pf)             ## effective sample size
logLik(pfilter(pf))      	## run it again with 1000 particles
## run it again with 2000 particles
pf <- pfilter(pf,Np=2000,filter.mean=TRUE)
fm <- filter.mean(pf)    	## extract the filtering means
}
\references{
  M. S. Arulampalam, S. Maskell, N. Gordon, & T. Clapp.
  A Tutorial on Particle Filters for Online Nonlinear, Non-Gaussian Bayesian Tracking.
  IEEE Trans. Sig. Proc. 50:174--188, 2002.
}
\author{Aaron A. King}
\seealso{
  \code{\link{pomp}}, \code{\link{pmcmc}}, \code{\link{bsmc2}},
  and the tutorials on the \href{https://kingaa.github.io/pomp}{package website}.
}
\keyword{ts}
