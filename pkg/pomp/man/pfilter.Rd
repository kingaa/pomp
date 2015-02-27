\name{Particle filter}
\title{Particle filter}
\alias{particle filter}
\alias{sequential Monte Carlo}
\alias{SMC}
\alias{pfilter}
\alias{pfilter,pomp-method}
\alias{pfilter-pomp}
\alias{pfilter,pfilterd.pomp-method}
\alias{pfilter-pfilterd.pomp}
\alias{pfilterd.pomp-class}
\alias{logLik,pfilterd.pomp-method}
\alias{logLik-pfilterd.pomp}
\alias{$,pfilterd.pomp-method}
\alias{$-pfilterd.pomp}
\alias{pred.mean}
\alias{pred.mean,pfilterd.pomp-method}
\alias{pred.mean-pfilterd.pomp}
\alias{pred.var}
\alias{pred.var,pfilterd.pomp-method}
\alias{pred.var-pfilterd.pomp}
\alias{filter.mean}
\alias{filter.mean,pfilterd.pomp-method}
\alias{filter.mean-pfilterd.pomp}
\alias{eff.sample.size}
\alias{eff.sample.size,pfilterd.pomp-method}
\alias{eff.sample.size-pfilterd.pomp}
\alias{cond.logLik}
\alias{cond.logLik,pfilterd.pomp-method}
\alias{cond.logLik-pfilterd.pomp}
\alias{as,pfilterd.pomp-method}
\alias{as.data.frame.pfilterd.pomp}
\alias{coerce,pfilterd.pomp,data.frame-method}
\description{
  A plain vanilla sequential Monte Carlo (particle filter) algorithm.
  Resampling is performed at each observation.
}
\usage{
\S4method{pfilter}{pomp}(object, params, Np, tol = 1e-17,
    max.fail = Inf, pred.mean = FALSE, pred.var = FALSE,
    filter.mean = FALSE, save.states = FALSE,
    save.params = FALSE, seed = NULL,
    verbose = getOption("verbose"), \dots)
\S4method{pfilter}{pfilterd.pomp}(object, params, Np, tol, \dots)
\S4method{logLik}{pfilterd.pomp}(object, \dots)
\S4method{cond.logLik}{pfilterd.pomp}(object, \dots)
\S4method{eff.sample.size}{pfilterd.pomp}(object, \dots)
\S4method{pred.mean}{pfilterd.pomp}(object, pars, \dots)
\S4method{pred.var}{pfilterd.pomp}(object, pars, \dots)
\S4method{filter.mean}{pfilterd.pomp}(object, pars, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp} or inheriting class \code{pomp}.
  }
  \item{params}{
    A \code{npars} x \code{Np} numeric matrix containing the parameters corresponding to the initial state values in \code{xstart}.
    This must have a \sQuote{rownames} attribute.
    If it desired that all particles should share the same parameter values, one one may supply \code{params} as a named numeric vector.
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
    When \code{object} is of class \code{mif}, this is by default the same number of particles used in the \code{mif} iterations.
  }
  \item{tol}{
    positive numeric scalar; particles with likelihood less than \code{tol} are considered to be \dQuote{lost}.
    A filtering failure occurs when, at some time point, all particles are lost.
    When all particles are lost, the conditional likelihood at that time point is set to \code{tol}.
  }
  \item{max.fail}{
    integer; the maximum number of filtering failures allowed.
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
  \item{save.states, save.params}{
    logical.
    If \code{save.states=TRUE}, the state-vector for each particle at each time is saved in the \code{saved.states} slot of the returned \code{\link[=pfilterd.pomp-class]{pfilterd.pomp}} object.
    If \code{save.params=TRUE}, the parameter-vector for each particle at each time is saved in the \code{saved.params} slot of the returned \code{\link[=pfilterd.pomp-class]{pfilterd.pomp}} object.
  }
  \item{seed}{
    optional; an object specifying if and how the random number generator should be initialized (\sQuote{seeded}).
    If \code{seed} is an integer, it is passed to \code{set.seed} prior to any simulation and is returned as the \dQuote{seed} element of the return list.
    By default, the state of the random number generator is not changed and the value of \code{.Random.seed} on the call is stored in the \dQuote{seed} element of the return list.
  }
  \item{verbose}{
    logical; if \code{TRUE}, progress information is reported as \code{pfilter} works.
  }
  \item{pars}{Names of parameters.}
  \item{\dots}{
    additional arguments that override the defaults.
  }
}
\value{
  An object of class \code{pfilterd.pomp}.
  This class inherits from class \code{\link[=pomp-class]{pomp}} and contains the following additional slots:
  \describe{
    \item{pred.mean, pred.var, filter.mean}{
      matrices of prediction means, variances, and filter means, respectively.
      In each of these, the rows correspond to states and parameters (if appropriate), in that order, the columns to successive observations in the time series contained in \code{object}.
    }
    \item{eff.sample.size}{
      numeric vector containing the effective number of particles at each time point.
    }
    \item{cond.loglik}{
      numeric vector containing the conditional log likelihoods at each time point.
    }
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
    \item{seed}{
      the state of the random number generator at the time \code{pfilter} was called.
      If the argument \code{seed} was specified, this is a copy;
      if not, this is the internal state of the random number generator at the time of call.
    }
    \item{Np, tol, nfail}{
      the number of particles used, failure tolerance, and number of filtering failures, respectively.
    }
    \item{loglik}{
      the estimated log-likelihood.
    }
  }
  These can be accessed using the \code{$} operator as if the returned object were a list.
  Note that if the argument \code{params} is a named vector, then these parameters are included in the \code{params} slot of the returned \code{pfilterd.pomp} object.
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
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{
  \code{\link{pomp}}, \code{\link{mif}}, \code{\link{pmcmc}}, \code{\link{bsmc2}},
  and the tutorials on the \href{http://pomp.r-forge.r-project.org}{package website}.
}
\keyword{ts}
