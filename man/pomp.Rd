\name{pomp constructor}
\title{Constructor of the basic POMP object}
\alias{pomp constructor}
\alias{pomp}
\alias{pomp,data.frame-method}
\alias{data.frame-pomp}
\alias{pomp,numeric-method}
\alias{numeric-pomp}
\alias{pomp,matrix-method}
\alias{matrix-pomp}
\alias{pomp,pomp-method}
\alias{pomp-pomp}
\alias{pomp-class}
\alias{accumulator variables}
\description{
  This function constructs a \code{pomp} object,
  encoding a partially-observed Markov process model together with a uni- or multivariate time series.
  One implements the model by specifying its \emph{components}, each of which can be written as \R functions or, for much greater computational efficiency, using C code.
  The preferred way to specify most components (as detailed below) is through the use of \code{\link{Csnippet}}s, snippets of C that are compiled and linked into a running \R session.
}
\usage{
\S4method{pomp}{data.frame}(data, times, t0, \dots, rprocess, dprocess, rmeasure, dmeasure,
       measurement.model,
       skeleton, skeleton.type, skelmap.delta.t,
       initializer, rprior, dprior, params, covar, tcovar,
       obsnames, statenames, paramnames, covarnames, zeronames,
       PACKAGE, fromEstimationScale, toEstimationScale, globals, cdir, cfile)
\S4method{pomp}{numeric}(data, times, t0, \dots, rprocess, dprocess, rmeasure, dmeasure,
       measurement.model,
       skeleton, skeleton.type, skelmap.delta.t,
       initializer, rprior, dprior, params, covar, tcovar,
       obsnames, statenames, paramnames, covarnames, zeronames,
       PACKAGE, fromEstimationScale, toEstimationScale, globals, cdir, cfile)
\S4method{pomp}{matrix}(data, times, t0, \dots, rprocess, dprocess, rmeasure, dmeasure,
       measurement.model,
       skeleton, skeleton.type, skelmap.delta.t,
       initializer, rprior, dprior, params, covar, tcovar,
       obsnames, statenames, paramnames, covarnames, zeronames,
       PACKAGE, fromEstimationScale, toEstimationScale, globals, cdir, cfile)
\S4method{pomp}{pomp}(data, times, t0, \dots, rprocess, dprocess, rmeasure, dmeasure,
       measurement.model, skeleton, skeleton.type, skelmap.delta.t,
       initializer, rprior, dprior, params, covar, tcovar,
       obsnames, statenames, paramnames, covarnames, zeronames,
       PACKAGE, fromEstimationScale, toEstimationScale, globals, cdir, cfile)
}
\arguments{
  \item{data, times}{
    The time series data and times at which observations are made.
    \code{data} can be specified as a vector, a matrix, a data-frame.
    Alternatively, a \code{pomp} object can be supplied in the \code{data} argument.

    If \code{data} is a numeric vector, \code{times} must be a numeric vector of the same length.

    If \code{data} is a matrix, it should have dimensions \code{nobs} x \code{ntimes}, where \code{nobs} is the number of observed variables and \code{ntimes} is the number of times at which observations were made (i.e., each column is a distinct observation of the \code{nobs} variables).
    In this case, \code{times} must be given as a numeric vector (of length \code{ntimes}).

    If \code{data} is a data-frame, \code{times} must name the column of observation times.
    Note that, in this case, \code{data} will be internally coerced to an array with storage-mode \sQuote{double}.

    \code{times} must be numeric and strictly increasing.
  }
  \item{t0}{
    The zero-time, at which the stochastic dynamical system is to be initialized.
    This must be no later than the time of the first observation, i.e., \code{t0 <= times[1]}.
  }
  \item{rprocess}{
    optional function; a function of prototype \preformatted{
      rprocess(xstart,times,params,...)
    } that simulates from the unobserved process.
    The form of this function is given below.
    \pkg{pomp} provides a number of \code{\link{plugins}} that construct appropriate \code{rprocess} arguments corresponding to several common stochastic simulation algorithms.
    See below for more details.
  }
  \item{dprocess}{
    optional function; a function of prototype \preformatted{
      dprocess(x,times,params,log,...)
    } that evaluates the likelihood of a sequence of consecutive state transitions.
    The form of this function is given below.
    It is not typically necessary (or even feasible) to define \code{dprocess}.
    This functionality is provided to support future algorithm development.
  }
  \item{rmeasure}{
    optional; the measurement model simulator.
    This can be specified in one of four ways:
    \enumerate{
      \item as a function of prototype \preformatted{
	rmeasure(x,t,params,...)
      } that makes a draw from the observation process given states \code{x}, time \code{t}, and parameters \code{params}.
      \item as the name of a native (compiled) routine of type \code{pomp_measure_model_simulator} as defined in the header file \file{pomp.h}.
      (To view the header file, execute \preformatted{
	file.show(system.file("include/pomp.h",package="pomp"))
      } in an \R session.)
      \item using the formula-based \code{measurement.model} facility (see below).
      \item as a snippet of C code (via \code{\link{Csnippet}}) that draws from the observation process as above.
      The last is typically the preferred option, as it results in much faster code execution.
    }
  }
  \item{dmeasure}{
    optional; the measurement model probability density function.
    This can be specified in one of four ways:
    \enumerate{
      \item as a function of prototype \preformatted{
	dmeasure(y,x,t,params,log,...)
      } that computes the p.d.f. of \code{y} given \code{x}, \code{t}, and \code{params}.
      \item as the name of a native (compiled) routine of type \code{pomp_measure_model_density} as defined in the header file \file{pomp.h}.
      (To view the header file, execute \preformatted{
	file.show(system.file("include/pomp.h",package="pomp"))
      } in an \R session.)
      \item using the formula-based \code{measurement.model} facility (see below).
      \item as a snippet of C code (via \code{\link{Csnippet}}) that computes the p.d.f. as above.
    }
    The last is typically the preferred option, as it results in much faster code execution.
    As might be expected, if \code{log=TRUE}, this function should return the log likelihood.
  }
  \item{measurement.model}{
    optional; a formula or list of formulae, specifying the measurement model.
    These formulae are parsed internally to generate \code{rmeasure} and \code{dmeasure} functions.
    If \code{measurement.model} is given it overrides any specification of \code{rmeasure} or \code{dmeasure}.
    \strong{NB:} This is a convenience function, primarily designed to facilitate exploration;
    it will typically be possible to acclerate measurement model computations by writing \code{dmeasure} and/or \code{rmeasure} using \code{\link{Csnippet}}s.
  }
  \item{skeleton}{
    optional.
    The \code{skeleton} argument specifies the deterministic skeleton of the unobserved Markov process.

    If the model is a discrete-time Markov process, its deterministic skeleton is a map:
    indicate this by \code{skeleton=map(f, delta.t)}, where
    \code{f} is specifies the map and \code{delta.t} is the time interval corresponding to one application of the map.
    Unless otherwise specified, \code{delta.t = 1}.

    If the model is a continuous-time Markov process, its deterministic skeleton is a vectorfield:
    indicate this by \code{skeleton=vectorfield(f)}, where \code{f} specifies the vectorfield.

    In either case, \code{f} can be specified in one of three ways:
    \enumerate{
      \item as an \R function of prototype \preformatted{
	f(x,t,params,...)
      } that evaluates the deterministic skeleton at state \code{x} and time \code{t} given the parameters \code{params},
      \item as the name of a native (compiled) routine of type \code{pomp_skeleton} as defined in the header file \file{pomp.h}.
      (To view the header file, execute \preformatted{
	file.show(system.file("include/pomp.h",package="pomp"))
      }
      in an \R session.)
      \item as a snippet of C code (via \code{\link{Csnippet}}) that performs this evaluation.
      The latter is typically the preferred option, for reasons of computational efficiency.
    }
  }
  \item{skeleton.type, skelmap.delta.t}{
    Deprecated.
  }
  \item{initializer}{
    optional.
    The initializer gives the parameterization of the initial state of the unobserved Markov process.
    Specifically, given a vector of parameters, \code{params} and an initial time, \code{t0}, the initializer determines the state vector at time \code{t0}.

    By default, any parameters in \code{params}, the names of which end in \dQuote{\code{.0}}, are assumed to be initial values of states.
    To initialize the unobserved state process, these are simply copied over as initial conditions.
    The names of the resulting state variables are obtained by dropping the \dQuote{\code{.0}} suffix.

    A custom initializer can be supplied here in one of two formats:
    \enumerate{
      \item as an \R function of prototype \preformatted{
	initializer(params,t0,...)
      } that yields initial conditions for the state process when given a vector, \code{params}, of parameters.
      \item as a snippet of C code (via \code{\link{Csnippet}}).
      See the \code{\link{Csnippet}} help for rules on writing a custom initializer.
    }

    One can sample from the initial-state distribution using the \code{\link{init.state}} function.
  }
  \item{rprior}{
    optional; function drawing a sample from a prior distribution on parameters.
    This can be specified in one of three ways:
    \enumerate{
      \item as an \R function of prototype \preformatted{
	rprior(params,...)
      } that makes a draw from the prior distribution given \code{params},
      \item as the name of a native (compiled) routine of type \code{pomp_rprior} as defined in the header file \file{pomp.h}, or
      (To view the header file, execute \preformatted{
	file.show(system.file("include/pomp.h",package="pomp"))
      }
      in an \R session.)
      \item as a snippet of C code (via \code{\link{Csnippet}}).
    }
    As above, the latter is typically preferable.
  }
  \item{dprior}{
    optional; function evaluating the prior distribution.
    This can be specified in one of three ways:
    \enumerate{
      \item as an \R function of prototype \preformatted{
	dprior(params,log=FALSE,...)
      } that evaluates the prior probability density,
      \item as the name of a native (compiled) routine of type \code{pomp_dprior} as defined in the header file \file{pomp.h}, or
      (To view the header file, execute \preformatted{
	file.show(system.file("include/pomp.h",package="pomp"))
      } in an \R session.)
      \item as a snippet of C code (via \code{\link{Csnippet}}).
    }
    As above, the latter is typically preferable.
  }
  \item{params}{
    optional named numeric vector of parameters.
    This will be coerced internally to storage mode \code{double}.
  }
  \item{covar, tcovar}{
    optional matrix or data frame of covariates:
    \code{covar} is the table of covariates (one column per variable);
    \code{tcovar} the corresponding times (one entry per row of \code{covar}).

    \code{covar} can be specified as either a matrix or a data frame.
    In either case the columns are taken to be distinct covariates.
    If \code{covar} is a data frame, \code{tcovar} can be either the name or the index of the time variable.

    If a covariate table is supplied, then the value of each of the covariates is interpolated as needed.
    The resulting interpolated values are passed to the corresponding functions as a numeric vector named \code{covars};
    see below.
  }
  \item{obsnames, statenames, paramnames, covarnames}{
    optional character vectors specifying the names of observables, state variables, parameters, and covariates, respectively.
    These are only used in the event that one or more of the basic functions (\code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, \code{skeleton}, \code{rprior}, \code{dprior}) are defined using \code{\link{Csnippet}} or native routines.
  }
  \item{zeronames}{
    optional character vector specifying the names of accumulator variables (see below).
  }
  \item{PACKAGE}{
    optional string giving the name of the dynamically loaded library in which any native routines are to be found.
  }
  \item{fromEstimationScale, toEstimationScale}{
    optional parameter transformations.
    Many algorithms for parameter estimation search an unconstrained space of parameters.
    When working with such an algorithm and a model for which the parameters are constrained, it can be useful to transform parameters.
    \code{toEstimationScale} and \code{fromEstimationScale} are transformations from the model scale to the estimation scale, and vice versa, respectively.
    These functions must have arguments \code{params} and \code{\dots}.
    \code{toEstimationScale} should transform parameters from the scale that \code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, and \code{skeleton} use internally to the scale used in estimation.
    \code{fromEstimationScale} should be the inverse of \code{toEstimationScale}.
    The parameter transformations can be defined (as above) using either \R functions, native routines, or \code{\link{Csnippet}}s.

    Note that it is the user's responsibility to make sure that these transformations are mutually inverse.
    If \code{obj} is the constructed \code{pomp} object, and \code{coef(obj)} is non-empty, a simple check of this property is \preformatted{
      x <- coef(obj,transform=TRUE)
      obj1 <- obj
      coef(obj1,transform=TRUE) <- x
      identical(coef(obj),coef(obj1))
      identical(coef(obj1,transform=TRUE),x).
    } By default, both functions are the identity transformation.
    See the demos, \preformatted{
      demo(package="pomp"),
    } \code{\link{pompExample}}, and the tutorials on the \href{http://kingaa.github.io/pomp}{package website} for examples.
  }
  \item{globals}{
    optional character;
    C code that will be included in the source for (and therefore hard-coded into) the shared-object library created when the call to \code{pomp} uses \code{\link{Csnippet}}s.
    If no \code{Csnippet}s are used, \code{globals} has no effect.
  }
  \item{cdir, cfile}{
    optional character.
    \code{cdir} specifies the name of the directory within which \code{Csnippet} code will be compiled.
    By default, this is in a temporary directory specific to the running instance of \R.
    \code{cfile} gives the name of the file (in directory \code{cdir}) into which \code{Csnippet} codes will be written.
    By default, a random filename is used.
  }
  \item{\dots}{
    Any additional arguments given to \code{pomp} will be stored in the \code{\link[=pomp-class]{pomp}} object and passed as arguments to each of the basic functions whenever they are evaluated.
  }
}
\value{
  \code{pomp} returns an object of class \code{pomp}.
  If \code{data} is an object of class \code{pomp}, then by default the returned \code{pomp} object is identical to \code{data}.
  If additional arguments are given, these override the defaults.
}
\section{Important note}{
  \strong{
    It is not typically necessary (or even feasible) to define all of the components \code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, and \code{skeleton} in any given problem.
    Each algorithm makes use of only a subset of these components.
    Any algorithm requiring a component that has not been defined will return an informative error.
  }
}
\section{The state process model}{
  Specification of process-model codes \code{rprocess} and/or \code{dprocess} in most cases is facilitated by \pkg{pomp}'s so-called \code{\link{plugins}}, which have been developed to handle common use-cases.
  Currently, if one's process model evolves in discrete time or one is willing to make such an approximation (e.g., via an Euler approximation), then the \code{\link{euler.sim}}, \code{\link{discrete.time.sim}}, and \code{\link{onestep.sim}} plugins for \code{rprocess} and \code{\link{onestep.dens}} plugin for \code{dprocess} are available.
  In addition, for exact simulation of certain continuous-time Markov chains, an implementation of Gillespie's algorithm is available (see \code{\link{gillespie.sim}}).
  To learn more about the use of plugins, consult the help documentation (\code{\link{plugins}}) and the tutorials on the \href{http://kingaa.github.io/pomp}{package website}.
  Several of the demos and examples make use of these as well.

  In specific cases, it may be possible to obtain increased computational efficiency by writing custom versions of \code{rprocess} and/or \code{dprocess} instead of using the plugins.
  If such custom versions are desired, the following describes how these functions should be written.
  \describe{
    \item{\code{rprocess}}{
      If the plugins are not used \code{rprocess} must be an \R function with at least the following arguments:
      \code{xstart}, \code{times}, \code{params}, and \code{\dots}.
      It can also take additional arguments.
      It is guaranteed that these will be filled with the corresponding elements the user has included as additional arguments in the construction of the \code{pomp} object.
      
      In calls to \code{rprocess}, \code{xstart} can be assumed to be an \code{nvar} x \code{nrep} matrix;
      its rows correspond to components of the state vector and columns correspond to independent realizations of the process.
      \code{params} will similarly be an \code{npar} x \code{nrep} matrix with rows corresponding to parameters and columns corresponding to independent realizations.
      Note that the columns of \code{params} correspond to those of \code{xstart};
      in particular, they will agree in number.
      Both \code{xstart} and \code{params} are guaranteed to have rownames.
      
      \code{rprocess} must return a rank-3 array with rownames.
      Suppose \code{x} is the array returned.
      Then \code{dim(x)=c(nvars,nrep,ntimes)}, where \code{ntimes} is the length of the vector \code{times}.
      \code{x[,j,k]} is the value of the state process in the \code{j}-th realization at time \code{times[k]}.
      In particular, \code{x[,,1]} must be identical to \code{xstart}.
      The rownames of \code{x} must correspond to those of \code{xstart}.
    }
    \item{\code{dprocess}}{
      If the plugins are not used, \code{dprocess} must have at least the following arguments:
      \code{x}, \code{times}, \code{params}, \code{log}, and \code{\dots}.
      It may take additional arguments: again, these will be filled with the corresponding elements the user defines when the \code{pomp} object is constructed.
      
      In calls to \code{dprocess}, \code{x} may be assumed to be an \code{nvars} x \code{nrep} x \code{ntimes} array, where these terms have the same meanings as above.
      \code{params} will be a matrix with rows corresponding to individual parameters and columns corresponding to independent realizations.
      The columns of \code{params} correspond to those of \code{x}; in particular, they will agree in number.
      Both \code{x} and \code{params} are guaranteed to have rownames.
      
      \code{dprocess} must return a \code{nrep} x \code{ntimes-1} matrix.
      Suppose \code{d} is the array returned.
      \code{d[j,k]} is the probability density of the transition from state \code{x[,j,k-1]} at time \code{times[k-1]} to state \code{x[,j,k]} at time \code{times[k]}.
      If \code{log=TRUE}, then the log of the pdf must be returned.

      \strong{In writing this function, you may assume that the transitions are consecutive.}
      It should be clear that, but for this assumption, it will in general be impossible to write the transition probabilities explicitly.
      In such cases, algorithms that make no use of \code{dprocess}, which are said to have the \dQuote{plug and play} property, are useful.
      Most of the algorithms in \pkg{pomp} have this property.
      In particular, \strong{at present, no methods in \pkg{pomp} make use of \code{dprocess}.}
    }
  }
}
\section{The observation process model}{
  The following is a guide to writing the measurement model components as \R functions.
  For a description on how to write these components using \code{\link{Csnippet}}s, see the tutorials on the \href{http://kingaa.github.io/pomp}{package website}.
  \describe{
    \item{\code{rmeasure}}{
      if provided, must take at least the arguments \code{x}, \code{t}, \code{params}, and \code{\dots}.
      It may take additional arguments, which will be filled with user-specified data as above.
      \code{x} will be a named numeric vector of length \code{nvars} (which has the same meaning as above).
      \code{t} will be a scalar quantity, the time at which the measurement is made.
      \code{params} will be a named numeric vector of length \code{npars}.
      The \code{rmeasure} function may take additional arguments which will be filled with user-specified data as above.
      
      \code{rmeasure} must return a named numeric vector of length \code{nobs}, the number of observable variables.
    }
    \item{\code{dmeasure}}{
      if provided, must take at least the arguments \code{y}, \code{x}, \code{t}, \code{params}, \code{log}, and \code{\dots}.
      \code{y} will be a named numeric vector of length \code{nobs} containing (actual or simulated) values of the observed variables;
      \code{x} will be a named numeric vector of length \code{nvar} containing state variables;
      \code{params} will be a named numeric vector containing parameters;
      and \code{t} will be a scalar, the corresponding observation time.
      The \code{dmeasure} function may take additional arguments which will be filled with user-specified data as above.

      \code{dmeasure} must return a single numeric value, the probability density of \code{y} given \code{x} at time \code{t}.
      If \code{log=TRUE}, then \code{dmeasure} should return the log of the probability density.
    }
  }
}
\section{The deterministic skeleton}{
  The following describes how to specify the deterministic skeleton as an \R function.
  For a description on how to write this component using \code{\link{Csnippet}}s, see the tutorials on the \href{http://kingaa.github.io/pomp}{package website} and the \code{\link{Csnippet}} help.

  If \code{skeleton} if provided, must have at least the arguments \code{x}, \code{t}, \code{params}, and \code{\dots}.
  \code{x} is a numeric vector containing the coordinates of a point in state space at which evaluation of the skeleton is desired.
  \code{t} is a numeric value giving the time at which evaluation of the skeleton is desired.
  Of course, these will be irrelevant in the case of an autonomous skeleton.
  \code{params} is a numeric vector holding the parameters.
  \code{skeleton} may take additional arguments, which will be filled, as above, with user-specified data.

  \code{skeleton} must return a numeric vector of the same length as \code{x}, which contains the value vectorfield (if the dynamical system is continuous) or the value of the map (if the dynamical system is discrete), at the point \code{x} at time \code{t}.
}
\section{The state-process initializer}{
  if provided, must have at least the arguments \code{params}, \code{t0}, and \code{\dots}.
  \code{params} will be a named numeric vector of parameters.
  \code{t0} will be the time at which initial conditions are desired.
  \code{initializer} must return a named numeric vector of initial states.
}
\section{Covariates}{
  If the \code{pomp} object contains covariates (via the \code{covar} argument; see above), then whenever any of the \R functions described above are called, they will each be supplied with an additional argument \code{covars}.
  This will be a named numeric vector containing the (interpolated) values of the covariates at the time \code{t}.
  In particular, \code{covars} will have one value for each column of the covariate table.
}
\section{Accumulator variables}{
  In formulating models, one often wishes to define a state variable that will accumulate some quantity over the interval between successive observations.
  \pkg{pomp} provides a facility to make such features more convenient.
  Specifically, variables named in the \code{pomp}'s \code{zeronames} argument will be set to zero immediately following each observation.
  See \code{euler.sir} and the tutorials on the \href{http://kingaa.github.io/pomp}{package website} for examples.
}
\examples{
\dontrun{
pompExample()
pomp.home <- system.file("examples",package="pomp")
pomp.examples <- list.files(pomp.home)
file.show(
          file.path(pomp.home,pomp.examples),
          header=paste("======",pomp.examples,"=======")
         )
}
}
\section{Warning}{
  Some error checking is done by \code{pomp}, but complete error checking is impossible.  
  If the user-specified functions do not conform to the above specifications, then the results may be invalid.
  In particular, if both \code{rmeasure} and \code{dmeasure} are specified, the user should verify that these two functions correspond to the same probability distribution.
  If \code{skeleton} is specified, the user is responsible for verifying that it corresponds to a deterministic skeleton of the model.
  Each \pkg{pomp}-package algorithm uses some subset of the five basic functions (\code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, \code{skeleton}).
  If an algorithm requires a component that has not been specified, an informative error will be generated.
}
\author{Aaron A. King}
\seealso{
  \link{pomp methods},
  \link{pomp low-level interface},
  \link{process model plugins}
 }
\keyword{models}
\keyword{interface}
\keyword{programming}
\keyword{ts}
