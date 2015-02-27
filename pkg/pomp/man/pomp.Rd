\name{pomp}
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
\alias{pompLoad}
\alias{pompLoad-pomp}
\alias{pompLoad,pomp-method}
\alias{pompUnload}
\alias{pompUnload-pomp}
\alias{pompUnload,pomp-method}
\title{Partially-observed Markov process object.}
\description{
  This function creates a \code{pomp} object,
  encoding a partially-observed Markov process model together with a uni- or multi-variate time series.
}
\usage{
\S4method{pomp}{data.frame}(data, times, t0, \dots, rprocess, dprocess, rmeasure, dmeasure,
       measurement.model,
       skeleton, skeleton.type = c("map","vectorfield"), skelmap.delta.t = 1,
       initializer, rprior, dprior, params, covar, tcovar,
       obsnames, statenames, paramnames, covarnames, zeronames,
       PACKAGE, parameter.transform, parameter.inv.transform, globals)
\S4method{pomp}{numeric}(data, times, t0, \dots, rprocess, dprocess, rmeasure, dmeasure,
       measurement.model,
       skeleton, skeleton.type = c("map","vectorfield"), skelmap.delta.t = 1,
       initializer, rprior, dprior, params, covar, tcovar,
       obsnames, statenames, paramnames, covarnames, zeronames,
       PACKAGE, parameter.transform, parameter.inv.transform, globals)
\S4method{pomp}{matrix}(data, times, t0, \dots, rprocess, dprocess, rmeasure, dmeasure,
       measurement.model,
       skeleton, skeleton.type = c("map","vectorfield"), skelmap.delta.t = 1,
       initializer, rprior, dprior, params, covar, tcovar,
       obsnames, statenames, paramnames, covarnames, zeronames,
       PACKAGE, parameter.transform, parameter.inv.transform, globals)
\S4method{pomp}{pomp}(data, times, t0, \dots, rprocess, dprocess, rmeasure, dmeasure,
       measurement.model, skeleton, skeleton.type, skelmap.delta.t,
       initializer, rprior, dprior, params, covar, tcovar,
       obsnames, statenames, paramnames, covarnames, zeronames,
       PACKAGE, parameter.transform, parameter.inv.transform, globals)
\S4method{pompLoad}{pomp}(object, \dots)
\S4method{pompUnload}{pomp}(object, \dots)
}
\arguments{
  \item{data, times}{
    The time series data and times at which observations are made.
    \code{data} can be specified as a vector, a matrix, a data-frame, or a \code{pomp} object.

    If \code{data} is a numeric vector, \code{times} must be a numeric vector of the same length.

    If \code{data} is a matrix, it should have dimensions \code{nobs} x \code{ntimes}, where \code{nobs} is the number of observed variables and \code{ntimes} is the number of times at which observations were made (i.e., each column is a distinct observation of the \code{nobs} variables).
    In this case, \code{times} must be given as a numeric vector (of length \code{ntimes}).

    If \code{data} is a data-frame, \code{times} must name the column of observation times.
    Note that, in this case, \code{data} will be internally coerced to an array with storage-mode \sQuote{double}.
    Note that the \code{times} must be numeric and strictly increasing.
  }
  \item{t0}{
    The zero-time: the stochastic dynamical system is initialized at time \code{t0}.
    This must be no later than the time of the first observation, i.e., \code{t0 <= times[1]}.
  }
  \item{rprocess}{
    optional function; a function of prototype \preformatted{rprocess(xstart,times,params,\dots)} that simulates from the unobserved process.
    The easiest way to specify \code{rprocess} is to use one of the \code{\link{plugins}} provided as part of the \pkg{pomp} package.
    See below for details.
  }
  \item{dprocess}{
    optional function; a function of prototype \preformatted{dprocess(x,times,params,log,\dots)} that evaluates the likelihood of a sequence of consecutive state transitions.
    The easiest way to specify \code{dprocess} is to use one of the \code{\link{plugins}} provided as part of the \pkg{pomp} package.
    It is not typically necessary (or even feasible) to define \code{dprocess}.
    See below for details.
  }
  \item{rmeasure}{
    optional; the measurement model simulator.
    This can be specified in one of four ways:
    (1) as a function of prototype \preformatted{rmeasure(x,t,params,\dots)} that makes a draw from the observation process given states \code{x}, time \code{t}, and parameters \code{params}.
    (2) as a snippet of C code (via \code{\link{Csnippet}}) that draws from the observation process as above.
    (3) as the name of a native (compiled) routine with prototype \dQuote{pomp_measure_model_simulator} as defined in the header file \dQuote{pomp.h}.
    In the above cases, if the measurement model depends on covariates, the optional argument \code{covars} will be filled with interpolated values at each call.
    (4) using the formula-based \code{measurement.model} facility (see below).
  }
  \item{dmeasure}{
    optional; the measurement model probability density function.
    This can be specified in one of four ways:
    (1) as a function of prototype \preformatted{dmeasure(y,x,t,params,log,\dots)} that computes the p.d.f. of \code{y} given \code{x}, \code{t}, and \code{params}.
    (2) as a snippet of C code (via \code{\link{Csnippet}}) that computes the p.d.f. as above.
    (3) as the name of a native (compiled) routine with prototype \dQuote{pomp_measure_model_density} as defined in the header file \dQuote{pomp.h}.
    In the above cases, if the measurement model depends on covariates, the optional argument \code{covars} will be filled with interpolated values at each call.
    (4) using the formula-based \code{measurement.model} facility (see below).
    As might be expected, if \code{log=TRUE}, this function should return the log likelihood.
  }
  \item{measurement.model}{
    optional; a formula or list of formulae, specifying the measurement model.
    These formulae are parsed internally to generate \code{rmeasure} and \code{dmeasure} functions.
    If \code{measurement.model} is given it overrides any specification of \code{rmeasure} or \code{dmeasure}.
    \strong{NB:} This is a convenience function, primarily designed to facilitate exploration;
    it will typically be possible to acclerate measurement model computations by writing \code{dmeasure} and/or \code{rmeasure} functions directly.
  }
  \item{skeleton, skeleton.type, skelmap.delta.t}{
    The function \code{skeleton} specifies the deterministic skeleton of the unobserved Markov process.
    If we are dealing with a discrete-time Markov process, its deterministic skeleton is a map:
    indicate this by specifying \code{skeleton.type="map"}.
    If we are dealing with a continuous-time Markov process, its deterministic skeleton is a vectorfield:
    indicate this by specifying \code{skeleton.type="vectorfield"}.

    The skeleton function can be specified in one of three ways:
    (1) as an \R function of prototype \preformatted{skeleton(x,t,params,\dots)} that evaluates the deterministic skeleton at state \code{x} and time \code{t} given the parameters \code{params},
    (2) as a snippet of C code (via \code{\link{Csnippet}}) that performs this evaluation, or
    (3) as the name of a native (compiled) routine with prototype \dQuote{pomp_skeleton} as defined in the header file \dQuote{pomp.h}.
    If the deterministic skeleton depends on covariates, the optional argument \code{covars} will be filled with interpolated values of the covariates at the time \code{t}.

    With a discrete-time skeleton, the default assumption is that time advances 1 unit per iteration of the map;
    to change this, set \code{skelmap.delta.t} to the appropriate time-step.
  }
  \item{initializer}{
    optional function of prototype \preformatted{initializer(params,t0,\dots)} that yields initial conditions for the state process when given a vector, \code{params}, of parameters.

    By default (i.e., if it is unspecified when \code{pomp} is called), the initializer assumes any parameters in \code{params} the names of which end in \dQuote{\code{.0}} are initial values.
    These are simply copied over as initial conditions when \code{init.state} is called (see \code{\link{init.state-pomp}}).
    The names of the state variables are the same as the corresponding initial value parameters, but with the \dQuote{\code{.0}} dropped.
  }
  \item{rprior}{
    optional; function drawing a sample from a prior distribution on parameters.
    This can be specified in one of three ways:
    (1) as an \R function of prototype \preformatted{rprior(params,\dots)} that makes a draw from the prior distribution given \code{params},
    (2) as a snippet of C code (via \code{\link{Csnippet}}) that draws from the prior, or
    (3) as the name of a native (compiled) routine with prototype \dQuote{pomp_rprior} as defined in the header file \dQuote{pomp.h}.
  }
  \item{dprior}{
    optional; function evaluating the prior distribution.
    This can be specified in one of three ways:
    (1) as an \R function of prototype \preformatted{dprior(params,log=FALSE,\dots)} that evaluates the prior probability density,
    (2) as a snippet of C code (via \code{\link{Csnippet}}) that computes this p.d.f.,
    (3) as the name of a native (compiled) routine with prototype \dQuote{pomp_dprior} as defined in the header file \dQuote{pomp.h}.
  }
  \item{params}{
    optional named numeric vector of parameters.
    This will be coerced internally to storage mode \code{double}.
  }
  \item{covar, tcovar}{
    An optional table of covariates:
    \code{covar} is the table (one column per variable);
    \code{tcovar} the corresponding times (one entry per row of \code{covar}).

    \code{covar} can be specified as either a matrix or a data frame.
    In either case the columns are taken to be distinct covariates.
    If \code{covar} is a data frame, \code{tcovar} can be either the name or the index of the time variable.

    If a covariate table is supplied, then the value of each of the covariates is interpolated as needed, i.e., whenever \code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, \code{skeleton}, \code{rprior}, \code{dprior}, or \code{init.state} is evaluated.
    The resulting interpolated values are passed to the corresponding functions as a numeric vector named \code{covars}.
  }
  \item{obsnames, statenames, paramnames, covarnames}{
    Optional character vectors specifying the names of observables, state variables, parameters, and covariates, respectively.
    These are only used in the event that one or more of the basic functions (\code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, \code{skeleton}, \code{rprior}, \code{dprior}) are defined using \code{\link{Csnippet}} or native routines.
    In these cases, these name vectors are matched against the corresponding names and the indices of the names are passed to the native routines.
    These vectors should list the parameters in the order that the native code expects them.
    Using this facility allows one to write one or more of the basic functions in native code in a way that does not depend on the order of states, parameters, and covariates at run time.
    See the \dQuote{Advanced topics in pomp} vignette for more on this topic, including examples.
  }
  \item{zeronames}{
    Optional character vector specifying the names of accumulator variables.
    See the \dQuote{Advanced topics in pomp} vignette for a discussion of accumulator variables.
  }
  \item{PACKAGE}{
    An optional string giving the name of the dynamically loaded library in which any native routines are to be found.
  }
  \item{parameter.transform, parameter.inv.transform}{
    Optional functions specifying parameter transformations.
    These functions must have arguments \code{params} and \code{\dots}.
    \code{parameter.transform} should transform parameters from the user's scale to the scale that \code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, \code{skeleton}, and \code{initializer} will use internally.
    \code{parameter.inv.transform} should be the inverse of \code{parameter.transform}.
    The parameter transformations can be defined (as in the case of \code{rmeasure}, \code{dmeasure}, \code{skeleton}, \code{rprior}, \code{dprior} above) using \R functions, \code{\link{Csnippet}}s, or native routines.
    Note that it is the user's responsibility to make sure this holds.
    If \code{obj} is the constructed \code{pomp} object, and \code{coef(obj)} is non-empty, a simple check of this is

    \preformatted{
      x <- coef(obj,transform=TRUE)
      obj1 <- obj
      coef(obj1,transform=TRUE) <- x
      identical(coef(obj),coef(obj1))
      identical(coef(obj1,transform=TRUE),x).
    }
    
    By default, both functions are the identity transformation.
    See the \dQuote{introduction_to_pomp} vignette for an example.
  }
  \item{globals}{
    optional character;
    C code that will be included in the source for (and therefore hard-coded into) the shared-object library created when the call to \code{pomp} uses \code{Csnippets}.
    If no \code{Csnippet}s are used, \code{globals} has no effect.
  }
  \item{object}{
    a \code{pomp} object.
  }
  \item{\dots}{
    Any additional arguments given to \code{pomp} will be stored in the \code{\link[=pomp-class]{pomp}} object and passed as arguments to each of the basic functions whenever they are evaluated.
  }
}
\value{
  An object of class \code{pomp}.
  If \code{data} is an object of class \code{\link[=pomp-class]{pomp}}, then by default the returned \code{pomp} object is identical to \code{data}.
  If additional arguments are given, these override the defaults.
}
\details{
  \strong{
    It is not typically necessary (or even feasible) to define all of the functions \code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, and \code{skeleton} in any given problem.
    Each algorithm makes use of a different subset of these functions.
  }
  Specification of process-model codes \code{rprocess} and/or \code{dprocess} in most cases is facilitated by so-called \code{\link{plugins}} that have been developed to handle common use-cases.
  Currently, if one's process model evolves in discrete time or one is willing to make such an approximation (e.g., via an Euler approximation), then the \code{\link{euler.sim}}, \code{\link{discrete.time.sim}}, or \code{\link{onestep.sim}} plugin for \code{rprocess} and \code{\link{onestep.dens}} plugin for \code{dprocess} are available.
  For exact simulation of certain continuous-time Markov chains, an implementation of Gillespie's algorithm is available (see \code{\link{gillespie.sim}}).
  To use the plugins, consult the help documentation (\code{?\link{plugins}}) and the vignettes.

  In specific cases, it may be possible to obtain increased computational efficiency by writing custom versions of \code{rprocess} and/or \code{dprocess}.
  See the \dQuote{Advanced topics in pomp} vignette for a discussion of this.
  If such custom versions are desired, the following describes how each of these functions should be written in this case.
  \describe{
    \item{\code{rprocess}}{
      If the plugins are not used \code{rprocess} must have at least the following arguments:
      \code{xstart}, \code{times}, \code{params}, and \code{\dots}.
      It can also take additional arguments.
      It is guaranteed that these will be filled with the corresponding elements the user has included as additional arguments in the construction of the \code{pomp} object.
      
      In calls to \code{rprocess}, \code{xstart} can be assumed to be a rank-2 array (matrix) with rows corresponding to state variables and columns corresponding to independent realizations of the process.
      \code{params} will similarly be a rank-2 array with rows corresponding to parameters and columns corresponding to independent realizations.
      The columns of \code{params} correspond to those of \code{xstart};
      in particular, they will agree in number.
      Both \code{xstart} and \code{params} will have rownames, which are available for use by the user.
      
      \code{rprocess} must return a rank-3 array with rownames.
      Suppose \code{x} is the array returned.
      Then \code{dim(x)=c(nvars,nreps,ntimes)}, where \code{nvars} (=\code{nrow(xstart)}) is the number of state variables, \code{nreps} (=\code{ncol(xstart)}) is the number of independent realizations simulated, and \code{ntimes} is the length of the vector \code{times}.
      \code{x[,j,k]} is the value of the state process in the \code{j}-th realization at time \code{times[k]}.
      In particular, \code{x[,,1]} must be identical to \code{xstart}.
      The rownames of \code{x} must correspond to those of \code{xstart}.
    }
    \item{\code{dprocess}}{
      If the plugins are not used, \code{dprocess} must have at least the following arguments:
      \code{x}, \code{times}, \code{params}, \code{log}, and \code{\dots}.
      It may take additional arguments.
      It is guaranteed that these will be filled with the corresponding elements the user has included as additional arguments in the construction of the \code{pomp} object.
      
      In calls to \code{dprocess}, \code{x} may be assumed to be an \code{nvars} x \code{nreps} x \code{ntimes} array, where these terms have the same meanings as above.
      \code{params} will be a rank-2 array with rows corresponding to individual parameters and columns corresponding to independent realizations.
      The columns of \code{params} correspond to those of \code{x}; in particular, they will agree in number.
      Both \code{x} and \code{params} will have rownames, available for use by the user.
      
      \code{dprocess} must return a rank-2 array (matrix).
      Suppose \code{d} is the array returned.
      Then \code{dim(d)=c(nreps,ntimes-1)}.
      \code{d[j,k]} is the probability density of the transition from state \code{x[,j,k-1]} at time \code{times[k-1]} to state \code{x[,j,k]} at time \code{times[k]}.
      If \code{log=TRUE}, then the log of the pdf is returned.

      \strong{In writing this function, you may assume that the transitions are consecutive.}
      It should be quite clear that, but for this assumption, it would be quite difficult in general to write the transition probabilities.
      In fact, from one perspective, the algorithms in \pkg{pomp} are designed to overcome just this difficulty.

      \strong{At present, no methods in \pkg{pomp} make use of \code{dprocess}.}
    }
  }

  The following is a guide to writing the measurement-model, deterministic skeleton, and initializer components. 
  \describe{
    \item{\code{rmeasure}}{
      if provided, must take at least the arguments \code{x}, \code{t}, \code{params}, and \code{\dots}.
      It may take additional arguments, which will be filled with user-specified data as above.
      \code{x} may be assumed to be a named numeric vector of length \code{nvars}, (which has the same meanings as above).
      \code{t} is a scalar quantity, the time at which the measurement is made.
      \code{params} may be assumed to be a named numeric vector of length \code{npars}.
      
      \code{rmeasure} must return a named numeric vector.
      If \code{y} is the returned vector, then \code{length(y)=nobs}, where \code{nobs} is the number of observable variables.
    }
    \item{\code{dmeasure}}{
      if provided, must take at least the arguments \code{y}, \code{x}, \code{t}, \code{params}, \code{log}, and \code{\dots}.
      \code{y} may be assumed to be a named numeric vector of length \code{nobs} containing (actual or simulated) values of the observed variables;
      \code{x} will be a named numeric vector of length \code{nvar} containing state variables
      \code{params}, a named numeric vector containing parameters;
      and \code{t}, a scalar, the corresponding observation time.
      It may take additional arguments which will be filled with user-specified data as above.
      \code{dmeasure} must return a single numeric value, the pdf of \code{y} given \code{x} at time \code{t}.
      If \code{log=TRUE}, then the log of the pdf is returned.
    }
    \item{\code{skeleton}}{
      If \code{skeleton} is an R function, it must have at least the arguments \code{x}, \code{t}, \code{params}, and \code{\dots}.
      \code{x} is a numeric vector containing the coordinates of a point in state space at which evaluation of the skeleton is desired.
      \code{t} is a numeric value giving the time at which evaluation of the skeleton is desired.
      Of course, these will be irrelevant in the case of an autonomous skeleton.
      \code{params} is a numeric vector holding the parameters.
      The optional argument \code{covars} is a numeric vector containing the values of the covariates at the time \code{t}.
      \code{covars} will have one value for each column of the covariate table specified when the \code{pomp} object was created.
      \code{covars} is constructed from the covariate table (see \code{covar}, below) by interpolation.
      \code{skeleton} may take additional arguments, which will be filled, as above, with user-specified data.
      \code{skeleton} must return a numeric vector of the same length as \code{x}.
      The return value is interpreted as the vectorfield (if the dynamical system is continuous) or the value of the map (if the dynamical system is discrete), at the point \code{x} at time \code{t}.

      If \code{skeleton} is the name of a native routine, this routine must be of prototype \dQuote{pomp_skeleton} as defined in the header \dQuote{pomp.h} (see the \dQuote{include} directory in the installed package directory).
    }
    \item{\code{initializer}}{
      if provided, must have at least the arguments \code{params}, \code{t0}, and \code{\dots}.
      \code{params} is a named numeric vector of parameters.
      \code{t0} will be the time at which initial conditions are desired.
      \code{initializer} must return a named numeric vector of initial states.
    }
  }
}
\examples{
## For examples, see the "Introduction to pomp"
## and "Advanced topics in pomp" documents,
## the demos, and the examples provided with the package, e.g.:
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
\section{Dynamically linked codes}{
  \code{pompLoad} and \code{pompUnload} cause any compiled codes associated with \code{object} to be dynamically linked or unlinked, respectively.
  When 
  \code{Csnippet}s are used in the construction of a \code{pomp} object, the resulting shared-object library is dynamically loaded (linked) at the time of its construction.
  A subsequent call to \code{pompUnload} will unload the shared object library;
  it can afterward be reloaded via a call to \code{pompLoad}.
}
\section{Warning}{
  Some error checking is done by \code{pomp}, but complete error checking is impossible.  
  If the user-specified functions do not conform to the above specifications (see Details), then the results may be invalid.
  In particular, if both \code{rmeasure} and \code{dmeasure} are specified, the user should verify that these two functions correspond to the same probability distribution and if \code{skeleton} is specified, the user is responsible for verifying that it corresponds to the true deterministic skeleton of the model.
  Each \pkg{pomp}-package algorithm uses some subset of the five basic functions (\code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, \code{skeleton}).
  If an algorithm requires a component that has not been specified, an informative error will be generated.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{
  \link{pomp methods},
  \link{pomp low-level interface},
  \link{process model plugins}
 }
\keyword{models}
\keyword{ts}
