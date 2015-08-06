\name{NEWS}
\title{News for package `pomp'}
\section{Changes in \pkg{pomp} version 0.78-1}{
  \itemize{
    \item There is a new, adaptive MCMC proposal function, \code{mvn.rw.adaptive} that can be used in \code{pmcmc}.
  }
}
\section{Changes in \pkg{pomp} version 0.77-1}{
  \itemize{
    \item There is a new, adaptive MCMC proposal function, \code{mvn.rw.adaptive} that can be used in \code{pmcmc}. 
    \item The deprecated \code{transform.params} argument from \code{nlf} has been removed.
    \item The deprecated \code{pars} and \code{rw.sd} arguments from \code{abc} and \code{pmcmc} have been removed.
    \item The deprecated \code{data.array} method has been removed.
  }
}
\section{Changes in \pkg{pomp} version 0.76-2}{
  \itemize{
    \item The default value of \code{max.fail} for \code{pmcmc} is now \code{Inf}.
    That is, by default \code{pmcmc} will not halt (but will issue a warning) when filtering failure occurs (i.e., no particles are consistent with the data, as defined by the likelihood threshold \code{tol}).
  }
}
\section{Changes in \pkg{pomp} version 0.76-1}{
  \itemize{
    \item \code{logmeanexp} has used a delta-method estimate of the variance in \eqn{log(mean(exp(x)))} when \code{se = TRUE}.
    This gave misleadingly small values when the variance in \eqn{exp(x)} was large.
    The new \code{logmeanexp} uses a jackknife estimate of the variance.
    This is still conservative, in that it is biased upwards, but is more robust.
    \item \pkg{pomp} now includes \pkg{magrittr}, \pkg{plyr}, \pkg{reshape2}, \pkg{ggplot2}, and \pkg{knitr} in its \code{Suggests} field:
    these packages are available for use in \code{pompExample}s and package tests.
    \item There is now a \code{coef} method for \code{mifList} and \code{mif2List} objects.
  }
}
\section{Changes in \pkg{pomp} version 0.75-1}{
  \itemize{
    \item By default, when \code{Csnippet}s are used, the source code is now stored internally to the \code{pomp} object.
    When a \code{pomp} object is used, a test is first performed to see if the needed shared-object file exists.
    If it does not, the source code is written to a file in the session's temporary directory, which is then compiled and loaded.
    This feature allows \code{pomp} objects to be stored and reused across \R sessions.
    The name of the file is constructed using a hash of the file's contents.
  }
}
\section{Changes in \pkg{pomp} version 0.74-1}{
  \itemize{
    \item The functions \code{bake} and \code{stew} have arguments that allow the control of the RNG.
    These functions can now set the state (and kind) of the RNG and restore it when they return.
    \item The new function \code{freeze} allows control of the RNG similarly.
    \item The \code{seed} argument of \code{bsmc} and \code{bsmc2} has been removed.
    Its use now generates a warning, stating that it has no effect.
    \item The \code{seed} argument to \code{pfilter} is now deprecated and will be removed soon.
    Equivalent functionality is provided via \code{freeze}, \code{bake}, or \code{stew}.
  }
}
\section{Changes in \pkg{pomp} version 0.73-1}{
  \itemize{
    \item The new functions \code{bake} and \code{stew} faciliate storing expensive calculations for later retrieval in a transparent fashion.
    \item \code{pompBuilder}, which has not been actively maintained for time due to the fact that its facilities are all provided by \code{pomp} via \code{Csnippet}s, is now no longer exported.
  }
}
\section{Changes in \pkg{pomp} version 0.72-1}{
  \itemize{
    \item \code{pmcmc} and \code{pfilter} now have the capability of saving filtered trajectories.
    These can be extracted using the new method \code{filter.traj}.
    The principal use will be in conjunction with \code{pmcmc}, where, upon convergence to the posterior, samples from the filtered trajectories are draws from the posterior \eqn{P[x_{1:T} | y_{1:T}]}{P[x[1:T] | y[1:T]]}.
    Thanks to Sebastian Funk for initiating this development.
  }
}
\section{Changes in \pkg{pomp} version 0.71-1}{
  \itemize{
    \item \code{partrans} and \code{parmat} have been jiggered to make them play better with \code{reshape2::melt}.
  }
}
\section{Changes in \pkg{pomp} version 0.70-1}{
  \itemize{
    \item Some of the underpinnings of \code{nlf} have been streamlined for speed.
  }
}
\section{Changes in \pkg{pomp} version 0.69-2}{
  \itemize{
    \item Measles data from England and Wales, kindly provided by Ben Bolker, are now included as \code{ewmeas} and \code{ewcitmeas}.
  }
}
\section{Changes in \pkg{pomp} version 0.68-2}{
  \itemize{
    \item When using \code{Csnippet}s, by default, the C codes and shared-object files are stored in the R session's temporary directory.
    One can override this behavior by setting the global \code{pomp.cache} option.
  }
}
\section{Changes in \pkg{pomp} version 0.67-1}{
  \itemize{
    \item A \code{pomp} object's initializer can now be specified as a Csnippet.
    \item The default initializer is now considerably faster.
  }
}
\section{Changes in \pkg{pomp} version 0.66-4}{
  \itemize{
    \item A full-featured version of IF2, an improved iterated filtering algorithm, is now available as \code{mif2}.
    This allows a more general structure for the random perturbations applied to parameters and is generally more efficient and stable than IF1.
  }
}
\section{Changes in \pkg{pomp} version 0.65-2}{
  \itemize{
    \item In the \code{pomp} constructor, unrecognized arguments are stored internally for use with user-supplied model functions.
    As of this version, an informative message is generated when this happens.
    Previously, this silent behavior was a source of confusion for some users.
  }
}
\section{Changes in \pkg{pomp} version 0.65-1}{
  \itemize{
    \item The unfortunately named \code{pomp} arguments \code{parameter.transform} and \code{parameter.inv.transform} have been deprecated in favor of the more mnemonic \code{toEstimationScale} and \code{fromEstimationScale}.
    Using the old arguments will result in the same behavior as previously, but a warning will be generated.
  }
}
\section{Changes in \pkg{pomp} version 0.64-1}{
  \itemize{
    \item The \code{mif} option \code{cooling.fraction} has been replaced by the more descriptive \code{cooling.fraction.50}.
    Recall that this parameter specifies the fraction of random-walk remaining after 50 \code{mif} iterations.
  }
}
\section{Changes in \pkg{pomp} version 0.63-6}{
  \itemize{
    \item Fixed bug arising when using \code{pomp} with \code{Csnippets} under windows.
    Thanks to Hamish McCallum for reporting it.
  }
}
\section{Changes in \pkg{pomp} version 0.63-3}{
  \itemize{
    \item Improvements to documentation.
  }
}
\section{Changes in \pkg{pomp} version 0.63-1}{
  \itemize{
    \item The optimization packages \pkg{nloptr} and \pkg{subplex} are no longer dependencies of \pkg{pomp} that will be loaded at the time \pkg{pomp} is.
    Instead, they are imported as needed.
  }
}
\section{Changes in \pkg{pomp} version 0.62-2}{
  \itemize{
    \item Bug fix in \code{partrans} when transformations are implemented using \code{Csnippet}s.
    \item Remove the \code{verhulst} example.
  }
}
\section{Changes in \pkg{pomp} version 0.62-1}{
  \itemize{
    \item \code{mif} now has the default value \code{var.factor=1}.
    Since \code{var.factor} has long been a required parameter, this change is backward compatible.
    \item The long-deprecated \code{cooling.factor} option to \code{mif} has now been removed.
  }
}
\section{Changes in \pkg{pomp} version 0.61-3}{
  \itemize{
    \item \code{mif}'s \code{pars} argument, which is redundant, is now deprecated.
    In a future release, this argument will be removed.
  }
}
\section{Changes in \pkg{pomp} version 0.61-1}{
  \itemize{
    \item When using \code{Csnippet}s or \code{pompBuilder}, it is no longer necessary to manually load and unload the dynamically loadable library.
    This is now handled automatically.
    Accordingly, the \code{link} argument of \code{pompBuilder} has been removed.
    The \code{pompLoad} and \code{pompUnload} functions introduced in 0.58-1 remain available, though they should rarely be needed.
  }
}
\section{Changes in \pkg{pomp} version 0.60-1}{
  \itemize{
    \item \code{pmcmc} and \code{abc} can now use arbitrary symmetric proposal distributions via the \code{proposal} argument.
    For the moment, these are constrained to be symmetric.
    Two new functions, \code{mvn.diag.rw} and \code{mvn.rw} generate suitable proposal functions.
    The first generates a multivariate normal random-walk proposal with diagonal variance-covariance matrix; this duplicates the old behavior of both \code{abc} and \code{pmcmc}.
    The second, \code{mvn.rw}, corresponds to a multivariate normal random-walk proposal with arbitrary variance-covariance matrix.
    \item In \code{pmcmc} and \code{abc}, the arguments \code{pars} and \code{rw.sd} are now unneeded and have been deprecated.
    Use of \code{rw.sd} will generate a warning and result in behavior equivalent to choosing \code{proposal=mvn.diag.rw(rw.sd)}.
    Use of \code{pars} will be ignored, with a warning.
    \item In \code{nlf}, the \code{transform.params} argument is now deprecated;
    use instead the \code{transform} argument, as in the other inference methods.
  }
}
\section{Changes in \pkg{pomp} version 0.59-1}{
  \itemize{
    \item A bug in \code{spect.match} has been fixed.
    Thanks to Karsten Hempel for finding it and contributing the fix.
  }
}
\section{Changes in \pkg{pomp} version 0.58-1}{
  \itemize{
    \item Compiler messages are now suppressed by default when \code{Csnippet}s are used in construction of a \code{pomp} object.
    These can be enabled by setting \code{options(verbose=TRUE)}.
  }
}
\section{Changes in \pkg{pomp} version 0.58-1}{
  \itemize{
    \item A \code{pomp} object created with \code{Csnippet}s or \code{pompBuilder} will now remember its associated shared-object file, which can be unloaded and reloaded using the new \code{pompUnload} and \code{pompLoad} methods.
  }
}
\section{Changes in \pkg{pomp} version 0.57-1}{
  \itemize{
    \item \code{simulate} has a new argument, \code{include.data}, that when used in conjuncion with \code{as.data.frame=TRUE} will include the actual data in the resulting data frame along with the simulations.
  }
}
\section{Changes in \pkg{pomp} version 0.56-1}{
  \itemize{
    \item Revamped the \code{pompExample} function.
    A search path for example directories is now stored in global option "pomp.examples".
  }
}
\section{Changes in \pkg{pomp} version 0.55-1}{
  \itemize{
    \item New \code{values} method extracts simulated probe values on \code{probed.pomp} object.
    \item Better documentation in the \file{pomp.h} header file.
    \item Vignettes now feature \code{bsmc2} instead of \code{bsmc}.
  }
}
\section{Changes in \pkg{pomp} version 0.54-1}{
  \itemize{
    \item A modified version of the Liu and West (2001) algorithm is included as \code{bsmc2}.
    \item By default, when \code{B} is of class \code{bsmcd.pomp}, \code{plot(B)} ignores fixed parameters.
  }
}
\section{Changes in \pkg{pomp} version 0.53-6}{
  \itemize{
    \item A bug having to do with paths to the temporary C files, encountered when using \code{Csnippet}s under windows, has been fixed.
    \item Package demos are now checked by \code{R CMD check} on all platforms.
  }
}
\section{Changes in \pkg{pomp} version 0.53-1}{
  \itemize{
    \item \code{nlf} now returns an S4 object of class \code{nlfd.pomp} with a \code{logLik} method for extracting the log quasi likelihood and a \code{$} method for extracting arbitrary components.
    \item The \code{transform} argument (for providing a function to transform the data) has been removed in favor of \code{transform.data}.
  }
}
\section{Changes in \pkg{pomp} version 0.52-1}{
  \itemize{
    \item The new \code{mifList} class facilitates approaches based on multiple \code{mif} runs.
    The \code{c} method constructs \code{mifList}s from \code{mif}s.
    The \code{plot} method produces diagnostic plots.
    \item \code{compare.mif} is now deprecated in favor of \code{plot} applied to a \code{mif} or \code{mifList}.
  }
}
\section{Changes in \pkg{pomp} version 0.51-3}{
  \itemize{
    \item \code{transform} argument for \code{pmcmc} and \code{abc} has been removed.
  }
}
\section{Changes in \pkg{pomp} version 0.51-2}{
  \itemize{
    \item More informative error messages are generated when a needed variable cannot be found.
    \item Features that facilitate multiple ABC chains have been added.
    There is a \code{c} method to combine \code{abc} objects into an \code{abcList} object and to combine \code{abcList} objects.
    Running \code{conv.rec} on an \code{abc} or \code{abcList} object produces \pkg{coda} \code{mcmc} and \code{mcmc.list} objects, respectively.
    Diagnostic plots are produced by running \code{plot} on such objects.
    \item Basic structure of \code{pmcmcList} and \code{abcList} objects is changed: these objects now inherit directly from \code{list}.
  }
}
\section{Changes in \pkg{pomp} version 0.51-1}{
  \itemize{
    \item Package \pkg{coda} is now a dependency of \pkg{pomp}.
    This package provides various generic MCMC diagnostics.
    \item Features that facilitate multiple PMCMC chains have been added.
    There is a \code{c} method to combine \code{pmcmc} objects into a \code{pmcmcList} object and to combine \code{pmcmcList} objects.
    Running \code{conv.rec} on a \code{pmcmc} or \code{pmcmcList} object now produces \pkg{coda} \code{mcmc} and \code{mcmc.list} objects, respectively.
    Diagnostic plots are produced by running \code{plot} on such objects.
    \item The \code{compare.pmcmc} and \code{compare.abc} methods have been removed as obsolete.
  }
}
\section{Changes in \pkg{pomp} version 0.50-9}{
  \itemize{
    \item Fix bug in \code{compare.pmcmc}.
  }
}
\section{Changes in \pkg{pomp} version 0.50-8}{
  \itemize{
    \item Introduced new \code{compare.pmcmc} and \code{compare.abc} functions, comparable to \code{compare.mif}.
  }
}
\section{Changes in \pkg{pomp} version 0.50-6}{
  \itemize{
    \item The package manual and tutorials are no longer included with the package source.
    Instead, links are provided to the \pkg{pomp} website.
    \item The package help pages have been revisited for accuracy and completeness.
  }
}
\section{Changes in \pkg{pomp} version 0.50-1}{
  \itemize{
    \item \code{pomp} can now take snippets of C code directly as definitions of the various basic functions.
    This is accomplished through the use of a new \code{Csnippet} function.
    See the \sQuote{sir} demo for a demonstration.
    \item It is now possible to specify global variables and typedefs using the \code{globals} argument to \code{pompBuilder} and \code{pomp}.
    See the \sQuote{sir} demo for a demonstration.
    Thanks to Carl Pearson for the suggestion.
  }
}
\section{Changes in \pkg{pomp} version 0.49-1}{
  \itemize{
    \item A new method, \code{probe.match.objfun}, constructs objective functions suitable for use in optimizers.
    \item \code{pomp} now depends on \code{nloptr}, which provides a suite of optimization algorithms.
    This package can now be used in various methods for optimization of an objective function.
    \item New inline C functions \code{to_log_barycentric} and \code{from_log_barycentric} are provided in \file{pomp.h} to facilitate log-barycentric transformations.
    These have proven very useful in dealing with parameters constrained to sum to one (e.g., initial conditions of compartmental models).
    \item Bug fix in \code{periodic.bspline.basis}.
    \item Package vignettes are no longer available via \code{vignette()} but can be reached from the help index.
  }
}
\section{Changes in \pkg{pomp} version 0.48-3}{
  \itemize{
    \item Correct a bug in \code{abc} to do with parameter transformation.
  }
}
\section{Changes in \pkg{pomp} version 0.48-1}{
  \itemize{
    \item Create new \code{dprior} and \code{rprior} slots for \code{pomp} objects.
    These will be used by the Bayesian methods (currently \code{abc}, \code{bsmc}, and \code{pmcmc}).
    The \code{hyperparams} argument to these methods has been removed;
    hyperparameters are now treated just as any others.
    \item \code{bsmc} can now be called so that the particles are drawn from the prior distribution (as set by the \code{rprior} slot).
    \item Some of the longer tests will only be run if an environment variable is set.
  }
}
\section{Changes in \pkg{pomp} version 0.47-4}{
  \itemize{
    \item Revisit \code{dacca} bug.  Set negative compartments to zero along with compartments immediately downstream.
  }
}
\section{Changes in \pkg{pomp} version 0.47-3}{
  \itemize{
    \item Fix bug that arises only very occasionally in the \code{dacca} cholera example.
    \item Modify \code{dacca} model so that rare positivity violations are punished in a different way.
    This modifies the formal model slightly.
  }
}
\section{Changes in \pkg{pomp} version 0.47-2}{
  \itemize{
    \item By default, \code{pompBuilder} now creates all files in the temporary directory.
  }
}
\section{Changes in \pkg{pomp} version 0.47-1}{
  \itemize{
    \item \code{abc} implements Approximate Bayesian Computation for pomp models.
  }
}
\section{Changes in \pkg{pomp} version 0.46-1}{
  \itemize{
    \item \code{pompExample} now has an optional argument, \code{envir}, determining which environment the \code{pomp} object will be loaded into.
    \item \code{logmeanexp} now has an optional argument, \code{se}, allowing for computation of an approximate standard error.
  }
}
\section{Changes in \pkg{pomp} version 0.45-8}{
  \itemize{
    \item bug fix in \code{gompertz} example.
  }
}
\section{Changes in \pkg{pomp} version 0.45-7}{
  \itemize{
    \item bug fix in \code{pmcmc}: Metropolis-Hastings ratio was incorrect except for flat priors!
  }
}
\section{Changes in \pkg{pomp} version 0.45-6}{
  \itemize{
    \item a fix so that \code{pompBuilder} will find \file{pomp.h} header file on Windows machines.  Thanks to Dave Hayman for finding the problem.
    \item new \code{logmeanexp} function.
    \item some reorganization of the source package structure.
  }
}
\section{Changes in \pkg{pomp} version 0.45-4}{
  \itemize{
    \item changes in the way \code{reulermultinom}, \code{deulermultinom}, \code{dot_product} are exported to other packages.  Rather than just being exported as linkables, these are now also defined as static inline functions in the \file{pomp.h} header.
  }
}
\section{Changes in \pkg{pomp} version 0.45-3}{
  \itemize{
    \item Fix bug with \code{continue} and method \code{"mif2"}.
  }
}
\section{Changes in \pkg{pomp} version 0.45-2}{
  \itemize{
    \item new method to coerce \code{pfilterd.pomp} objects to data-frames.
  }
}
\section{Changes in \pkg{pomp} version 0.45-1}{
  \itemize{
    \item \code{profileDesign} can now handle variables of mixed type.
  }
}
\section{Changes in \pkg{pomp} version 0.44-1}{
  \itemize{
    \item New \code{"mif2"} option for \code{mif}.
    \item The default behavior for \code{mif} and \code{pfilter} with respect to filtering failures has been changed.
    Before, the default behavior has been to stop with an error on the first filtering failure (\code{max.fail=0}).
    Now, the default is \code{max.fail=Inf}, i.e., an error is never triggered.
    \item The implementation of MIF cooling schedules has been changed to make it more general.
    The cooling schedule is now specified by a \code{type} and a \code{fraction}.
    Currently, supported \code{cooling.type}s include \code{geometric} (the old behavior) and \code{hyperbolic}, i.e., a \eqn{1/(1+n)} schedule.
    The \code{cooling.fraction} argument specifies the cooling at 50 iterations.
    That is, if \code{s} is the intensity of the random-walk perturbation to parameters at the first iteration (\code{rw.sd}), then the intensity at iteration 50 is \code{s*cooling.fraction}.
    \item Remove all \code{data()}-loadable \code{pomp} objects.
    To load the prebuilt example \code{pomp} objects from previous versions, use the new \code{pompExample} function.
    E.g., instead of \code{data(euler.sir)}, do \code{pompExample("euler.sir")}.
  }
}
\section{Changes in \pkg{pomp} version 0.43-4}{
  \itemize{
    \item Clean up tests.
  }
}
\section{Changes in \pkg{pomp} version 0.43-3}{
  \itemize{
    \item Fix bug in \code{pompBuilder} on mac.
    \item Fix bug in \code{dmeasure} when \code{log=TRUE}.  Bug located by R. Reiner.
    \item Add \code{eff.sample.size} and \code{cond.logLik} methods for \code{pfilterd.pomp} objects.
  }
}
\section{Changes in \pkg{pomp} version 0.43-1}{
  \itemize{
    \item Modify \code{traj.match.objfun} to allow arbitrary parameters to be passed to the ODE integrator.
    \item Add \code{params} argument to \code{pomp} so that parameters can be specified at the time the \code{pomp} object is created.
  }
}
\section{Changes in \pkg{pomp} version 0.42-6}{
  \itemize{
    \item Improve error checking in \code{pomp}.
    \item Add support for covariates to \code{pompBuilder}.
    \item Add native parameter transformation routines to the demos.
    \item Change the warning behavior in \code{lookup_table} so that warnings are not generated unless true extrapolation is performed (i.e., eliminate the warnings at the boundary).
    \item Some small modifications to the \file{pomp.h} header file.
  }
}
\section{Changes in \pkg{pomp} version 0.42-4}{
  \itemize{
    \item When \code{trajectory} calls \code{ode}, there is a need to preserve information across calls to the vectorfield-evaluating function.
    The code has been modified to do this in a cleaner way.
  }
}
\section{Changes in \pkg{pomp} version 0.42-2}{
  \itemize{
    \item Several minor bugs have been fixed.
    See the ChangeLog for details.
  }
}
\section{Changes in \pkg{pomp} version 0.42-1}{
  \itemize{
    \item An EXPERIMENTAL facility for constructing \code{pomp} objects out of native C routines is now included (pompBuilder).
    This facility is being actively developed and future changes may not be backward compatible.
    \item A new facility allowing access to the \code{userdata} slot of a \code{pomp} object from compiled \code{rprocess}, \code{rmeasure}, \code{dprocess}, \code{dmeasure}, \code{skeleton}, and \code{partrans} codes is now available.
    Calls to the new C routines \code{get_pomp_userdata}, \code{get_pomp_userdata_int}, and \code{get_pomp_userdata_double} allow retrieval of these elements.
    \item Some of the \code{data()}-loadable examples have been reworked to make use of the above facility.
    The parameterization of these examples has changed.
    These changes are not backward compatible: codes that depend on the \code{data()}-loadable examples may be broken.
  }
}
\section{Changes in \pkg{pomp} version 0.41-8}{
  \itemize{
    \item A demonstration of \code{pompBuilder} has been put into the \dQuote{advanced topics} vignette.
  }
}
\section{Changes in \pkg{pomp} version 0.41-7}{
  \itemize{
    \item A bug in the \code{blowflies} example has been fixed.
    Thanks to Greg Minshall for discovering it.
  }
}
\section{Changes in \pkg{pomp} version 0.41-6}{
  \itemize{
    \item The \code{gompertz} example parameter transformations have been changed.
    No longer are the names of the parameter vector changed in the transformation.
    This change is not backward-compatible, but only codes using this example are affected.
    \item The \code{euler.sir} and \code{gillespie.sir} examples have been changed.
    The transmission rate \eqn{\beta(t)} is now the arithmetic sum of the seasonality basis functions.
    Before, it was the geometric sum.
    \eqn{R_0} is now given by the arithmetic average of the \eqn{\beta} parameters divided by \eqn{\gamma+\mu}.
    This change is not backward-compatible, but only codes using these examples are affected.
  }
}
\section{Changes in \pkg{pomp} version 0.41-5}{
  \itemize{
    \item An experimental facility for constructing \code{pomp} objects with native C routines is now included.
  }
}
\section{Changes in \pkg{pomp} version 0.41-4}{
  \itemize{
    \item The \code{blowflies}, \code{euler.sir}, \code{gillespie.sir}, \code{bbs}, and \code{ricker} \code{data()}-loadable examples have been changed to make the parameterization simpler and more natural.
    This change is not backward-compatible, but only codes using this example are affected.
  }
}
\section{Changes in \pkg{pomp} version 0.41-3}{
  \itemize{
    \item In \code{trajectory}, all vectorfield and map evaluation is now done in C for speed.
    For continuous-time dynamical systems, this gives an approximately 4-fold speedup.
    \item Demos have been created to show some examples.
  }
}
\section{Changes in \pkg{pomp} version 0.41-2}{
  \itemize{
    \item A segfault bug in \code{bbs} example has been fixed.
    \item The posterior medians (not means) are now stored in the \code{params} slot of the \code{bsmcd.pomp} object.
  }
}
\section{Changes in \pkg{pomp} version 0.41-1}{
  \itemize{
    \item New arguments in \code{mif}, \code{nlf}, \code{bsmc}, \code{pmcmc}, \code{probe-match}, and \code{traj-match} allow the estimation to be done on a transformed parameter space.
    When \code{transform=TRUE} in these commands (\code{transform.params=TRUE} for \code{nlf}), estimation is performed on the transformed parameter space.
    This is described and demonstrated in the \dQuote{intro_to_pomp} vignette.
    \item The \code{data()}-loadable examples have been re-implemented to make use of the above-mentioned facility.
    Note that this new functionality makes it unnecessary to un-transform model parameters within the user-specified \code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, \code{skeleton}, and \code{initializer} codes.
    This change is not backward-compatible, but only codes using these \code{data()}-loadable example are affected.
    \item The Bayesian sequential Monte Carlo command \code{bsmc} now returns not a list but an object of class \code{bsmcd.pomp}.
    An experimental \code{plot} method for objects of this class now exists.
    Also, the parameter posterior means are now stored in the \code{params} slot of the \code{bsmcd.pomp} object:
    access them with the \code{coef} command as usual.
    \item A new example, using data from an influenza outbreak in a British boarding school and an SIR model, has been included.
    Do \code{data(bbs)} to load it.
  }
}
\section{Changes in \pkg{pomp} version 0.40-9}{
  \itemize{
    \item Setting the new argument \code{as.data.frame} to \code{TRUE} in \code{simulate} and \code{trajectory} causes the results to be returned as a data-frame.
  }
}
\section{Changes in \pkg{pomp} version 0.40-8}{
  \itemize{
    \item A new method \code{partrans} allows transformation of vectors or matrices of parameters.
    The parameter transformations have been pushed into C for speed.
    It is possible to specify native parameter transformation routines in addition to C functions, but this last facility has not yet been extensively tested.
    A new slot \code{has.trans} has been introduced into the \code{pomp} class, and the types of slots \code{par.trans} and \code{par.untrans} have changed.
  }
}
\section{Changes in \pkg{pomp} version 0.40-7}{
  \itemize{
    \item \code{parmat} now gracefully handles the case when \code{params} is already a matrix.
    \item Error in the documentation of \code{pfilter} repaired.
    Thanks to Alex Smith for pointing this out.
    \item Improvements in the documentation of the Nicholson blowflies data \code{data(blowflies)}.
  }
}
\section{Changes in \pkg{pomp} version 0.40-6}{
  \itemize{
    \item When a \code{pomp} object is updated using \code{pomp}, the \code{params} slot will now be copied over (appropriately transformed).
    Before, it was simply dropped.
    \item The \code{gillespie.sir} example has been changed to make the reporting rate equal to \code{rho}.
    Before, it was fixed at 10\%.
    \item A new error-trap has been put into \code{probe} and \code{probe.match} for the case \code{nsim < nprobes}.
    Thanks to Dani Jonas for finding this error.
  }
}
\section{Changes in \pkg{pomp} version 0.40-5}{
  \itemize{
    \item More informative error messages when dimension of state space or data space disagree internally.
    \item The \code{weighted} argument to \code{mif} and \code{continue} is now deprecated in favor of a new argument \code{method}.
    The old \code{weighted=TRUE} corresponds to \code{method="mif"} while the old \code{weighted=FALSE} corresponds to \code{method="unweighted"}.
  }
}
\section{Changes in \pkg{pomp} version 0.40-4}{
  \itemize{
    \item New functions \code{traj.match.objfun} and \code{probe.match.objfun} have been added.
    These functions construct functions of one argument suitable for use as objective functions in \code{optim}-like optimizers (iincluding \code{subplex} and \code{sannbox}).
    Minimizing these functions solves the trajectory-matching problem and probe-matching problem, respectively.
    This allows the user much greater flexibility in choice of optimization algorithm than is afforded by \code{traj.match} and \code{probe.match}.
    \item The \code{sannbox} optimizer, which performs simulated annealing with box constraints, is now exported and available for general use.
  }
}
\section{Changes in \pkg{pomp} version 0.40-3}{
  \itemize{
    \item The \code{bsmc} method has been improved, due to the contributions of Pierre Jacob.
    Specifically, \code{bsmc} no longer reports a log-likelihood (which it never actually computed anyway) but a log-evidence.
    The latter computation was supplied by Pierre Jacob.
    \item A new helper function, \code{exp2geom_rate_correction} has been added to the \file{pomp.h} file.
    This function computes the rate \eqn{r} such that if \eqn{N \sim geometric(prob=1-\exp(-r dt))} and \eqn{T \sim exponential(rate=R)}, then \eqn{E[N dt] = E[T]}.
    This is useful in approximating a continuous-time death process by a discrete time (Euler) process.
    In particular, in such a case, \code{T} is the waiting time to death in the former and \code{N*dt} is the waiting time to death in the latter.
    An Euler binomial or multinomial process with rate \code{r=exp2geom_rate_correction(R,dt)} will have the same mean waiting time as the exponential process with rate \code{R}.
    Thanks to Sebastien Ballesteros for suggesting this.
    \item A new helper function, \code{rgammawn} has been added, with both R and C interfaces.
    This function draws an increment of a Gamma white-noise process with a given intensity.
    Specifically, \code{dw=rgammwn(sigma,dt)} is Gamma distributed with mean \code{dt} and variance \code{sigma^2*dt}.
    In this case, \code{mu*dw/dt} is suitable for use as a random rate in an Euler-multinomial process.
    Thanks to Sebastien Ballesteros for suggesting this.
  }
}
\section{Changes in \pkg{pomp} version 0.40-2}{
  \itemize{
    \item A bug to do with computation of the number of steps needed in discrete-time simulation and trajectory computations has been fixed.
    This bug was introduced in version 0.40-1.
  }
}
\section{Changes in \pkg{pomp} version 0.40-1}{
  \itemize{
    \item It is now possible for discrete-time determistic skeletons to have a general step-size.
    Before, the stepsize was always forced to be 1.
    The new argument \code{skelmap.delta.t} in \code{pomp} (default value 1) sets this step-size.
    NB: in implementing this feature, a new slot has been added to the basic \code{pomp} class.
    \code{pomp} objects created with versions <0.40-1 will not therefore work with this version and must be re-created.
    \item The deterministic skeleton trajectory functions now look for variables named in \code{zeronames} and treat these variables as accumulators.
    This brings the behavior of the deterministic skeleton functions into line with \code{rprocess}.
    A new section of the \dQuote{advanced topics} vignette discusses accumulator variables and \code{zeronames}.
    \item Include the user manual PDF among the vignettes.
    \item The new function \code{parmat} creates a suitable matrix of parameters in one line.
  }
}
\section{Changes in \pkg{pomp} version 0.39-5}{
  \itemize{
    \item Bug fix in \code{bsmc} that arises when the state space is 1-D.
    Thanks to Arnab Bhattacharya for finding the bug.
    \item Some inconsistencies in the arguments of \code{bsmc} have been fixed.
  }
}
\section{Changes in \pkg{pomp} version 0.39-4}{
  \itemize{
    \item Minor corrections to the documentation for \code{spect}, \code{probe}, and \code{traj.match}.
    \item Fix infelicity in \code{traj.match} and \code{nlf} that arises when the \code{data} slot has no rownames.
    Rownames are not strictly required for the \code{data} slot.
    \item The \code{sannbox} algorithm now should work even when the objective function throws a non-finite value.
    The behavior now is to replace non-finite values with a large number.
    This behavior follows that of the \code{SANN} method in \code{optim}.
  }
}
\section{Changes in \pkg{pomp} version 0.39-3}{
  \itemize{
    \item The help files for \code{pfilter} and \code{mif} have been updated to explain how to use a variable number of particles in the particle filtering algorithm.
    \item Inside \code{init.state}, a check has been added to make sure that the user's initializer function returns a vector of uniform size, as per the algorithms' assumptions.
    Thanks to Micaela Martinez-Bakker for catching this.
    \item \code{mif} has been modified so as to allow estimation of initial-value parameters alone via fixed-lag smoothing.
  }
}
\section{Changes in \pkg{pomp} version 0.39-2}{
  \itemize{
    \item When \code{po} is a \code{pomp}-class object with covariates, \code{as.data.frame(po)} and \code{as(po,"data.frame")} now contain columns for the covariates.
    The latter are interpolated, if necessary, at the observation times.
    \item Bug fix: when \code{pfilter} is called with a single parameter vector supplied in the \code{params} argument, these parameters are copied into the \code{params} slot of the resulting \code{pfilterd.pomp} object.
    Before, this failed when the \code{params} slot of the \code{pomp} object was empty.
    A similar bug existed in \code{traj.match} and has likewise been fixed.
    \item To eliminate warnings and confusion during package check, \code{profile.design} has been replaced by \code{profileDesign}.
    For parallelism, \code{slice.design} and \code{sobol.design} have been replaced by \code{sliceDesign} and \code{sobolDesign}, respectively.
    The older forms are now deprecated and will be removed in a future release.
  }
}
\section{Changes in \pkg{pomp} version 0.39-1}{
  \itemize{
    \item New facilities for parameter transformation are provided.
    New optional arguments \code{parameter.transform} and \code{parameter.inv.transform} to \code{pomp} allow the user to specify transformations between natural and internal parameter scalings.
    A new option \code{transfom} to \code{coef} and \code{coef<-} allows the user to get and set parameters on the natural scale.
    See \code{pomp?coef} and the \dQuote{intro_to_pomp} vignette for details and examples.
    Similarly, the \code{conv.rec} method now takes an optional \code{transform} argument.
    \item Minor improvements to the documentation.
  }
}
\section{Changes in \pkg{pomp} version 0.38-5}{
  \itemize{
    \item \code{pfilter} will now give an error if ever a non-finite likelihood (dmeasure value) is returned.
    Before, errors were only generated when dmeasure returned NA.
    \item \code{spect} now returns frequencies reflective of the observation times.
    Before, the frequencies were given as fractions of the sampling frequency.
  }
}
\section{Changes in \pkg{pomp} version 0.38-4}{
  \itemize{
    \item Improvements to the documentation.
  }
}
\section{Changes in \pkg{pomp} version 0.38-3}{
  \itemize{
    \item Clarifications and minor improvements in the help pages.
  }
}
\section{Changes in \pkg{pomp} version 0.38-2}{
  \itemize{
    \item A new section describing and demonstrating nonlinear forecasting has been added to the \dQuote{intro to pomp} vignette.
    \item Fixed bug in \code{nlf}: desire for transient period (nconverge>0) was unrecognized.
  }
}
\section{Changes in \pkg{pomp} version 0.38-1}{
  \itemize{
    \item Fixed bug in \code{pomp} when \code{data} is of class \code{pomp}: userdata should be preserved.
    \item Improved discussion of SIR model in \dQuote{intro to pomp} vignette.
    \item In \code{pfilterd.pomp} objects, \code{saved.states} and \code{saved.params} slots are now length-ntimes lists of arrays.
    \item It is now possible to use variable numbers of particles in \code{pfilter}.
    The \code{Np} argument may be specified as a single number as before, or now as either a vector of numbers or a function.
  }
}
\section{Changes in \pkg{pomp} version 0.37-2}{
  \itemize{
    \item The \dQuote{advanced topics} vignette has been much augmented and improved.
    \item The \code{euler.sir} example (load it with \code{data(euler.sir)}) has a different measurement model: a discretized normal distribution with the same mean and variance as before.
  }
}
\section{Changes in \pkg{pomp} version 0.37-1}{
  \itemize{
    \item The arguments \code{skeleton.map} and \code{skeleton.vectorfield} are now deprecated.
    Specify a discrete-time skeleton using \code{skeleton.type="map"} and a continuous-time skeleton via \code{skeleton.type="vectorfield"}.
    \item \code{pomp} is now a generic function.
    Calling it on a \code{pomp} object copies the object; one can optionally modify some of the slots.
    \item Two new files \file{gompertz.R} and \file{gompertz.c} are in the \file{examples} directory.
    These contain codes needed to reproduce the Gompertz example in the \dQuote{intro_to_pomp} vignette.
    \item The warning messages to do with the default behaviors of \code{simulate} and \code{trajectory} with respect to the arguments \code{times} and \code{t0} have been removed.
    These warnings were introduced in version 0.34-1 to alert users to changes in the default behavior of \code{simulate} and \code{trajectory} introduced with that version.
    See the help (\code{pomp?simulate}, \code{pomp?trajectory}) for details.
    \item It is now possible to save both the states and particles portions of the filtering distributions computed in \code{pfilter}.
    This is accomplished by means of the logical \code{save.states} and \code{save.params} arguments to \code{pfilter}.
    The particles are stored in rank-3 arrays in the slots \code{saved.states} and \code{saved.params}.
    This results in the need for a change in the structure of \code{pfilterd.pomp} objects: note that \code{pfilterd.pomp} objects (including \code{mif} and \code{pmcmc} objects) created under pre-0.37-1 version of the package will not load under versions >= 0.37-1.
  }
}
\section{Changes in \pkg{pomp} version 0.36-7}{
  \itemize{
    \item \code{trajectory} now gives a more informative error when no skeleton is present
    \item \code{pred.mean}, \code{pred.var}, and \code{filter.mean} methods are now defined for \code{pfilterd.pomp} objects (and not just \code{mif} objects)
    \item A new \code{data()}-loadable example is included: do \code{?blowflies} to read about it.
  }
}
\section{Changes in \pkg{pomp} version 0.36-6}{
  \itemize{
    \item The B-spline basis construction functions \code{bspline.basis} and \code{periodic.bspline.basis} now take a \code{names} argument that allows one to give names to the basis functions.
  }
}
\section{Changes in \pkg{pomp} version 0.36-5}{
  \itemize{
    \item The \dQuote{intro_to_pomp} vignette has been updated.  The first part of the vignette centers on a new, more transparent first example.
    \item A bug in \code{pfilter} has been fixed.  This bug was introduced in version 0.32-6 (SVN revision 315).  Upon filtering failure, state variables were left un-updated.  Thanks to Karim Khader for reporting circumstances that led to the identification of this bug.
    \item There is now a method for \code{loglik} on \code{probed.pomp} objects.  This returns the synthetic likelihood.
    \item The box-constrained simulated annealing optimization method \code{sannbox} is now an option in \code{probe.match}.
    \item Calling \code{states} on a \code{pomp} object having no unobserved states now fails gracefully (returning NULL rather than generating an error).
    \item The rprocess plugin \code{discrete.time.sim} now takes an optional \code{delta.t} argument to allow for arbitrary (but still necessarily regular) spacing of observations.  By default, delta.t=1.
  }
}
\section{Changes in \pkg{pomp} version 0.36-1}{
  \itemize{
    \item There are two new classes: \code{traj.matched.pomp} and \code{pfilterd.pomp} to hold the results of \code{traj.match} and \code{pfilter} computations, respectively.  The old classes \code{mif} and \code{pmcmc} now inherit from \code{pfilterd.pomp} in a natural way.
    \item When \code{traj.match} is called, the returned \code{traj.matched.pomp} object has its \code{states} slot filled with the model trajectory at the fitted parameters.
    \item More optimization methods are now provided in \code{traj.match}.  These include the new algorithm \dQuote{sannbox}, an optionally box-constrained simulated annealing algorithm.
    \item \code{probe.match} is now an S4 method (like \code{pfilter}, \code{mif}, \code{pmcmc}, and \code{traj.match}) with methods for classes \code{pomp}, \code{probed.pomp}, and \code{probe.matched.pomp}.
    \item There is a change to the interface to \code{pfilter} in that \code{save.states} results in the filling of the \code{last.states} slot.  Before version 0.36-1, the \code{pfilter} returned a list.  The element \code{states} of that list corresponds to the slot \code{last.states} in version 0.36-1.  It was necessary to make this name-change in order to avoid a conflict with the \code{states} slot inherited from the \code{pomp} class.
  }
}
\section{Changes in \pkg{pomp} version 0.35-1}{
  \itemize{
    \item Added a facility for computation of \sQuote{synthetic likelihood} (Wood 2010) using a robust estimate of the covariance.  A new slot in objects of class \code{probed.pomp} holds this value.
    \item The synthetic likelihood is now the quantity that is maximized in \code{probe.match}.
    \item Added the argument \code{type} to \code{probe.ccf} so that either cross-covariance or cross-correlation can be returned.  \code{probe.ccf} and \code{probe.acf} now have comparable interfaces.
    \item Changed the way that \code{probe.acf} accepts lags.  Now an arbitrary set of lags can be selected.
    \item A bug in the new method \code{obs} has been fixed.
    \item A protection-stack overflow bug in \code{probe} has been fixed.
    \item The \code{seed} slot in \code{probed.pomp} and \code{spect.pomp} objects has been removed, along with unwanted behavior to do with it.
  }
}
\section{Changes in \pkg{pomp} version 0.34-2}{
  \itemize{
    \item Added a method to coerce \code{probed.pomp} objects to \code{data.frame} objects.
    \item Fix bug when \code{simulate} is run on objects that inherit from \code{pomp} but are not themselves pomps.
    \code{simulate} should return \code{pomp} objects when \code{obs=FALSE} and \code{states=FALSE}, regardless of what class \code{object} is.
    \item \code{simulate} and \code{rprocess} are accelerated somewhat through judicious use of \code{memcpy}.
  }
}
\section{Changes in \pkg{pomp} version 0.34-1}{
  \itemize{
    \item The default behaviors of \code{simulate} and \code{trajectory} with respect to the argument \code{times} has changed.
    There is a new argument, \code{t0}, which specifies the time at which the initial conditions obtain.
    The default for \code{t0} is \code{timezero(object)}.
    The default for \code{times} is now \code{time(object,t0=FALSE)} (it was \code{time(object,t0=TRUE)}).
    This change eliminates an infelicity in which \code{simulate} and \code{trajectory} returned undesired values at the initial time.
    Users can ensure that their codes will continue to function as intended by specifying the values of the \code{times} and \code{t0} arguments to these functions, thus removing dependence on the defaults.
    To aid users in  modifying their codes, a warning is now displayed whenever a potentially affected call is made.
    These warnings will be removed in a future release.
    See the documentation (\code{pomp?simulate}, \code{pomp?trajectory}) for a description of the new default behaviors.
    \item The \code{states} slot now holds values of the state process only at the times given in the \code{times} slot.
    This is a change from the earlier behavior in which the \code{states} slot sometimes included the value of the state process at \code{t0}.
    \item \code{simulate} has been re-implemented in C for speed.
    \item A new method \code{obs} for \code{pomp} objects returns the contents of the \code{data} slot.
    It is synonymous with \code{data.array}.
    \item A banner is now displayed on package attachment.
    \item \code{coef<-} now behaves like ordinary vector assignment.
    \item \code{probe} and \code{spect} now take an argument \code{params}.
    \item \code{probe.cov} and \code{probe.cor} have been removed in favor of \code{probe.acf} and \code{probe.ccf}.
    \item The functions \code{euler.simulate}, \code{onestep.simulate}, and \code{onestep.density}, deprecated since version 0.29-1, have been removed.
  }
}
\section{Changes in \pkg{pomp} version 0.33-1}{
  \itemize{
    \item Major improvements in speed in the probe-matching routines have been realized by coding the computationally-intensive portions of these calculations in C.
    \item New probes recommended by Simon Wood (Nature 466:1102-1104, 2010) have been implemented.
    In particular, the new function \code{probe.marginal} regresses the marginal distributions of actual and simulated data against a reference distribution, returning as probes the regression coefficients.
    \code{probe.marginal} and \code{probe.acf} have been coded in C for speed.
  }
}
\section{Changes in \pkg{pomp} version 0.32-1}{
  \itemize{
    \item pomp now includes parameter estimation methods based on probe-matching and power-spectrum-matching as explained in papers by Kendall et al. (Ecology, 80:1789--1805, 1999) and Reuman et al. (Proc. Nat'l. Acad. Sci., U.S.A. 103:18860--18865, 2006).
  }
}
\section{Changes in \pkg{pomp} version 0.31-2}{
  \itemize{
    \item new \code{window}, \code{timezero}, and \code{timezero<-} methods are available for manipulating \code{pomp} objects.
  }
}
\section{Changes in \pkg{pomp} version 0.31-1}{
  \itemize{
    \item pomp now includes an implementation of the particle Markov chain Monte Carlo algorithm of Andrieu et al. (J. R. Stat. Soc. B 72:1--33, 2010).
    This algorithm is implemented in the method \code{pmcmc}: see \code{help(pmcmc)} for details.
    Implementation primarily due to Ed Ionides.
  }
}
\section{Changes in \pkg{pomp} version 0.30-1}{
  \itemize{
    \item pomp now includes an implementation of the approximate Bayesian Sequential Monte Carlo algorithm of Liu & West (Liu, J. & West, M. (2001) Combining Parameter and State Estimation in Simulation-Based Filtering. In Doucet, A.; de Freitas, N. & Gordon, N. J. (eds.) Sequential Monte Carlo Methods in Practice, Springer, New York, pp. 197--224).
    This algorithm is implemented in the method \code{bsmc}: see \code{help(bsmc)} for details.
    Thanks to Matt Ferrari and Michael Lavine for the implementation.
  }
}
\section{Changes in \pkg{pomp} version 0.29-5}{
  \itemize{
    \item The cholera transmission model from King, A. A., Ionides, E. L., Pascual, M., and Bouma, M. J. (2008) Inapparent infections and cholera dynamics, Nature 454:877--880 together with data from the Dacca district of historic Bengal are now included.
    See \code{help(dacca)} for details.
  }
}
\section{Changes in \pkg{pomp} version 0.29-1}{
  \itemize{
    \item This is a major revision.
    \item Big improvements to the documentation, including worked examples in vignettes and improvements to the manual pages.
    \item Encoding models is now much easier through the use of \dQuote{plugin} modules.
    These are functions that write the \code{rprocess} and \code{dprocess} slots given minimal information on the structure of the model.
    \item A plugin for continuous-time Markov process models using the Gillespie algorithm is available.
    \item A plugin for discrete-time Markov processes is now available.
    \item Several example pomps are now loadable via the \code{data()} command.
    \item There are changes in the C interface to measurement models that are NOT BACKWARD COMPATIBLE.
    Existing user measurement model codes will need to be modified, but these changes are minor.
    \item The changes are documented in detail in the file \file{CHANGES_0.29-1.txt}.
  }
}
\section{Changes in \pkg{pomp} version 0.22-1}{
  \itemize{
    \item The nonlinear forecasting method of Kendall, Ellner, et al. is now implemented in the package.  This is an \dQuote{indirect inference} method using quasi-likelihood as an objective function.
  }
}
\section{Changes in \pkg{pomp} version 0.21-4}{
  \itemize{
    \item Improved documentation in \file{euler.Rd}.
  }
}
\section{Changes in \pkg{pomp} version 0.21-3}{
  \itemize{
    \item Trajectories of the deterministic skeleton can now be computed via the method \code{trajectory}.
    \item Bug fix in stochastic Euler algorithm.
  }
}
\section{Changes in \pkg{pomp} version 0.20-8}{
  \itemize{
    \item Better error handling in mif
    \item Fix inaccuracy in documentation of \code{pomp} arguments \code{rmeasure} and \code{dmeasure}
  }
}
\section{Changes in \pkg{pomp} version 0.20-7}{
  \itemize{
    \item Bug fix in convergence record storage
  }
}
\section{Changes in \pkg{pomp} version 0.20-6}{
  \itemize{
    \item Bug fix in \dQuote{intro_to_pomp} vignette
  }
}
\section{Changes in \pkg{pomp} version 0.20-5}{
  \itemize{
    \item Better error handling and improved documentation for \code{mif}
    \item Testing of examples in the \file{examples} directory, including compilation, linking, and running of dynamically-loadable library on unix machines
  }
}
\section{Changes in \pkg{pomp} version 0.20-4}{
  \itemize{
    \item A bug in the Euler simulator was repaired.
    \item The Euler simulator no longer calls \code{GetRNGState()} and \code{PutRNGState()}.
    When a native routine is used for the Euler step, it must contain these calls.
    \item A potential bug in the Euler-multinomial simulator was repaired.
  }
}
\section{Changes in \pkg{pomp} version 0.20-2}{
  \itemize{
    \item The lookup_table (see file \file{include/pomp.h}) interface for C codes has been streamlined.
    It is now slightly more straightforward to set up a \code{struct lookup_table}.
    See the examples (\file{include/examples}).
    \item Added \file{logistic.R} to the \file{tests} directory.
    This is a one-variable stochastic differential equation simulated using Euler-Maruyama.
    \item The SIR example (\file{examples/sir.R}) now includes two pomps: one with native routines implementing rprocess, dprocess, and skeleton, the other using R files.
    The comparison may be instructive.
    \item Several bugs fixed.
  }
}
\section{Changes in \pkg{pomp} version 0.20-1}{
  \itemize{
    \item This is a major revision.
    \item The basic structure of the \code{pomp} class has been changed to better support the use of native routines as well as R functions.
    In particular, the \code{rmeasure}, \code{dmeasure}, and \code{skeleton} methods can now be specified either as R functions or as native routines.
    Native routines are specified by name; an optional \code{PACKAGE} argument is provided to narrow the search.
    The search is done not at the time of the call to \code{pomp} but at the time the function is used.
    \item Support has been added for covariates, in the form of lookup tables.
    This has been accomplished through the addition of slots \code{tcovar}, \code{covar}, and \code{covarnames} to the \code{pomp} class.
    \code{rmeasure}, \code{dmeasure}, \code{skeleton}, and \code{init.state} now perform lookup when covariates exist and pass this information to the user-defined function.
    \code{rprocess} and \code{dprocess} receive the full lookup table (\code{tcovar} plus \code{covar}).
    \item \code{rmeasure}, \code{dmeasure}, and \code{skeleton} are now specified by a function (or native routine) that does the computation for ONE state vector and ONE parameter vector.
    This is analogous to the way that \code{init.state} has operated from the beginning.
    It greatly simplifies model specification from the user's point of view without changing the low-level methods interfaces (pomp methods \code{rmeasure}, \code{dmeasure}, \code{skeleton}, \code{init.states}).
    However, \strong{this change is not backward-compatible}: users will need to rewrite code for the \code{rmeasure}, \code{dmeasure}, and \code{skeleton} slots.
    In practice, this should not be difficult.
    It involves removing loops over \code{times} and columns of \code{params} and (if covariates are used in the measurement model or deterministic skeleton) making use of the \code{covars} argument passed to the user function.
    \item The interface to the \code{euler.simulate} and \code{euler.density} facilities have changed to make them consistent with the new pomp class.
    \item The header files have been reorganized: \file{include/euler.h} and \file{include/lookup_table.h} have been replaced with a single \file{include/pomp.h}.
    Significant additions have been made to support the new structure of the pomp class.
    Specifically, new prototypes have been defined for native routines to implement the measurement model, the process model, and the deterministic skeleton.
  }
}
