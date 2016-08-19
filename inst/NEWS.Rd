\name{NEWS}
\title{News for package `pomp'}
\section{Changes in \pkg{pomp} version 1.7.3}{
  \itemize{
    \item Fixed a bug in \code{plot} on \code{pomp} objects with more than 10 variables to plot.
    \item Added the \code{shlib.args} option to \code{pomp}, allowing arbitrary options to be passed to the compiler when C snippets are used.
  }
}  
\section{Changes in \pkg{pomp} version 1.7.2}{
  \itemize{
    \item Modified the demos to use best practices.
  }
}  
\section{Changes in \pkg{pomp} version 1.7.1}{
  \itemize{
    \item Improved the test suite.
  }
}  
\section{Changes in \pkg{pomp} version 1.6.7}{
  \itemize{
    \item Improve the test suite and error messages.
    \item Replace the automatic \dQuote{verbose} setting during compilation under Windows.
    This was dropped in version 1.6.4.3, when it no longer appeared to be necessary.
    Apparently, it is, at least under some Windows OS/R combinations.
  }
}  
\section{Changes in \pkg{pomp} version 1.6.6}{
  \itemize{
    \item Fix \code{probe} bug introduced in version 1.6.5.9.
  }
}  
\section{Changes in \pkg{pomp} version 1.6.5}{
  \itemize{
    \item Several obscure bug-fixes.
    \item Improvement of the test suite and of error messages.
  }
}
\section{Changes in \pkg{pomp} version 1.6.4}{
  \itemize{
    \item Remove the automatic \dQuote{verbose} setting during compilation under Windows.
    This was introduced in version 0.63-6 as a workaround.
    It appears no longer to be necessary.
    \item Improve the \code{show} method for \code{pomp} objects.
    \item Change the \code{print} method for \code{pomp} objects.
    \item Improve messages for plugin errors.
    \item Some internal refactoring.
  }
}
\section{Changes in \pkg{pomp} version 1.6.3}{
  \itemize{
    \item Improve error messages for missing basic components.
  }
}
\section{Changes in \pkg{pomp} version 1.6.2}{
  \itemize{
    \item Fix typos in one of the help-page examples.
    \item Add some references.
  }
}
\section{Changes in \pkg{pomp} version 1.6.1}{
  \itemize{
    \item Major revisions to the package help pages and manual.
    The separate documentation on C snippets and process-model plugins has been merged into the main \code{pomp} help page, which now has comprehensive instructions on building a \code{pomp} object.
    \item Update to the package help page.
  }
}
\section{Changes in \pkg{pomp} version 1.5.10}{
  \itemize{
    \item The matrix and numeric vector options for the \code{data} argument to \code{pomp} are now deprecated and will be removed in a future release.
    In calls to \code{pomp}, \code{data} should be either a data frame or a \code{pomp}-class object.
    \item Many more errors are now trapped and reported with more informative messages.
  }
}
\section{Changes in \pkg{pomp} version 1.5.9}{
  \itemize{
    \item Redundant \code{bootstrap} argument to \code{nlf} has been removed.
    \item More improvements to the testing suite.
  }
}
\section{Changes in \pkg{pomp} version 1.5.8}{
  \itemize{
    \item Major improvements to error and warning handling.
    \item Fix bugs that might arise in rare circumstances in \code{simulate} and \code{eulermultinom} functions.
    \item More improvements to the testing suite.
  }
}
\section{Changes in \pkg{pomp} version 1.5.7}{
  \itemize{
    \item Many improvements to the test suite.
    \item Fixed a bug in \code{pompExample}.
  }
}
\section{Changes in \pkg{pomp} version 1.5.6}{
  \itemize{
    \item Major overhaul of error handling system and error messages.
    \item There is a new plugin for \code{rprocess} called \code{kleap.sim}.
    This implements the so-called \dQuote{K-leap method} of Cai & Xu (2007).
  }
}
\section{Changes in \pkg{pomp} version 1.5.5}{
  \itemize{
    \item Many improvements to the test suite.
  }
}
\section{Changes in \pkg{pomp} version 1.5.4}{
  \itemize{
    \item More improvements to the test suite.
    \item Fixed bugs found as a result of improved testing.
    \item In \code{mif2}, it is now required that \code{Nmif>0}.
  }
}
\section{Changes in \pkg{pomp} version 1.5.3}{
  \itemize{
    \item More improvements to the test suite.
    \item Fixed bugs found as a result of improved testing.
    \item Fixed bug in specification of skeleton in pompExamples.
  }
}
\section{Changes in \pkg{pomp} version 1.5.2}{
  \itemize{
    \item Improved the test suite to obtain better code coverage.
    \item \code{pomp} is no longer an S4 method.
    Its new status as an ordinary function abbreviates and consolidates the code and simplifies the documentation.
    \item Some minor improvements to error messages.
    \item A bug restricting access in C code to \code{userdata} has been fixed.
    \item Fixed bug in \code{pompExample(...,show=TRUE)}.
  }
}
\section{Changes in \pkg{pomp} version 1.5.1}{
  \itemize{
    \item New methods for the Ensemble Kalman Filter (\code{enkf}) and Ensemble Adjustment Kalman Filter (\code{eakf}) have been added to the package.
  }
}
\section{Changes in \pkg{pomp} version 1.4.9}{
  \itemize{
    \item A better interface for specifying a model's deterministic skeleton is provided.
    One specifies \code{skeleton=map(f,delta.t)} for a discrete-time skeleton (a map) and \code{skeleton=vectorfield(f)} for a continuous-time skeleton (a vectorfield).
    The old arguments \code{skeleton.type} and \code{skelmap.delta.t} are deprecated and will be removed in a future release.
    \item The \code{pfilter} help page has been improved.
    Specifically, the discussion of filtering failures is better.
    \item Miscellaneous improvements to the help pages.
  }
}
\section{Changes in \pkg{pomp} version 1.4.8}{
  \itemize{
    \item \code{mif2} no longer computes filter means.
  }
}
\section{Changes in \pkg{pomp} version 1.4.7}{
  \itemize{
    \item The new argument \code{show} of \code{pompExamples} allows one to display the example code instead of executing it.
    \item The minimum version of \R supported is now 3.1.2.
  }
}
\section{Changes in \pkg{pomp} version 1.4.6}{
  \itemize{
    \item \code{init.state} now has the optional argument \code{nsim}.
    Using this, one can request multiple initial state vectors per parameter vector.
    \item \code{pfilter} now uses less memory when it is asked to run at a single point in parameter space.
  }
}
\section{Changes in \pkg{pomp} version 1.4.5}{
  \itemize{
    \item A bug associated with the \code{rw.sd} argument to \code{mif2} on Windows platforms has been fixed.
  }
}
\section{Changes in \pkg{pomp} version 1.4.4}{
  \itemize{
    \item The \code{method="mif2"} option to \code{mif} has been removed.
    Use \code{mif2} instead.
  }
}
\section{Changes in \pkg{pomp} version 1.4.3}{
  \itemize{
    \item The \code{particles} method (rarely if ever used), has been removed to streamline the \code{mif} codes.
  }
}
\section{Changes in \pkg{pomp} version 1.4.2}{
  \itemize{
    \item The \code{method="mif2"} option to \code{mif} is now deprecated.
    It will be removed in a future release.
    Use \code{mif2} instead.
  }
}
\section{Changes in \pkg{pomp} version 1.4.1}{
  \itemize{
    \item Following the publication of the Journal of Statistical Software announcement paper, the \file{CITATION} has been modified to include the appropriate citation.
    Please use \code{citation("pomp")} to learn about how to cite \pkg{pomp} correctly.
  }
}
\section{Changes in \pkg{pomp} version 1.3.4}{
  \itemize{
    \item \code{conv.rec} on the results of a \code{mif2} computation now play nicely with \code{reshape2::melt}.
  }
}
\section{Changes in \pkg{pomp} version 1.3.3}{
  \itemize{
    \item It is now possible to direct \code{pomp} to place \code{Csnippet} codes in a specified location.
    This is accomplished via the new \code{cdir} and \code{cfile} arguments to \code{pomp}.
  }
}
\section{Changes in \pkg{pomp} version 1.3.2}{
  \itemize{
    \item When several concurrent \R processes attempt to simultaneously compile \code{Csnipet} code from copies of a \code{pomp} object on the same machine, a race condition exists.
    This release fixes this by ensuring that \code{Csnippet} codes are written and compiled in a process-specific directory.
  }
}
\section{Changes in \pkg{pomp} version 1.3.1}{
  This release removes several long-deprecated features.
  \itemize{
    \item The \code{seed} argument to \code{pfilter} has been removed.
    Use \code{freeze} to obtain similar functionality.
    \item The redundant \code{pars} argument to \code{mif} has been removed.
    Names of parameters to be estimated are taken from the nonzero entries of \code{rw.sd}.
    Initial-value parameters are designated using the \code{ivps} argument.
    \item The old \code{parameter.transform} and \code{parameter.inv.transform} arguments to \code{pomp} have been removed.
    They were superseded by the arguments \code{fromEstimationScale} and \code{toEstimationScale} in version 0.65-1.
  }
}
\section{Changes in \pkg{pomp} version 1.2.2}{
  \itemize{
    \item Some internal changes to \code{mif2} by C. Breto, in support of panel methods under development.
  }
}
\section{Changes in \pkg{pomp} version 1.2.1}{
  \itemize{
    \item Bug fix: for POMPs with a scalar state variable, the \code{filter.mean}, \code{pred.mean}, and \code{pred.var} methods dropped the singleton dimension of the respective arrays, in conflict with the assumptions of certain \code{plot} methods.
    This resolves issue #3.
  }
}
\section{Changes in \pkg{pomp} version 1.1.4}{
  \itemize{
    \item Documentation improvements.
    Specifically, the details of how step-sizes are chosen in the \code{euler.sim}, \code{discrete.time.sim}, and \code{onestep.sim} plugins has been spelled out explicitly.
    Also, incompleteness in the \code{sannbox} documentation has been corrected.
  }
}
\section{Changes in \pkg{pomp} version 1.1.3}{
  \itemize{
    \item Provide \code{filter.traj} method for \code{pmcmcList} objects.
  }
}
\section{Changes in \pkg{pomp} version 1.1.2}{
  \itemize{
    \item Fixed a bug in \code{mif} and \code{mif2} when \code{cooling.type="hyperbolic"} and \code{cooling.fraction.50=1}.
  }
}
