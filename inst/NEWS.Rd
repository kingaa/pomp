\name{NEWS}
\title{News for package `pomp'}
\section{Changes in \pkg{pomp} version 1.5.1}{
  \itemize{
    \item New methods for the Ensemble Kalman Filter (\code{enkf}) and Ensemble Adjustment Kalman Filter (\code{eakf}) have been added to the package.
    \item Improved the test suite to obtain better code coverage.
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
