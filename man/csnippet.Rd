\name{Csnippet}
\title{C code snippets for accelerating computations}
\alias{Csnippet}
\alias{Csnippet-class}
\description{
  For including snippets of C code in \code{pomp} objects.
}
\usage{
  Csnippet(text)
}
\arguments{
  \item{text}{
    character; a snippet of C code.
  }
}
\value{
  An object of class \code{Csnippet}.
}
\section{Using C snippets to accelerate computations}{
  From version 0.50, \pkg{pomp} provides a facility whereby users can define their model's components using inline C code.
  Furnishing one or more \code{Csnippet}s as arguments to the \code{pomp} constructor causes them to be written to a C file stored in the \R session's temporary directory, which is then compiled (via \code{\link[=SHLIB]{R CMD SHLIB}}) into a dynamically loadable shared object file.
  This is then loaded as needed.

  \strong{Note to Windows and Mac users:}
  By default, your \R installation may not support \code{\link[=SHLIB]{R CMD SHLIB}}.
  The package website contains installation instructions that show how to enable this powerful feature of \R.

  In writing a \code{Csnippet} one must bear in mind both the \emph{goal} of the snippet, i.e., what computation it is intended to perform, and the \emph{context} in which it will be executed.
  Details of both of these are given below in the form of rules governing the use of \code{Csnippet}s.
  Illustrative examples are given in the tutorials on the \href{http://kingaa.github.io/pomp}{package website}.
}
\section{General rules for writing C snippets}{
  \enumerate{
    \item C snippets must be valid C.  They will embedded verbatim in a template file which will then be compiled by a call to \code{R CMD SHLIB}.  If the resulting file does not compile, an error message wil be generated.  No attempt is made by \pkg{pomp} to interpret this message.  Typically, compilation errors are due to either invalid C syntax or undeclared variables.
    \item State variables, parameters, observables, and covariates must be left undeclared in the snippet.  State variables and parameters must be declared in either the \code{statenames} or \code{paramnames} argument to \code{pomp}, as appropriate. Compiler errors that complain about undeclared state variables or parameters are usually due to failure to include these parameters in the appropriate vector of names.
    \item A C snippet can declare local variables.  Be careful not to use names that match those of state variables, observables, or parameters.  The latter must never be declared within a snippet.
    \item Names of observables are determined by their names in the data.  They must be referred to in measurement model snippets (rmeasure or dmeasure) by those names.
    \item If the \code{pomp} object contains a table of covariates (see \code{\link{pomp}}), then the variables in the covariate table will be available, by their names, in the context within which the snippet is executed.
    \item \R variables with names containing dots (\sQuote{.}) are replaced in the C codes by variable names in which all dots have been replaced by underscores (\sQuote{_}).
    \item The header \file{R.h}, provided with \R, will be included in the generated C file, making all of the \href{http://cran.r-project.org/doc/manuals/r-release/R-exts.html#The-R-API}{\R C API} available for use in the snippet.
    \item The header \file{pomp.h}, provided with \pkg{pomp}, will also be included, making all of the \pkg{pomp} C API available for use in every snippet.  Do \preformatted{file.show(system.file("include/pomp.h",package="pomp"))} to view this header file.
    \item Snippets of C code passed to the \code{globals} argument of \code{pomp} will be included at the head of the generated C file.  This can be used to declare global variables, define useful functions, and include arbitrary header files.
  }
}
\section{Rules for writing \code{rmeasure} snippets}{
  \enumerate{
    \item The goal of such a snippet is to fill the observables with random values drawn from the measurement model distribution.  Accordingly, each observable should be assigned a new value.
    \item In addition to the states, parameters, covariates (if any), and observables, the variable \code{t}, containing the time of the observation, will be defined in the context in which the snippet is executed.
  }
}
\section{Rules for writing \code{dmeasure} snippets}{
  \enumerate{
    \item The goal of such a snippet is to fill the \code{lik} variable with the likelihood of the data given the state.
    Alternatively, if \code{give_log=1}, \code{lik} should be filled with the log likelihood.
    \item In addition to the states, parameters, covariates (if any), and observables, the variable \code{t}, containing the time of the observation, will be defined in the context in which the snippet is executed.
  }
}
\section{Rules for writing \code{euler.sim} and \code{discrete.time.sim} snippets}{
  \enumerate{
    \item The goal of such a snippet is to replace the state variables with their new random values at the end of the time interval.  Accordingly, each state variable should be over-written with its new value.
    \item In addition to the states, parameters, covariates (if any), and observables, the variables \code{t} and \code{dt}, containing respectively the time at the beginning of the Euler step and the Euler step's duration, will be defined in the context in which the snippet is executed.
  }
}
\section{Rules for writing \code{skeleton} snippets}{
  \enumerate{
    \item For each state variable, there is a corresponding component of the deterministic skeleton.  The goal of such a snippet is to compute all the components.
    \item When the skeleton is a map, the component corresponding to state variable \code{x} is named \code{Dx} and is the new value of \code{x} after one iteration of the map.
    \item When the skeleton is a vectorfield, the component corresponding to state variable \code{x} is named \code{Dx} and is the value of $dx/dt$.
    \item As with the other C snippets, all states, parameters and covariates, as well as the current time, \code{t}, will be defined in the context within which the snippet is executed.
  }
}
\section{Rules for writing \code{rprior} snippets}{
  \enumerate{
    \item The goal of such a snippet is the replacement of parameters with values drawn from the prior distribution.
    \item Within the context in which the snippet will be evaluated, only the parameters will be defined.
    \item Hyperparameters can be included in the ordinary parameter list.
    Obviously, hyperparameters should not be replaced with random draws.
  }
}
\section{Rules for writing \code{dprior} snippets}{
  \enumerate{
    \item The goal of such a snippet is computation of the prior probability density, or the log of same, at a given point in parameter space.
    This scalar value should be returned in the variable \code{lik}.
    When \code{give_log=1}, the user should return the log of the prior probability density.
    \item Within the context in which the snippet will be evaluated, only the parameters and \code{give_log} will be defined.
    \item Hyperparameters can be included in the ordinary parameter list.
    Obviously, hyperparameters should not be replaced with random draws.
  }
}
\section{Rules for writing parameter transformation snippets}{
  \enumerate{
    \item The parameter transformation mapping a parameter vector from the scale used by the model codes to another scale is specified using the \code{toEstimationScale} argument whilst the transformation mapping a parameter vector from the alternative scale to that on which the model is defined is specified with the \code{fromEstimationScale} argument.
    \item The goal of these snippets is the computation of the values of the transformed parameters.  The value of transformed parameter \code{x} should be assigned to variable \code{Tx}.
    \item Time-, state-, and covariate-dependent transformations are not allowed.  Therefore, neither the time, nor any state variables, nor any of the covariates will be available in the context within which a parameter transformation snippet is executed.
  }
}
\section{Rules for writing \code{initializer} snippets}{
  \enumerate{
    \item The goal of a this snippet is the construction of a state vector, i.e., the setting of the dynamical states at time \eqn{t_0}{\code{t0}}.
    \item In addition to the parameters and covariates (if any), the variable \code{t}, containing the time of the observation, will be defined in the context in which the snippet is executed.
    \item \strong{NB:} The \code{statenames} argument plays a particularly important role when the initializer is specified using a \code{Csnippet}.
    In particular, every state variable must be named in \code{statenames}.
    \strong{Failure to follow this rule will result in undefined behavior.}
  }
}
\section{Viewing the generated C code}{
  It can be useful to view the C code generated by calling \code{pomp} with one or more \code{Csnippet} arguments.
  To do so, set \code{options(verbose=TRUE)} before calling \code{pomp}.
  A message giving the name of the generated C file (in the session temporary directory) will be printed.
}
\author{Aaron A. King}
\seealso{
  \code{\link{pomp}}, \code{\link{plugins}}, and the tutorials on the \href{http://kingaa.github.io/pomp}{package website}.
}
\keyword{programming}
\keyword{interface}
