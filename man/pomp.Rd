\name{pomp constructor}
\title{Constructor of the basic pomp object}
\alias{pomp}
\alias{pomp-class}
\alias{pomp constructor}
\alias{accumulator variables}
\alias{process model plug-ins}
\alias{plug-ins}
\alias{onestep.sim}
\alias{euler.sim}
\alias{discrete.time.sim}
\alias{onestep.dens}
\alias{gillespie.sim}
\alias{gillespie.hl.sim}
\alias{Csnippet}
\alias{Csnippet-class}
\alias{coerce,Csnippet,character-method}
\alias{vectorfield}
\alias{map}
\alias{parameter_trans}
\alias{parameter_trans,ANY,ANY-method}
\alias{parameter_trans,missing,missing-method}
\alias{parameter_trans,ANY,missing-method}
\alias{parameter_trans,missing,ANY-method}
\alias{parameter_trans,Csnippet,Csnippet-method}
\alias{parameter_trans,character,character-method}
\alias{parameter_trans,function,function-method}
\alias{parameter_trans,NULL,NULL-method}
\alias{parameter_trans,pomp_fun,pomp_fun-method}
\alias{covariate_table}
\alias{covariate_table,character-method}
\alias{covariate_table,numeric-method}
\alias{covariate_table,missing-method}
\description{
  This function constructs a \sQuote{pomp} object,
  encoding a partially-observed Markov process model together with a uni- or multi-variate time series.
  As such, it is central to all the package's functionality.
  One implements the model by specifying some or all of its \emph{basic components}.
  These include:
  \describe{
    \item{rprocess,}{the simulator of the unobserved Markov state process;}
    \item{dprocess,}{the evaluator of the probability density function for transitions of the unobserved Markov state process;}
    \item{rmeasure,}{the simulator of the observed process, conditional on the unobserved state;}
    \item{dmeasure,}{the evaluator of the measurement model probability density function;}
    \item{rinit,}{which samples from the distribution of the state process at the zero-time;}
    \item{rprior,}{which samples from a prior probability distribution on the parameters;}
    \item{dprior}{which evaluates the prior probability density function;}
    \item{skeleton}{which computes the deterministic skeleton of the unobserved state process.}
  }
  The basic structure and its rationale are described in the \emph{Journal of Statistical Software} paper, an updated version of which is to be found on the \href{https://kingaa.github.io/pomp/}{package website}.
}
\usage{
pomp(data, times, t0, \dots, rinit, rprocess, dprocess, rmeasure, dmeasure,
     skeleton, rprior, dprior, partrans, params, covar, zeronames,
     obsnames, statenames, paramnames, covarnames, PACKAGE, globals,
     cdir, cfile, shlib.args)
}
\arguments{
  \item{data, times}{
    required; the time series data and times at which observations are made.
    \code{data} should be given as a data-frame and \code{times} must indicate the column of observation times by name or index.
    \code{times} must be numeric and strictly increasing.
    Internally, \code{data} will be internally coerced to an array with storage-mode \code{double}.

    In addition, a \sQuote{pomp} object can be supplied in the \code{data} argument.
    In this case, the call to \code{pomp} will add element to, or replace elements of, the supplied \sQuote{pomp} object.
  }
  \item{t0}{
    The zero-time, at which the stochastic dynamical system is to be initialized.
    This must be no later than the time of the first observation, i.e., \code{t0 <= times[1]}.
    This argument is required whenever \code{data} is a data-frame.
  }
  \item{rinit}{
    optional; draws from the distribution of initial values of the unobserved Markov state process.
    Specifically, given a vector of parameters, \code{params} and an initial time, \code{t0}, the rinit determines the state vector at time \code{t0}.
    See below under \dQuote{The State-Process Initializer} for details.
  }
  \item{rprocess, dprocess}{
    optional; specification of the simulator and probability density evaluation function of the unobserved state process.
    See below under \dQuote{The Unobserved Markov State-Process Model} for details.

    \strong{Note:} it is not typically necessary (or even feasible) to define \code{dprocess}.
    In fact, no current \pkg{pomp} inference algorithm makes use of \code{dprocess}.
    This functionality is provided only to support future algorithm development.
  }
  \item{rmeasure, dmeasure}{
    optional; specifications of the measurement model.
    See below under \dQuote{The Measurement Model} for details.
  }
  \item{skeleton}{
    optional; the deterministic skeleton of the unobserved state process.
    See below under \dQuote{The Deterministic Skeleton} for details.
  }
  \item{rprior, dprior}{
    optional; specification of the prior distribution on parameters.
    See below under \dQuote{Specifying a Prior} for details.
  }
  \item{partrans}{
    optional parameter transformations.
    Many algorithms for parameter estimation search an unconstrained space of parameters.
    When working with such an algorithm and a model for which the parameters are constrained, it can be useful to transform parameters.
    One should supply the \code{partrans} argument via a call to \code{parameter_trans(toEst,fromEst,\dots,log,logit,barycentric)}.
    See below under \dQuote{Parameter Transformations} for more details.
  }
  \item{params}{
    optional; named numeric vector of parameters.
    This will be coerced internally to storage mode \code{double}.
  }
  \item{covar}{
    optional covariate table, constructed using \code{covariate_table}.

    If a covariate table is supplied, then the value of each of the covariates is interpolated as needed.
    The resulting interpolated values are made available to the appropriate basic components.
    Note that \code{covar} will be coerced internally to storage mode \code{double}.
    See below under \dQuote{Covariates} for more details.
  }
  \item{obsnames, statenames, paramnames, covarnames}{
    optional character vectors specifying the names of observables, state variables, parameters, and covariates, respectively.
    These are used only in the event that one or more of the basic components are defined using C snippets or native routines.
    It is usually unnecessary to specify \code{obsnames} or \code{covarnames}, as these will by default be read from \code{data} and \code{covars}, respectively.
  }
  \item{zeronames}{
    optional character vector specifying the names of accumulator variables (see below under \dQuote{Accumulator Variables}).
  }
  \item{PACKAGE}{
    optional string giving the name of the dynamically loaded library in which any native routines are to be found.
    This is only useful if one or more of the model components has been specified using a precompiled dynamically loaded library;
    it is not used for any component specified using C snippets.
  }
  \item{globals}{
    optional character;
    C code that will be included in the source for (and therefore hard-coded into) the shared-object library created when the call to \code{pomp} uses C snippets.
    If no C snippets are used, \code{globals} has no effect.
  }
  \item{cdir, cfile, shlib.args}{
    optional character variables.
    \code{cdir} specifies the name of the directory within which C snippet code will be compiled.
    By default, this is in a temporary directory specific to the running instance of \R.
    \code{cfile} gives the name of the file (in directory \code{cdir}) into which C snippet codes will be written.
    By default, a random filename is used.
    The \code{shlib.args} can be used to pass command-line arguments to the \code{R CMD SHLIB} call that will compile the C snippets.
  }
  \item{\dots}{
    Any additional arguments given to \code{pomp} will be made available to each of the user-defined functions.
    To prevent errors due to misspelling, a warning is issued if any such arguments are detected.
  }
}
\value{
  The \code{pomp} constructor function returns an object, call it \code{P}, of class \sQuote{pomp}.
  \code{P} contains, in addition to the data, any elements of the model that have been specified as arguments to the \code{pomp} constructor function.
  One can add or modify elements of \code{P} by means of further calls to \code{pomp}, using \code{P} as the first argument in such calls.
}
\section{Important note}{
  \strong{
    It is not typically necessary (or even feasible) to define all of the basic components for any given purpose.
    Each \pkg{pomp} algorithm makes use of only a subset of these components.
    Any algorithm requiring a component that is not present will generate an error letting you know that you have not provided a needed component.
  }
}
\section{Using C snippets to accelerate computations}{
  \pkg{pomp} provides a facility whereby users can define their model's components using inline C code.
  Furnishing one or more C snippets as arguments to the \code{pomp} constructor causes them to be written to a C file stored in the \R session's temporary directory, which is then compiled (via \code{\link[=SHLIB]{R CMD SHLIB}}) into a dynamically loadable shared object file.
  This is then loaded as needed.

  \strong{Note to Windows and Mac users:}
  By default, your \R installation may not support \code{\link[=SHLIB]{R CMD SHLIB}}.
  The \href{https://kingaa.github.io/pomp/install.html}{package website contains installation instructions} that explain how to enable this powerful feature of \R.
}
\section{General rules for writing C snippets}{
  In writing a C snippet one must bear in mind both the \emph{goal} of the snippet, i.e., what computation it is intended to perform, and the \emph{context} in which it will be executed.
  These are explained here in the form of general rules.
  Additional specific rules apply according to the function of the particular C snippet.
  Illustrative examples are given in the tutorials on the \href{https://kingaa.github.io/pomp/}{package website}.
  \enumerate{
    \item C snippets must be valid C.
    They will embedded verbatim in a template file which will then be compiled by a call to \code{R CMD SHLIB}.
    If the resulting file does not compile, an error message will be generated.
    No attempt is made by \pkg{pomp} to interpret this message.
    Typically, compilation errors are due to either invalid C syntax or undeclared variables.
    \item State variables, parameters, observables, and covariates must be left undeclared within the snippet.
    State variables and parameters are declared via the \code{statenames} or \code{paramnames} arguments to \code{pomp}, respectively.
    Compiler errors that complain about undeclared state variables or parameters are usually due to failure to declare these in \code{statenames} or \code{paramnames}, as appropriate.
    \item A C snippet can declare local variables.
    Be careful not to use names that match those of state variables, observables, or parameters.
    The latter must never be declared within a C snippet.
    \item Names of observables are determined by their names in the data.
    They must be referred to in measurement model C snippets (\code{rmeasure} and \code{dmeasure}) by those names.
    \item If the \sQuote{pomp} object contains a table of covariates (see above), then the variables in the covariate table will be available, by their names, in the context within which the C snippet is executed.
    \item Because the dot \sQuote{.} has syntactic meaning in C, \R variables with names containing dots (\sQuote{.}) are replaced in the C codes by variable names in which all dots have been replaced by underscores (\sQuote{_}).
    \item The header \file{R.h}, provided with \R, will be included in the generated C file, making all of the \href{https://CRAN.R-project.org/doc/manuals/r-release/R-exts.html#The-R-API}{\R C API} available for use in the C snippet.
    This makes a great many useful functions available, including all of \R's \href{https://CRAN.R-project.org/doc/manuals/r-release/R-exts.html#Distribution-functions}{statistical distribution functions}.
    \item The header \href{https://github.com/kingaa/pomp/blob/master/inst/include/pomp.h}{\file{pomp.h}}, provided with \pkg{pomp}, will also be included, making all of the \pkg{pomp} C API available for use in every C snippet.
    Do \preformatted{file.show(system.file("include/pomp.h",package="pomp"))} to view this header file.
    \item Snippets of C code passed to the \code{globals} argument of \code{pomp} will be included at the head of the generated C file.
    This can be used to declare global variables, define useful functions, and include arbitrary header files.
  }
}
\section{The Unobserved Markov State-Process Model}{
  Specification of process-model codes \code{rprocess} and/or \code{dprocess} is facilitated by \pkg{pomp}'s so-called plug-ins, which allow one to easily specify the most common kinds of process model.

  \subsection{Discrete-time processes}{
    If the state process evolves in discrete time, specify \code{rprocess} using the \code{discrete.time.sim} plug-in.
    Specifically, provide \preformatted{
      rprocess = discrete.time.sim(step.fun = f, delta.t)}
    to \code{pomp}, where \code{f} is a C snippet or \R function that simulates one step of the state process.
    The former is the preferred option, due to its much greater computational efficiency.
    The goal of such a C snippet is to replace the state variables with their new random values at the end of the time interval.
    Accordingly, each state variable should be over-written with its new value.
    In addition to the states, parameters, covariates (if any), and observables, the variables \code{t} and \code{dt}, containing respectively the time at the beginning of the step and the step's duration, will be defined in the context in which the C snippet is executed.
    See above under \dQuote{General rules for C snippet writing} for more details.
    Examples are to be found in the tutorials on the \href{https://kingaa.github.io/pomp/}{package website}.

    If \code{f} is given as an \R function, it should have prototype \preformatted{
      f(x, t, params, delta.t, ...)}
    When \code{f} is called, \code{x} will be a named numeric vector containing the value of the state process at time \code{t},
    \code{params} will be a named numeric vector containing parameters,
    and \code{delta.t} will be the time-step.
    It should return a named vector of the same length, and with the same set of names, as \code{x}, representing a draw from the distribution of the state process at time \code{t+delta.t}, conditional on its having value \code{x} at time \code{t}.
  }

  \subsection{Continuous-time processes}{
    If the state process evolves in continuous time, but you can use an Euler approximation, specify \code{rprocess} using the \code{euler.sim} plug-in.
    Furnish \preformatted{
      rprocess = euler.sim(step.fun = f, delta.t)}
    to \code{pomp} in this case.
    As before, \code{f} can be provided either as a C snippet or as an \R function, the former resulting in much quicker computations.
    The form of \code{f} will be the same as above (in the discrete-time case).

    If you have a procedure that allows you, given the value of the state process at any time, to simulate it at an arbitrary time in the future, use the \code{onestep.sim} plug-in.
    To do so, furnish \preformatted{
      rprocess = onestep.sim(step.fun = f)}
    to \code{pomp}.
    Again, \code{f} can be provided either as a C snippet or as an \R function, the former resulting in much quicker computations.
    The form of \code{f} should be as above (in the discrete-time or Euler cases).

    If you desire exact simulation of certain continuous-time Markov chains, an implementation of Gillespie's algorithm (Gillespie 1977) is available, via the \code{gillespie.sim} and \code{gillespie.hl.sim} plug-ins.
    The former allows for the rate function to be provided as an \R function or a single C snippet, while the latter provides a means of specifying the elementary events via a list of C snippets.

    A high-level interface to the simulator is provided by \code{gillespie.hl.sim}.
    To use it, supply \preformatted{
      rprocess = gillespie.hl.sim(..., .pre = "", .post = "", hmax = Inf)}
    to \code{pomp}.
    Each argument in \code{...} corresponds to a single elementary event and should be a list containing two elements.
    The first should be a string or C snippet;
    the second should be a named integer vector.
    The string contains C code that assigns to the variable \code{rate} the corresponding elementary event rate.
    The variable \code{rate} will exist in the context of the code, as will the covariates (as named in \code{covar}) and the time \code{t}.
    The state variables and parameters will also be available provided that the \code{statenames} and \code{paramnames} arguments are provided to \code{pomp}.
    The named integer vector contains the changes to the state variables corresponding to the elementary event.
    There should be named value for each of the state variables returned by the rinit.
    Typically, these changes will be integers but floating point values are permissible.
    The arguments \code{.pre} and \code{.post} can be used to provide C code that will run either before or after the code in the event arguments.
    These hooks could be useful for avoiding duplication of code that performs calculations needed to obtain several of the different event rates.
    Here's how a simple birth-death model might be specified: \preformatted{
      gillespie.hl.sim(
        birth=list("rate = b*N;",c(N=1)),
        death=list("rate = m*N;",c(N=-1)))}
    In the above, the state variable \code{N} represents the population size and parameters \code{b}, \code{m} are the birth and death rates, respectively.

    To use the lower-level \code{gillespie.sim} interface, furnish \preformatted{
      rprocess = gillespie.sim(rate.fun = f, v, hmax = Inf)}
    to \code{pomp}, where \code{f} gives the rates of the elementary events.
    Here, \code{f} may be an \R function of the form \preformatted{
      f(j, x, t, params, ...)}
    When \code{f} is called,
    the integer \code{j} will be the number of the elementary event (corresponding to the column the matrix \code{v}, see below),
    \code{x} will be a named numeric vector containing the value of the state process at time \code{t} and
    \code{params} is a named numeric vector containing parameters.
    \code{f} should return a single numerical value, representing the rate of that elementary event at that point in state space and time.

    Matrix \code{v} specifies the continuous-time Markov process in terms of its elementary events.
    It should have dimensions \code{nvar} x \code{nevent}, where \code{nvar} is the number of state variables and \code{nevent} is the number of elementary events.
    \code{v} describes the changes that occur in each elementary event:
    it will usually comprise the values 1, -1, and 0 according to whether a state variable is incremented, decremented, or unchanged in an elementary event.
    The rows of \code{v} may be unnamed or named.
    If the rows are unnamed, they are assumed to be in the same order as the vector of state variables returned by the rinit.
    If the rows are named, the names of the state variables returned by the rinit will be matched to the rows of \code{v} to ensure a correct mapping.
    If any of the row names of \code{v} cannot be found among the state variables or if any row names of \code{v} are duplicated, an error will occur.

    It is also possible to provide a C snippet via the \code{rate.fn} argument to \code{gillespie.sim}.
    Such a snippet should assign the correct value to a \code{rate} variable depending on the value of \code{j}.
    The same variables will be available as for the C code provided to \code{gillespie.hl.sim}.
    This lower-level interface may be preferable if it is easier to write code that calculates the correct rate based on \code{j} rather than to write a snippet for each possible value of \code{j}.
    For example, if the number of possible values of \code{j} is large and the rates vary according to a few simple rules, the lower-level interface may provide the easier way of specifying the model.

    When the process is non-autonomous (i.e., the event rates depend explicitly on time), it can be useful to set \code{hmax} to the maximum step that will be taken.
    By default, the elementary event rates will be recomputed at least once per observation interval.
  }
  \subsection{Size of time step}{
    The simulator plug-ins \code{discrete.time.sim}, \code{euler.sim}, and \code{onestep.sim} all work by taking discrete time steps.
    They differ as to how this is done.
    Specifically,
    \enumerate{
      \item \code{onestep.sim} takes a single step to go from any given time \code{t1} to any later time \code{t2} (\code{t1 < t2}).
      Thus, this plug-in is designed for use in situations where a closed-form solution to the process exists.
      \item To go from \code{t1} to \code{t2}, \code{euler.sim} takes \code{n} steps of equal size, where \preformatted{
	n = ceiling((t2-t1)/delta.t).}
      \item \code{discrete.time.sim} assumes that the process evolves in discrete time, where the interval between successive times is \code{delta.t}.
      Thus, to go from \code{t1} to \code{t2}, \code{discrete.time.sim} takes \code{n} steps of size exactly \code{delta.t}, where \preformatted{
	n = floor((t2-t1)/delta.t).}
    }
  }
  \subsection{Specifying \code{dprocess}}{
    If you have a procedure that allows you to compute the probability density of an arbitrary transition from state \eqn{x_1}{x1} at time \eqn{t_1}{t1} to state \eqn{x_2}{x2} at time \eqn{t_2>t_1}{t2}, assuming that the state remains unchanged between \eqn{t_1}{t1} and \eqn{t_2}{t2}, then you can furnish \preformatted{
      dprocess = f}
    to \code{pomp}, where \code{f} is a C snippet or \R function that computes the \emph{log} probability density.
    Using a C snippet is preferred, due to its much greater computational efficiency.
    The goal of such a C snippet is to fill the variable \code{loglik} with the log probability density.
    In the context of such a C snippet, the parameters, and covariates will be defined, as will the times \code{t_1} and \code{t_2}.
    The state variables at time \code{t_1} will have their usual name (see \code{statenames}) with a \dQuote{\code{_1}} appended.
    Likewise, the state variables at time \code{t_2} will have a \dQuote{\code{_2}} appended.
    See above under \dQuote{General rules for C snippet writing} for more details.
    Examples are to be found in the tutorials on the \href{https://kingaa.github.io/pomp/}{package website}.

    If \code{f} is given as an \R function, it should have prototype \preformatted{
      f(x1, x2, t1, t2, params, ...)}
    When \code{f} is called, \code{x1} will be a named numeric vector containing the value of the state process at time \code{t1},
    \code{x2} will be a named numeric vector containing the value of the state process at time \code{t2},
    and \code{params} will be a named numeric vector containing parameters.
    \code{f} should return the \emph{log} likelihood of a transition from \code{x1} at time \code{t1} to \code{x2} at time \code{t2}, assuming that no intervening transitions have occurred.
  }

  To see examples, consult the tutorials on the \href{https://kingaa.github.io/pomp/}{package website}.

}
\section{The Measurement Model}{
  The measurement model is the link between the data and the unobserved state process.
  It can be specified either by using one or both of the \code{rprocess} and \code{dprocess} arguments.

  The best way to specify the measurement model is by giving C snippets for \code{rmeasure} and \code{dmeasure}.
  In writing an \code{rmeasure} C snippet, bear in mind that:
  \enumerate{
    \item The goal of such a snippet is to fill the observables with random values drawn from the measurement model distribution.
    Accordingly, each observable should be assigned a new value.
    \item In addition to the states, parameters, covariates (if any), and observables, the variable \code{t}, containing the time of the observation, will be defined in the context in which the snippet is executed.
  }
  General rules for writing C snippets are provided above.
  The tutorials on the \href{https://kingaa.github.io/pomp/}{package website} give examples as well.

  It is also possible, though far less efficient, to specify \code{rmeasure} using an \R function.
  In this case, specify the measurement model simulator by furnishing \preformatted{
    rmeasure = f}
  to \code{pomp}, where \code{f} is an \R function with prototype \preformatted{
    f(x, t, params, \dots)}
  It can also take any additional arguments if these are passed along with it in the call to \code{pomp}.
  When \code{f} is called,
  \itemize{
    \item \code{x} will be a named numeric vector of length \code{nvar}, the number of state variables.
    \item \code{t} will be a scalar quantity, the time at which the measurement is made.
    \item \code{params} will be a named numeric vector of length \code{npar}, the number of parameters.
  }
  \code{f} must return a named numeric vector of length \code{nobs}, the number of observable variables.

  In writing a \code{dmeasure} C snippet, observe that:
  \enumerate{
    \item In addition to the states, parameters, covariates (if any), and observables, the variable \code{t}, containing the time of the observation,
    and the Boolean variable \code{give_log} will be defined in the context in which the snippet is executed.
    \item The goal of such a snippet is to set the value of the \code{lik} variable to the likelihood of the data given the state.
    Alternatively, if \code{give_log == 1}, \code{lik} should be set to the log likelihood.
  }

  If \code{dmeasure} is to be provided instead as an \R function, this is accomplished by supplying \preformatted{
    dmeasure = f}
  to \code{pomp}, where \code{f} has prototype \preformatted{
    f(y, x, t, params, log, \dots)}
  Again, it can take additional arguments that are passed with it in the call to \code{pomp}.
  When \code{f} is called,
  \itemize{
    \item \code{y} will be a named numeric vector of length \code{nobs} containing values of the observed variables;
    \item \code{x} will be a named numeric vector of length \code{nvar} containing state variables;
    \item \code{params} will be a named numeric vector of length \code{npar} containing parameters;
    \item \code{t} will be a scalar, the corresponding observation time.
  }
  \code{f} must return a single numeric value, the probability density of \code{y} given \code{x} at time \code{t}.
  If \code{log == TRUE}, then \code{f} should return instead the log of the probability density.
  \strong{Note: it is a common error to fail to account for both \code{log = TRUE} and \code{log = FALSE} when writing the \code{dmeasure} C snippet or function.}
}
\section{The Deterministic Skeleton}{
  The skeleton is a dynamical system that expresses the central tendency of the unobserved Markov state process.
  As such, it is not uniquely defined, but can be both interesting in itself and useful in practice.
  In \pkg{pomp}, the skeleton is used by \code{\link{trajectory}} and \code{\link{traj.match}}.

  If the state process is a discrete-time stochastic process, then the skeleton is a discrete-time map.
  To specify it, provide \preformatted{
    skeleton = map(f, delta.t)}
  to \code{pomp}, where \code{f} implements the map and \code{delta.t} is the size of the timestep covered at one map iteration.

  If the state process is a continuous-time stochastic process, then the skeleton is a vectorfield (i.e., a system of ordinary differential equations).
  To specify it, supply \preformatted{
    skeleton = vectorfield(f)}
  to \code{pomp}, where \code{f} implements the vectorfield, i.e., the right-hand-size of the differential equations.

  In either case, \code{f} can be furnished either as a C snippet (the preferred choice), or an \R function.
  In writing a \code{skeleton} C snippet, be aware that:
  \enumerate{
    \item For each state variable, there is a corresponding component of the deterministic skeleton.
    The goal of such a snippet is to compute all the components.
    \item When the skeleton is a map, the component corresponding to state variable \code{x} is named \code{Dx} and is the new value of \code{x} after one iteration of the map.
    \item When the skeleton is a vectorfield, the component corresponding to state variable \code{x} is named \code{Dx} and is the value of \eqn{dx/dt}.
    \item As with the other C snippets, all states, parameters and covariates, as well as the current time, \code{t}, will be defined in the context within which the snippet is executed.
  }
  The tutorials on the \href{https://kingaa.github.io/pomp/}{package website} give some examples.

  If \code{f} is an \R function, it must be of prototype \preformatted{
    f(x, t, params, \dots)}
  where, as usual,
  \itemize{
    \item \code{x} is a numeric vector (length \code{nvar}) containing the coordinates of a point in state space at which evaluation of the skeleton is desired.
    \item \code{t} is a scalar value giving the time at which evaluation of the skeleton is desired.
    \item \code{params} is a numeric vector (length \code{npar}) holding the parameters.
  }
  As with the other basic components, \code{f} may take additional arguments, provided these are passed along with it in the call to \code{pomp}.
  The function \code{f} must return a numeric vector of the same length as \code{x}, which contains the value of the map or vectorfield at the required point and time.
}
\section{The State-Process Initializer}{
  To fully specify the unobserved Markov state process, one must give its distribution at the zero-time (\code{t0}).
  By default, \code{pomp} assumes that this initial distribution is concentrated on a single point.
  In particular, any parameters in \code{params}, the names of which end in \dQuote{\code{.0}}, are assumed to be initial values of states.
  When the state process is initialized, these are simply copied over as initial conditions.
  The names of the resulting state variables are obtained by dropping the \dQuote{\code{.0}} suffix.

  One can override this default behavior by furnishing a value for the \code{rinit} argument of \code{pomp}.
  As usual, this can be provided either as a C snippet or as an \R function.
  In the former case, bear in mind that:
  \enumerate{
    \item The goal of a this snippet is the construction of a state vector, i.e., the setting of the dynamical states at time \eqn{t_0}{\code{t0}}.
    \item In addition to the parameters and covariates (if any), the variable \code{t}, containing the zero-time, will be defined in the context in which the snippet is executed.
    \item \strong{NB:} The \code{statenames} argument plays a particularly important role when the rinit is specified using a C snippet.
    In particular, every state variable must be named in \code{statenames}.
    \strong{Failure to follow this rule will result in undefined behavior.}
  }

  If an \R function is to be used, pass \preformatted{
    rinit = f}
  to \code{pomp}, where \code{f} is a function with prototype \preformatted{
    f(params, t0, \dots)}
  When \code{f} is called,
  \itemize{
    \item \code{params} will be a named numeric vector of parameters.
    \item \code{t0} will be the time at which initial conditions are desired.
  }
  As usual, \code{f} may take additional arguments, provided these are passed along with it in the call to \code{pomp}.
  \code{f} must return a named numeric vector of initial states.
  It is of course important that the names of the states match the expectations of the other basic components.

  Note that the state-process rinit can be either deterministic (the default) or stochastic.
  In the latter case, it samples from the distribution of the state process at the zero-time, \code{t0}.
}
\section{Specifying a Prior}{
  A prior distribution on parameters is specified by means of the \code{rprior} and/or \code{dprior} arguments to \code{pomp}.
  As with the other basic model components, it is preferable to specify these using C snippets.
  In writing a C snippet for the prior sampler (\code{rprior}), keep in mind that:
  \enumerate{
    \item Within the context in which the snippet will be evaluated, only the parameters will be defined.
    \item The goal of such a snippet is the replacement of parameters with values drawn from the prior distribution.
    \item Hyperparameters can be included in the ordinary parameter list.
    Obviously, hyperparameters should not be replaced with random draws.
  }
  In writing a C snippet for the prior density function (\code{dprior}), observe that:
  \enumerate{
    \item Within the context in which the snippet will be evaluated, only the parameters and \code{give_log} will be defined.
    \item The goal of such a snippet is computation of the prior probability density, or the log of same, at a given point in parameter space.
    This scalar value should be returned in the variable \code{lik}.
    When \code{give_log == 1}, \code{lik} should contain the log of the prior probability density.
    \item Hyperparameters can be included in the ordinary parameter list.
  }

  Alternatively, one can furnish \R functions for one or both of these arguments.
  In this case, \code{rprior} must be a function of prototype \preformatted{
    f(params, \dots)}
  that makes a draw from the prior distribution given \code{params} and returns a named vector of the same length and with the same set of names, as \code{params}.
  The \code{dprior} function must be of prototype \preformatted{
    f(params, log = FALSE, \dots).}
  Its role is to evaluate the prior probability density (or log density if \code{log == TRUE}) and return that single scalar value.
}
\section{Covariates}{
  If the \sQuote{pomp} object contains covariates (specified via the \code{covar} argument; see above), then interpolated values of the covariates will be available to each of the model components whenever it is called.
  In particular, variables with names as they appear in the \code{covar} covariate table will be available to any C snippet.
  When a basic component is defined using an \R function, that function will be called with an extra argument, \code{covars}, which will be a named numeric vector containing the interpolated values from the covariate table.

  An exception to this rule is the prior (\code{rprior} and \code{dprior}): covariate-dependent priors are not allowed.
  Nor are parameter transformations permitted to depend upon covariates.
}
\section{Parameter Transformations}{
  When parameter transformations are desired, they can be integrated into the \sQuote{pomp} object via the \code{partrans} arguments using the \code{parameter_trans(toEst,fromEst,\dots,log,logit,barycentric)} function.
  As with the basic model components, these should ordinarily be specified using C snippets.
  When doing so, note that:
  \enumerate{
    \item The parameter transformation mapping a parameter vector from the scale used by the model codes to another scale, and the inverse transformation, are specified via a call to \code{parameter_trans(toEst,fromEst)}.
    \item The goal of these snippets is the transformation of the parameters from the natural scale to the estimation scale, and vice-versa.
    If \code{p} is the name of a variable on the natural scale, its value on the estimation scale is \code{T_p}.
    Thus the \code{toEst} snippet computes \code{T_p} given \code{p} whilst the \code{fromEst} snippet computes \code{p} given \code{T_p}.
    \item Time-, state-, and covariate-dependent transformations are not allowed.
    Therefore, neither the time, nor any state variables, nor any of the covariates will be available in the context within which a parameter transformation snippet is executed.
  }

  These transformations can also be specified using \R functions with arguments \code{params} and \code{\dots}.
  In this case, \code{toEst} should transform parameters from the scale that the basic components use internally to the scale used in estimation.
  \code{fromEst} should be the inverse of \code{toEst}.

  Note that it is the user's responsibility to make sure that the transformations are mutually inverse.
  If \code{obj} is the constructed \sQuote{pomp} object, and \code{coef(obj)} is non-empty, a simple check of this property is \preformatted{
    x <- coef(obj, transform = TRUE)
    obj1 <- obj
    coef(obj1, transform = TRUE) <- x
    identical(coef(obj), coef(obj1))
    identical(coef(obj1, transform=TRUE), x)}

  By default, both functions are the identity transformation.

  One can use the \code{log} and \code{logit} arguments of \code{parameter_trans} to name variables that should be log-transformed or logit-transformed, respectively.
  The \code{barycentric} argument can name sets of parameters that should be log-barycentric transformed.
}
\section{Accumulator Variables}{
  In formulating models, one sometimes wishes to define a state variable that will accumulate some quantity over the interval between successive observations.
  \pkg{pomp} provides a facility to make such features more convenient.
  Specifically, variables named in the \code{pomp}'s \code{zeronames} argument will be set to zero immediately following each observation.
  See \code{\link{sir}} and the tutorials on the \href{https://kingaa.github.io/pomp/}{package website} for examples.
}
\section{Viewing generated C code}{
  It can be useful to view the C code generated by calling \code{pomp} with one or more C snippet arguments.
  You can set \code{cdir} and \code{cfile} to control where this code is written.
  Alternatively, set \code{options(verbose=TRUE)} before calling \code{pomp}.
  This will cause a message giving the name of the generated C file (in the session temporary directory) to be printed.
}
\examples{

## pomp encoding a stochastic Ricker model with a covariate:

pomp(data = data.frame(t = 1:100, y = NA),
     times = "t", t0 = 0,
     covar = covariate_table(
       data.frame(t=0:100,K=seq(from=50,to=200,length=101)),
       times = "t"),
     rprocess = discrete.time.sim(Csnippet("double e = rnorm(0,sigma);
                                            n = n*exp(r*(1-n/K)+e);"), delta.t = 1),
     skeleton = map(Csnippet("Dn = n*exp(r*(1-n/K));"), delta.t = 1),
     rmeasure = Csnippet("y = rnbinom_mu(theta,n);"),
     dmeasure = Csnippet("lik = dpois(y,n,give_log);"),
     rprior = Csnippet("r = rgamma(1,1); sigma = rgamma(1,1);"),
     dprior = Csnippet("lik = dgamma(r,1,1,1) + dgamma(sigma,1,1,1);
                        if (!give_log) lik = exp(lik);"),
     rinit = Csnippet("n = n_0;"),
     partrans=parameter_trans(log=c("r","sigma")),
     paramnames = c("n_0", "r", "sigma","theta"),
     statenames = "n") -> rick

## fill it with simulated data:

rick <- simulate(rick, params = c(r=log(17), sigma = 0.1, theta=1, n_0 = 50))
plot(rick)

\dontrun{
    pompExample()
    demo(package="pomp")
}
}
\section{Warning}{
  Some error checking is done by \code{pomp}, but complete error checking for arbitrary models is impossible.
  If the user-specified functions do not conform to the above specifications, then the results may be invalid.
  In particular, if both \code{rmeasure} and \code{dmeasure} are specified, the user should verify that these two functions correspond to the same probability distribution.
  If \code{skeleton} is specified, the user is responsible for verifying that it corresponds to a deterministic skeleton of the model.
}
\author{Aaron A. King}
\references{
  A. A. King, D. Nguyen, and E. L. Ionides (2016)
  Statistical Inference for Partially Observed Markov Processes via the \R Package \pkg{pomp}.
  Journal of Statistical Software 69(12): 1--43.

  D. T. Gillespie (1977)
  Exact stochastic simulation of coupled chemical reactions.
  Journal of Physical Chemistry 81:2340--2361.
}
\seealso{
  \link{pomp methods},
  \link{pomp low-level interface}
}
\keyword{models}
\keyword{interface}
\keyword{programming}
\keyword{ts}
