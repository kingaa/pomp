\name{Process model plugins}
\title{Plug-ins for state-process models}
\alias{process model plugins}
\alias{plugins}
\alias{onestep.sim}
\alias{euler.sim}
\alias{discrete.time.sim}
\alias{onestep.dens}
\alias{gillespie.sim}
\description{
  The functions documented here can be used to specify the \code{rprocess} and \code{dprocess} slots for a \code{pomp} model.
  There are options for discrete- and continuous-time Markov processes.
}
\usage{
onestep.sim(step.fun, PACKAGE)
euler.sim(step.fun, delta.t, PACKAGE)
discrete.time.sim(step.fun, delta.t = 1, PACKAGE)
gillespie.sim(rate.fun, v, d, PACKAGE)
onestep.dens(dens.fun, PACKAGE)
}
\arguments{
  \item{step.fun}{
    This can be either an \R function, the name of a compiled, dynamically loaded native function containing the model simulator, or a \code{\link{Csnippet}}.
    The latter is the preferred option.
    It should be written to take a single Euler step from a single point in state space.

    For an explanation and examples on the use of \code{\link{Csnippet}} to write fast simulators, see the tutorials on the \href{http://kingaa.github.io/pomp}{package website}.

    If it is an \R function, it should have prototype \preformatted{step.fun(x,t,params,delta.t,...).}
    Here, \code{x} is a named numeric vector containing the value of the state process at time \code{t},
    \code{params} is a named numeric vector containing parameters,
    and \code{delta.t} is the length of the Euler time-step.

    If \code{step.fun} is the name of a native function, it must be of type \preformatted{pomp_onestep_sim} as defined in the header file \file{pomp.h}, which is included with the \pkg{pomp} package.
    Do \preformatted{file.show(system.file("include/pomp.h",package="pomp"))} to view this header file.
    For details on how to write such codes, see Details.
  }
  \item{rate.fun}{
    This can be either an \R function, a \code{\link{Csnippet}}, or the name of a compiled, dynamically loaded native function that computes the transition rates.
    If it is an \R function, it should be of the form \code{rate.fun(j,x,t,params,...)}.
    Here, \code{j} is the number of the event,
    \code{x} is a named numeric vector containing the value of the state process at time \code{t} and 
    \code{params} is a named numeric vector containing parameters.

    For examples on the use of \code{\link{Csnippet}} to write fast simulators easily, see tutorials on the \href{http://kingaa.github.io/pomp}{package website}.

    If \code{rate.fun} is a native function, it must be of type \preformatted{pomp_ssa_rate_fn} as defined in the header \file{pomp.h}, which is included with the package.
    For details on how to write such codes, see Details.
  }
  \item{v, d}{
    Matrices that specify the continuous-time Markov process in terms of its elementary events.
    Each should have dimensions \code{nvar} x \code{nevent}, where \code{nvar} is the number of state variables and \code{nevent} is the number of elementary events.
    \code{v} describes the changes that occur in each elementary event:
    it will usually comprise the values 1, -1, and 0 according to whether a state variable is incremented, decremented, or unchanged in an elementary event.
    \code{d} is a binary matrix that describes the dependencies of elementary event rates on state variables:
    \code{d[i,j]} will have value 1 if event rate \code{j} must be updated as a result of a change in state variable \code{i} and 0 otherwise
    }
  \item{dens.fun}{
    This can be either an R function, a \code{\link{Csnippet}}, or a compiled, dynamically loaded native function containing the model transition log probability density function.

    If it is an R function, it should be of the form \preformatted{dens.fun(x1,x2,t1,t2,params,...).}
    Here, \code{x1} and \code{x2} are named numeric vectors containing the values of the state process at times \code{t1} and \code{t2},
    \code{params} is a named numeric vector containing parameters.

    If \code{dens.fun} is the name of a native function, it should be of type \preformatted{pomp_onestep_pdf} as defined in the header \file{pomp.h}, which is included with the \pkg{pomp} package.
    This function should return the log likelihood of a transition from \code{x1} at time \code{t1} to \code{x2} at time \code{t2}, assuming that no intervening transitions have occurred.
    For details on how to write such codes, see Details.
  }
  \item{delta.t}{
    Size of time-steps.
    See below.
  }
  \item{PACKAGE}{
    an optional argument that specifies to which dynamically loaded library we restrict the search for the native routines.
    If this is \dQuote{base}, we search in the R executable itself.
    This argument is ignored if \code{step.fun}, \code{rate.fn}, or \code{dens.fun} is provided as an \R function or a \code{\link{Csnippet}}.
  }
}
\details{
  \code{onestep.sim} is the appropriate choice when it is possible to simulate the change in state from one time to another, regardless of how large the interval between them is.
  To use \code{onestep.sim}, you must write a function \code{step.fun} that will advance the state process from one arbitrary time to another.
  \code{euler.sim} is appropriate when one cannot do this but can compute the change in state via a sequence of smaller steps.
  This is desirable, for example, if one is simulating a continuous time process but is willing to approximate it using an Euler approach.
  \code{discrete.time.sim} is appropriate when the process evolves in discrete time.

  To use \code{euler.sim} or \code{discrete.time.sim}, you must write a function \code{step.fun} that will take a single Euler step, of size at most \code{delta.t}.
  The functions \code{euler.sim} and \code{discrete.time.sim} will create simulators that take as many steps as needed to get from one time to another.
  See below for information on how \code{euler.sim} chooses the actual step size it uses.

  \code{gillespie.sim} allows exact simulation of a continuous-time, discrete-state Markov process using Gillespie's algorithm.
  This is an \dQuote{event-driven} approach: correspondingly, to use \code{gillespie.sim}, you must write a function \code{rate.fun} that computes the rates of each elementary kind of event and specify two matrices (\code{d,v}) that describe, respectively, the dependencies of each rate and the consequences of each event.

  \code{onestep.dens} will generate a suitable \code{dprocess} function when one can compute the likelihood of a given state transition simply by knowing the states at two times under the assumption that the state has not changed between the times.
  This is typically possible, for instance, when the \code{rprocess} function is implemented using \code{onestep.sim}, \code{euler.sim}, or \code{discrete.time.sim}.
  [NB: currently, there are no high-level algorithms in \pkg{pomp} that use \code{dprocess}.
  This function is provided for completeness only, with an eye toward future development.]

  If \code{step.fun} is written as an \R function, it must have at least the arguments \code{x}, \code{t}, \code{params}, \code{delta.t}, and \code{\dots}.
  On a call to this function, \code{x} will be a named vector of state variables, \code{t} a scalar time, and \code{params} a named vector of parameters.
  The length of the Euler step will be \code{delta.t}.
  If the argument \code{covars} is included and a covariate table has been included in the \code{pomp} object, then on a call to this function, \code{covars} will be filled with the values, at time \code{t}, of the covariates.
  This is accomplished via interpolation of the user-supplied covariate table.
  Additional arguments may be given: these will be filled by the correspondingly-named elements in the \code{userdata} slot of the \code{pomp} object (see \code{\link{pomp}}).

  If \code{step.fun} is written in a native language, it must be a function of type \preformatted{pomp_onestep_sim} as specified in the header \file{pomp.h} included with the package.
  Execute \preformatted{file.show(system.file("include/pomp.h",package="pomp.h"))} to view this file.

  If \code{rate.fun} is written as an \R function, it must have at least the arguments \code{j}, \code{x}, \code{t}, \code{params}, and \code{\dots}.
  Here, \code{j} is the an integer that indicates for which of the elementary events the current rate is desired.
  \code{x} is a named vector containing the value of the state process at time \code{t}, and
  \code{params} is a named vector containing parameters.
  If the argument \code{covars} is included and a covariate table has been included in the \code{pomp} object, then on a call to this function, \code{covars} will be filled with the values, at time \code{t}, of the covariates.
  This is accomplished via interpolation of the covariate table.
  If \code{rate.fun} is a native function, it must be of type \preformatted{pomp_ssa_rate_fn} as defined in the header \file{pomp.h}; see above for instructions on how to view this file.
  
  In writing \code{dens.fun}, you must assume that no state transitions have occurred between \code{t1} and \code{t2}.
  If \code{dens.fun} is written as an \R function, it must have at least the arguments \code{x1}, \code{x2}, \code{t1}, \code{t2}, \code{params}, and \code{\dots}.
  On a call to this function, \code{x1} and \code{x2} will be named vectors of state variables at times \code{t1} and \code{t2}, respectively.
  The named vector \code{params} contains the parameters.
  If the argument \code{covars} is included and a covariate table has been included in the \code{pomp} object, then on a call to this function, \code{covars} will be filled with the values, at time \code{t1}, of the covariates.
  If the argument \code{covars} is included and a covariate table has been included in the \code{pomp} object, then on a call to this function, \code{covars} will be filled with the values, at time \code{t1}, of the covariates.
  This is accomplished via interpolation of the covariate table.
  As above, any additional arguments will be filled by the correspondingly-named elements in the \code{userdata} slot of the \code{pomp} object (see \code{\link{pomp}}).
  If \code{dens.fun} is written in a native language, it must be a function of type \preformatted{pomp_onestep_pdf} as defined in the header \file{pomp.h} included with the package; see above for instructions on how to view this file.
}
\section{Size of time step}{
  All the simulator plug-ins documented here work by taking discrete time steps.
  They differ as to how this is done.
  
  \code{onestep.sim} takes a single step to go from any given time \eqn{t_1}{t1} to any other time \eqn{t_2}{t2} (\eqn{t_1 < t_2}{t1 < t2}).
  Thus, this plugin is only useful if a closed-form solution to the process exists.

  To go from \eqn{t_1}{t1} to \eqn{t_2}{t2}, \code{euler.sim} takes \eqn{n} steps of equal size, where \eqn{n = ceiling((t_2-t_1)/delta.t)}{n = ceiling((t2-t1)/delta.t)}.

  \code{discrete.time.sim} assumes that the process evolves in discrete time, where the interval between successive times is \code{delta.t}.
  Thus, to go from \eqn{t_1}{t1} to \eqn{t_2}{t2}, \code{discrete.time.sim} takes \eqn{n} steps of size exactly \eqn{delta.t}, where \eqn{n = floor((t_2-t_1)/delta.t)}{n = floor((t2-t1)/delta.t)}.
}
\value{
  \code{onestep.sim}, \code{euler.sim}, \code{discrete.time.sim}, and \code{gillespie.sim} each return functions for use as the arguments \code{rprocess} argument in \code{\link{pomp}}.

  \code{onestep.dens} returns a function for use as the argument \code{dprocess} in \code{\link{pomp}}.
}
\author{Aaron A. King}
\seealso{
  \code{\link{pomp}} and the tutorials on the \href{http://kingaa.github.io/pomp}{package website}.
}
\keyword{models}
