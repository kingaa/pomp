##' The latent state process simulator
##'
##' Specification of the latent state process simulator, rprocess.
##'
##' @name rprocess specification
##' @rdname rprocess_spec
##' @docType methods
##' @include pomp_fun.R csnippet.R
##' @family implementation information
##' @seealso \code{\link{rprocess}}
##' 
##' @section Discrete-time processes:
##' If the state process evolves in discrete time, specify \code{rprocess} using the \code{discrete_time} plug-in.
##' Specifically, provide
##' \preformatted{
##'     rprocess = discrete_time(step.fun = f, delta.t),
##' }
##' where \code{f} is a C snippet or \R function that simulates one step of the state process.
##' The former is the preferred option, due to its much greater computational efficiency.
##' The goal of such a C snippet is to replace the state variables with their new random values at the end of the time interval.
##' Accordingly, each state variable should be over-written with its new value.
###' In addition to the states, parameters, covariates (if any), and observables, the variables \code{t} and \code{dt}, containing respectively the time at the beginning of the step and the step's duration, will be defined in the context in which the C snippet is executed.
##' See \code{\link{Csnippet}} for general rules on writing C snippets.
##' Examples are to be found in the tutorials on the \href{https://kingaa.github.io/pomp/}{package website}.
##'
##' If \code{f} is given as an \R function, its arguments should come from the state variables, parameters, covariates, and time.
##' It may also take the argument \sQuote{\code{delta.t}};
##' when called, the latter will be the timestep.
##' It must also have the argument \sQuote{\code{...}}.
##' It should return a named vector of length equal to the number of state variables, representing a draw from the distribution of the state process at time \code{t+delta.t} conditional on its value at time \code{t}.
##'
##' @section Continuous-time processes:
##' If the state process evolves in continuous time, but you can use an Euler approximation, implement \code{rprocess} using the \code{euler} plug-in.
##' Specify
##' \preformatted{
##'     rprocess = euler(step.fun = f, delta.t)
##' }
##' in this case.
##' As before, \code{f} can be provided either as a C snippet or as an \R function, the former resulting in much quicker computations.
##' The form of \code{f} will be the same as above (in the discrete-time case).
##'
##' If you have a procedure that allows you, given the value of the state process at any time,
##' to simulate it at an arbitrary time in the future, use the \code{onestep} plug-in.
##' To do so, specify
##' \preformatted{
##'     rprocess = onestep(step.fun = f).
##' }
##' Again, \code{f} can be provided either as a C snippet or as an \R function, the former resulting in much quicker computations.
##' The form of \code{f} should be as above (in the discrete-time or Euler cases).
##' @section Size of time step:
##' The simulator plug-ins \code{discrete_time}, \code{euler}, and \code{onestep} all work by taking discrete time steps.
##' They differ as to how this is done.
##' Specifically,
##' \enumerate{
##' \item \code{onestep} takes a single step to go from any given time \code{t1} to any later time \code{t2} (\code{t1 < t2}).
##' Thus, this plug-in is designed for use in situations where a closed-form solution to the process exists.
##' \item To go from \code{t1} to \code{t2}, \code{euler} takes \code{n} steps of equal size, where
##' \preformatted{
##'     n = ceiling((t2-t1)/delta.t).
##' }
##' \item \code{discrete_time} assumes that the process evolves in discrete time, where the interval between successive times is \code{delta.t}.
##' Thus, to go from \code{t1} to \code{t2}, \code{discrete_time} takes \code{n} steps of size exactly \code{delta.t}, where
##' \preformatted{
##'     n = floor((t2-t1)/delta.t).
##' }
##' }
##'
##' @section Exact (event-driven) simulations:
##' If you desire exact simulation of certain continuous-time Markov chains, an implementation of Gillespie's algorithm (Gillespie 1977) is available,
##' via the \code{gillespie} and \code{gillespie_hl} plug-ins.
##' The former allows for the rate function to be provided as an \R function or a single C snippet,
##' while the latter provides a means of specifying the elementary events via a list of C snippets.
##'
##' A high-level interface to the simulator is provided by \code{gillespie_hl}.
##' To use it, supply
##' \preformatted{
##'     rprocess = gillespie_hl(..., .pre = "", .post = "", hmax = Inf)
##' }
##' to \code{pomp}.
##' Each argument in \code{...} corresponds to a single elementary event and should be a list containing two elements.
##' The first should be a string or C snippet;
##' the second should be a named integer vector.
##' The variable \code{rate} will exist in the context of the C snippet, as will the parameter, state variables, covariates, and the time \code{t}.
##' The C snippet should assign to the variable \code{rate} the corresponding elementary event rate.
##'
##' The named integer vector specifies the changes to the state variables corresponding to the elementary event.
##' There should be named value for each of the state variables returned by \code{rinit}.
##' The arguments \code{.pre} and \code{.post} can be used to provide C code that will run respectively before and after the elementary-event snippets.
##' These hooks can be useful for avoiding duplication of code that performs calculations needed to obtain several of the different event rates.
##'
##' Here's how a simple birth-death model might be specified:
##' \preformatted{
##'     gillespie_hl(
##'         birth=list("rate = b*N;",c(N=1)),
##'         death=list("rate = m*N;",c(N=-1))
##'     )
##' }
##' In the above, the state variable \code{N} represents the population size and parameters \code{b}, \code{m} are the birth and death rates, respectively.
##'
##' To use the lower-level \code{gillespie} interface, furnish
##' \preformatted{
##'     rprocess = gillespie(rate.fun = f, v, hmax = Inf)
##' }
##' to \code{pomp}, where \code{f} gives the rates of the elementary events.
##' Here, \code{f} may be an \R function of the form
##' \preformatted{
##'     f(j, x, t, params, ...)
##' }
##' When \code{f} is called,
##' the integer \code{j} will be the number of the elementary event (corresponding to the column the matrix \code{v}, see below),
##' \code{x} will be a named numeric vector containing the value of the state process at time \code{t} and
##' \code{params} is a named numeric vector containing parameters.
##' \code{f} should return a single numerical value, representing the rate of that elementary event at that point in state space and time.
##'
##' Here, the stoichiometric matrix \code{v} specifies the continuous-time Markov process in terms of its elementary events.
##' It should have dimensions \code{nvar} x \code{nevent}, where \code{nvar} is the number of state variables and \code{nevent} is the number of elementary events.
##' \code{v} describes the changes that occur in each elementary event:
##' it will usually comprise the values 1, -1, and 0 according to whether a state variable is incremented, decremented, or unchanged in an elementary event.
##' The rows of \code{v} should have names corresponding to the state variables.
##' If any of the row names of \code{v} cannot be found among the state variables or if any row names of \code{v} are duplicated, an error will occur.
##'
##' It is also possible to provide a C snippet via the \code{rate.fun} argument to \code{gillespie}.
##' Such a snippet should assign the correct value to a \code{rate} variable depending on the value of \code{j}.
##' The same variables will be available as for the C code provided to \code{gillespie_hl}.
##' This lower-level interface may be preferable if it is easier to write code that calculates the correct rate based on \code{j} rather than to write a snippet for each possible value of \code{j}.
##' For example, if the number of possible values of \code{j} is large and the rates vary according to a few simple rules, the lower-level interface may provide the easier way of specifying the model.
##'
##' When the process is non-autonomous (i.e., the event rates depend explicitly on time), it can be useful to set \code{hmax} to the maximum step that will be taken.
##' By default, the elementary event rates will be recomputed at least once per observation interval.
##'
##' @section Default behavior:
##' The default \code{rprocess} is undefined.
##' It will yield missing values (\code{NA}) for all state variables.
##'
##' @param step.fun a C snippet, an R function, or
##' the name of a native routine in a shared-object library.
##' This gives a procedure by which one simulates a single step of the latent state process.
##' @param delta.t positive numerical value; for \code{euler} and \code{discrete_time}, the size of the step to take
##' @param rate.fun a C snippet, an R function, or
##' the name of a native routine in a shared-object library.
##' This gives a procedure by which one computes the event-rate of the elementary events in the continuous-time latent Markov chain.
##' @param v integer matrix; giving the stoichiometry of the continuous-time latent Markov process.
##' It should have dimensions \code{nvar} x \code{nevent}, where \code{nvar} is the number of state variables and \code{nevent} is the number of elementary events.
##' \code{v} describes the changes that occur in each elementary event:
##' it will usually comprise the values 1, -1, and 0 according to whether a state variable is incremented, decremented, or unchanged in an elementary event.
##' The rows of \code{v} may be unnamed or named.
##' If the rows are unnamed, they are assumed to be in the same order as the vector of state variables returned by \code{rinit}.
##' If the rows are named,
##' the names of the state variables returned by \code{rinit} will be matched
##' to the rows of \code{v} to ensure a correct mapping.
##' If any of the row names of \code{v} cannot be found among the state variables or if any row names of \code{v} are duplicated, an error will occur.
##' @param .pre,.post C snippets (see Details)
##' @param \dots individual C snippets corresponding to elementary events
##' @param hmax maximum time step allowed (see below)
##'
##' @inheritSection pomp Note for Windows users
##'
NULL

## also defined in 'pomp_internal.h'
rprocmode <- list(default=0L,onestep=1L,discrete=2L,euler=3L,gillespie=4L)

setClass(
  "rprocPlugin",
  slots=c(
    csnippet='logical',
    slotname='character',
    type='integer',
    step.fn="ANY",
    rate.fn="ANY"
  ),
  prototype=prototype(
    csnippet=FALSE,
    slotname=character(0),
    type=rprocmode$default,
    step.fn=NULL,
    rate.fn=NULL
  )
)

##' @rdname undefined
setMethod(
  "undefined",
  signature=signature(object="rprocPlugin"),
  definition=function (object, ...) {
    undefined(object@step.fn) && undefined(object@rate.fn)
  }
)

setClass(
  "onestepRprocPlugin",
  contains="rprocPlugin",
  prototype=prototype(
    type=rprocmode$onestep
  )
)

setClass(
  "discreteRprocPlugin",
  contains="rprocPlugin",
  slots=c(
    delta.t="numeric"
  ),
  prototype=prototype(
    type=rprocmode$discrete,
    delta.t=1.0
  )
)

setClass(
  "eulerRprocPlugin",
  contains="rprocPlugin",
  slots=c(
    delta.t="numeric"
  ),
  prototype=prototype(
    type=rprocmode$euler,
    delta.t=NA_real_
  )
)

setClass(
  "gillespieRprocPlugin",
  contains="rprocPlugin",
  slots=c(
    hmax="numeric",
    v="matrix"
  ),
  prototype=prototype(
    type=rprocmode$gillespie
  )
)

rproc_plugin <- function (object, step.fn, rate.fn) {
  if (missing(object)) {
    new("rprocPlugin")
  } else {
    if (!missing(step.fn)) object@step.fn <- step.fn
    if (!missing(rate.fn)) object@rate.fn <- rate.fn
    object
  }
}

##' @name onestep
##' @rdname rprocess_spec
##' @export
onestep <- function (step.fun) {
  new("onestepRprocPlugin",
    step.fn=step.fun,
    slotname="step.fun",
    csnippet=is(step.fun,"Csnippet")
  )
}

##' @rdname rprocess_spec
##' @name discrete_time
##' @export
discrete_time <- function (step.fun, delta.t = 1) {
  new("discreteRprocPlugin",
    step.fn=step.fun,
    delta.t=delta.t,
    slotname="step.fun",
    csnippet=is(step.fun,"Csnippet")
  )
}

##' @rdname rprocess_spec
##' @name euler
##' @export
euler <- function (step.fun, delta.t) {
  new("eulerRprocPlugin",
    step.fn=step.fun,
    delta.t=delta.t,
    slotname="step.fun",
    csnippet=is(step.fun,"Csnippet")
  )
}

##' @rdname rprocess_spec
##' @name gillespie
##' @export
gillespie <- function (rate.fun, v, hmax = Inf) {
  ep <- "gillespie"
  if (!is.matrix(v)) pStop(ep,sQuote("v")," must be a matrix.")
  if (anyDuplicated(rownames(v)))
    pStop(ep,"duplicates in row names of ",sQuote("v"),".")
  if (!is.null(colnames(v)) && anyDuplicated(colnames(v)))
    pStop(ep,"duplicates in column names of ",sQuote("v"),".")

  new("gillespieRprocPlugin",
    rate.fn=rate.fun,
    v=v,
    hmax=hmax,
    slotname="rate.fun",
    csnippet=is(rate.fun,"Csnippet")
  )
}

##' @rdname rprocess_spec
##' @name gillespie_hl
##' @export
gillespie_hl <- function (..., .pre = "", .post = "", hmax = Inf) {

  ep <- "gillespie_hl"

  args <- list(...)

  if (!all(vapply(args,inherits,what="list",logical(1L))) ||
      !all(vapply(args,length,integer(1L)) == 2L))
    pStop(ep,"each event should be specified using a length-2 list.")

  codeChunks <- lapply(args,"[[",1)
  stoich <- lapply(args,"[[",2)

  checkCode <- function (x) {
    inh <- inherits(x,what=c("Csnippet", "character"))
    if (!any(inh) || length(x) != 1)
      pStop(ep,"for each event, the first list-element should be a C snippet or string.")
    as(x,"character")
  }

  codeChunks <- lapply(codeChunks,checkCode)

  if (!inherits(.pre,what=c("character","Csnippet")) ||
      !inherits(.post,what=c("character","Csnippet")))
    pStop(ep,sQuote(".pre")," and ",sQuote(".post")," must be C snippets or strings.")

  .pre <- paste(as(.pre,"character"),collapse="\n")
  .post <- paste(as(.post,"character"),collapse="\n")

  if (!all(vapply(stoich,is.numeric,logical(1L))) ||
      any(vapply(stoich,\(x)invalid_names(names(x)),logical(1L))))
    pStop(ep,"for each event, the second list-element should be ",
      "a named numeric vector (without duplicate names).")

  ## Create C snippet of switch statement
  header <- paste0(.pre, "\nswitch (j) {")
  body <- paste0(
    sprintf("case %d:\n{\n%s\n}\nbreak;\n",seq_along(codeChunks),codeChunks),
    collapse="\n"
  )
  footer <- paste0("default:\nerror(\"unrecognized event %d\",j);\nbreak;\n}\n",.post)
  rate.fn <- Csnippet(paste(header, body, footer, sep="\n"))

  ## now put together the stoichiometry matrix
  if (anyDuplicated(names(stoich))) pStop(ep,"duplicated elementary event names.")
  vars <- unique(do.call(c,lapply(stoich,names)))
  v <- array(data=0,dim=c(length(vars),length(stoich)),
    dimnames=list(name=vars,event=names(stoich)))
  for (s in seq_along(stoich))
    v[names(stoich[[s]]),s] <- stoich[[s]]

  new("gillespieRprocPlugin",
    rate.fn=rate.fn,
    v=v,
    hmax=hmax,
    slotname="rate.fun",
    csnippet=TRUE
  )
}
