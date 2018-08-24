##' Simulated annealing with box constraints.
##'
##' A straightforward implementation of simulated annealing with box constraints.
##'
##' @name sannbox
##' @rdname sannbox
##' @importFrom stats rnorm
##'
##' @details
##' The \code{control} argument is a list that can supply any of the following components:
##' \describe{
##' \item{trace}{ Non-negative integer.
##' If positive, tracing information on the progress of the optimization is
##' produced.  Higher values may produce more tracing information.  }
##' \item{fnscale}{ An overall scaling to be applied to the value of
##' \code{fn} during optimization.  If negative, turns the problem into a
##' maximization problem.  Optimization is performed on \code{fn(par)/fnscale}.}
##' \item{parscale}{ A vector of scaling values for the parameters.
##' Optimization is performed on \code{par/parscale} and these should be
##' comparable in the sense that a unit change in any element produces about a
##' unit change in the scaled value.  }
##' \item{maxit}{ The total number of function evaluations: there is no
##' other stopping criterion.  Defaults to \code{10000}.  }
##' \item{temp}{ starting temperature for the cooling
##' schedule.  Defaults to \code{1}.  }
##' \item{tmax}{ number of function evaluations at each temperature.
##' Defaults to \code{10}.  }
##' \item{candidate.dist}{ function to randomly select a new candidate
##' parameter vector.  This should be a function with three arguments, the
##' first being the current parameter vector, the second the temperature, and
##' the third the parameter scaling.  By default, \code{candidate.dist} is
##' \preformatted{function(par,temp,scale)
##'                 rnorm(n=length(par),mean=par,sd=scale*temp).} }
##' \item{sched}{ cooling schedule.  A function of a three arguments giving the
##' temperature as a function of iteration number and the control parameters
##' \code{temp} and \code{tmax}.
##' By default, \code{sched} is
##' \preformatted{function(k,temp,tmax) temp/log(((k-1)\%/\%tmax)*tmax+exp(1)).}
##' Alternatively, one can supply a numeric vector of temperatures.
##' This must be of length at least \code{maxit}. }
##' \item{lower,upper}{ optional
##' numeric vectors.  These describe the lower and upper box constraints,
##' respectively.  Each can be specified either as a single scalar (common to
##' all parameters) or as a vector of the same length as \code{par}.  By
##' default, \code{lower=-Inf} and \code{upper=Inf}, i.e., there are no
##' constraints.} }
##'
##' @param par Initial values for the parameters to be optimized over.
##' @param fn A function to be minimized, with first argument the vector of
##' parameters over which minimization is to take place.  It should return a
##' scalar result.
##' @param control A named list of control parameters.  See \sQuote{Details}.
##' @param \dots ignored.
##' @return \code{sannbox} returns a list with components:
##' \describe{
##' \item{counts}{ two-element integer vector.  The first number gives the
##' number of calls made to \code{fn}.  The second number is provided for
##' compatibility with \code{\link{optim}} and will always be NA.  }
##' \item{convergence}{provided for compatibility with \code{\link{optim}};
##' will always be 0.} \item{final.params}{last tried value of \code{par}.}
##' \item{final.value}{value of \code{fn} corresponding to
##' \code{final.params}.}
##' \item{par}{best tried value of \code{par}.}
##' \item{value}{value of \code{fn} corresponding to \code{par}.} }
##'
##' @author Daniel Reuman, Aaron A. King
##'
##' @seealso \code{\link{traj.match}}, \code{\link{probe.match}}.
##'
##' @keywords optimize
NULL

######################
## Simulated annealing minimizer with box constraints
##
## By default, the annealing schedule is
## temp / log(((k-1) %/% tmax)*tmax + exp(1)), where
## the parameters of this schedule can be changed via
## the 'control' argument, and k ranges from 0 to 'maxit'
##
## modified from code originally written by
## Daniel Reuman, Imperial College London

##' @rdname sannbox
##' @export
sannbox <- function (par, fn, control = list(), ...) {

  ep <- "sannbox"
  big <- 1e35  ## a very large number

  npar <- length(par)
  neval <- 0

  control.default <- list(
    maxit=10000,
    temp=1,
    tmax=10,
    sched=NULL,
    candidate.dist=NULL,
    fnscale=1,
    parscale=1,
    lower=-Inf,
    upper=Inf,
    trace=0
  )
  control.default[names(control)] <- control
  control <- control.default
  if (is.null(control$lower)) control$lower <- -Inf
  if (is.null(control$upper)) control$upper <- Inf

  if (is.null(control$sched))           # default cooling schedule
    control$sched <- function (k, temp, tmax) temp/log(((k-1)%/%tmax)*tmax+exp(1))

  if (is.function(control$sched))
    temps <- vapply(
      seq_len(control$maxit),
      FUN=control$sched,
      FUN.VALUE=numeric(1),
      temp=control$temp,
      tmax=control$tmax
    )
  else if (is.numeric(control$sched)) {
    temps <- control$sched
    if (length(temps)<control$maxit)
      pStop(ep,"insufficiently many temperatures supplied in ",
        sQuote("control$sched"))
  }

  if (is.null(control$candidate.dist))
    candidate.dist <- function (par, temp, scale)
      rnorm(n=npar,mean=par,sd=scale*temp)
  else if (is.function(control$candidate.dist)) {
    candidate.dist <- control$candidate.dist
    if (!all(c("par","temp","scale") %in% names(formals(candidate.dist))))
      pStop(ep,sQuote("candidate.dist")," must be a function of prototype ",
        sQuote("candidate.dist(par, temp, scale, ...)"),".")
  } else
    pStop(ep,sQuote("control$candidate.dist")," must be a function.")

  if (length(control$lower)<npar)
    control$lower <- rep(control$lower,npar)
  if (length(control$upper)<npar)
    control$upper <- rep(control$upper,npar)

  ## initialization for the algorithm
  thetabest <- thetacurrent <- par
  ycurrent <- fn(thetacurrent,...)/control$fnscale
  if (!is.finite(ycurrent)) ycurrent <- big
  ybest <- ycurrent
  neval <- 1

  if (control$trace>0)
    cat("initial evaluation: ",ycurrent,"\n")
  if (control$trace>2)
    cat("initial parameters: ",thetacurrent,"\n")

  ## main loop
  for (k in seq_len(control$maxit)) {
    ## get a candidate thetacand
    thetacand <- candidate.dist(thetacurrent,temps[k],control$parscale)
    ## enforce box constraints
    thetacand <- ifelse(
      thetacand<control$lower,
      control$lower,
      thetacand
    )
    thetacand <- ifelse(
      thetacand>control$upper,
      control$upper,
      thetacand
    )
    ycand <- fn(thetacand,...)/control$fnscale
    if (!is.finite(ycand)) ycand <- big
    neval <- neval+1

    ## see if you have a new best.params
    if (ycand<ybest) {
      ybest <- ycand
      thetabest <- thetacand
    }

    accept <- runif(1)<exp((ycurrent-ycand)/temps[k])
    if (accept) { # simulated annealing step
      thetacurrent <- thetacand
      ycurrent <- ycand
    }

    if (control$trace>1)
      cat("iter ",k," val=",ycurrent,", accept=",accept,"\n")
    if (control$trace>3)
      cat("proposed params: ",thetacand,"\n")
    if (control$trace>2)
      cat("current params: ",thetacurrent,"\n")

  }

  if (control$trace>0)
    cat("best val=",ybest,"\n")

  names(thetacurrent) <- names(thetabest) <- names(par)

  list(
    counts=c(neval,NA),
    convergence=0,
    final.params=thetacurrent,
    final.value=ycurrent,
    par=thetabest,
    value=ybest
  )
}
