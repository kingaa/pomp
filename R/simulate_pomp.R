##' Simulations of a partially-observed Markov process
##'
##' \code{simulate} generates simulations of the state and measurement
##' processes.
##'
##' Simulation of the state process and of the measurement process are each
##' accomplished by a single call to the user-supplied \code{rprocess} and
##' \code{rmeasure} functions, respectively.  This makes it possible for the
##' user to write highly optimized code for these potentially expensive
##' computations.
##'
##' @name simulate
##' @docType methods
##' @rdname simulate
##' @include workhorses.R pomp_class.R
##'
##' @param object An object of class \sQuote{pomp}.
##' @param nsim The number of simulations to perform.  Note that the number of
##' replicates will be \code{nsim} times \code{ncol(params)}.
##' @param seed optional; if set, the pseudorandom number generator (RNG) will
##' be initialized with \code{seed}.  the random seed to use.  The RNG will be
##' restored to its original state afterward.
##' @param params either a named numeric vector or a numeric matrix with
##' rownames.  The parameters to use in simulating the model.  If \code{params}
##' is not given, then the contents of the \code{params} slot of \code{object}
##' will be used, if they exist.
##' @param states Do we want the state trajectories?
##' @param obs Do we want data-frames of the simulated observations?
##' @param times,t0 \code{times} specifies the times at which simulated
##' observations will be made.  \code{t0} specifies the start time (the time at
##' which the initial conditions hold).  The default for \code{times} is is
##' \code{times=time(object,t0=FALSE)} and \code{t0=timezero(object)},
##' respectively.
##' @param as.data.frame,include.data logical; if \code{as.data.frame=TRUE},
##' the results are returned as a data-frame.  A factor variable, \sQuote{sim},
##' distinguishes one simulation from another.  If, in addition,
##' \code{include.data=TRUE}, the original data are included as an additional
##' \sQuote{simulation}.  If \code{as.data.frame=FALSE}, \code{include.data} is
##' ignored.
##' @param rprocess,rmeasure,obsnames Specifications of the model latent-state
##' simulator (\code{rprocess}), the model measurement simulator
##' (\code{rmeasure}), and names of observables (\code{obsnames}).  See
##' \code{\link{pomp}} for details.
##' @param \dots Additional arguments are passed to \code{\link{pomp}},
##' allowing one to supply new or modify existing model characteristics or
##' components.
##' @param verbose logical; setting \code{verbose = TRUE} will print more
##' information to the console.
##'
##' @return
##' If \code{states=FALSE} and \code{obs=FALSE} (the default), a list
##' of \code{nsim} \sQuote{pomp} objects is returned.  Each has a simulated
##' data set, together with the parameters used (in slot \code{params}) and the
##' state trajectories also (in slot \code{states}).  If \code{times} is
##' specified, then the simulated observations will be at times \code{times}.
##'
##' If \code{nsim=1}, then a single \sQuote{pomp} object is returned (and not a
##' singleton list).
##'
##' If \code{states=TRUE} and \code{obs=FALSE}, simulated state trajectories
##' are returned as a rank-3 array with dimensions \code{nvar} x
##' \code{(ncol(params)*nsim)} x \code{ntimes}.  Here, \code{nvar} is the
##' number of state variables and \code{ntimes} the length of the argument
##' \code{times}.  The measurement process is not simulated in this case.
##'
##' If \code{states=FALSE} and \code{obs=TRUE}, simulated observations are
##' returned as a rank-3 array with dimensions \code{nobs} x
##' \code{(ncol(params)*nsim)} x \code{ntimes}.  Here, \code{nobs} is the
##' number of observables.
##'
##' If both \code{states=TRUE} and \code{obs=TRUE}, then a named list is
##' returned.  It contains the state trajectories and simulated observations as
##' above.
##'
##' @author Aaron A. King
NULL

setGeneric(
  "simulate",
  function (object, nsim=1, seed=NULL, ...)
    standardGeneric("simulate")
)

##' @name simulate-pomp
##' @aliases simulate simulate,pomp-method
##' @rdname simulate
setMethod(
  "simulate",
  signature=signature(object="pomp"),
  definition=function (object, nsim = 1, seed = NULL, params,
    states = FALSE, obs = FALSE, times, t0, as.data.frame = FALSE,
    include.data = FALSE, ..., verbose = getOption("verbose", FALSE)) {

    simulate.internal(
      object=object,
      nsim=nsim,
      seed=seed,
      params=params,
      states=states,
      obs=obs,
      times=times,
      t0=t0,
      as.data.frame=as.data.frame,
      include.data=include.data,
      ...,
      verbose=verbose
    )

  }
)

##' @name simulate-missing
##' @aliases simulate,missing-method
##' @rdname simulate
setMethod(
  "simulate",
  signature=signature(object="missing"),
  definition=function (object, nsim = 1, seed = NULL, params,
    states = FALSE, obs = FALSE, times, t0, as.data.frame = FALSE,
    include.data = FALSE, rprocess, rmeasure, obsnames, ...,
    verbose = getOption("verbose", FALSE)) {

    ep <- paste0("in ",sQuote("simulate"),": ")

    if (missing(times))
      stop(ep,sQuote("times")," is a required argument.",call.=FALSE)
    if (missing(t0))
      stop(ep,sQuote("t0")," is a required argument.",call.=FALSE)

    object <- construct_pomp(data=NULL,times=times,t0=t0,
      rprocess=rprocess,rmeasure=rmeasure,obsnames=obsnames,...,
      verbose=verbose)

    simulate.internal(
      object=object,
      nsim=nsim,
      seed=seed,
      params=params,
      states=states,
      obs=obs,
      times=times,
      t0=t0,
      as.data.frame=as.data.frame,
      include.data=include.data,
      verbose=verbose
    )
  }
)

simulate.internal <- function (object, nsim = 1L, seed = NULL, params,
  states = FALSE, obs = FALSE, times, t0, as.data.frame = FALSE,
  include.data = FALSE, .getnativesymbolinfo = TRUE, verbose, ...) {

  ep <- paste0("in ",sQuote("simulate"),": ")

  object <- pomp(object,...,verbose=verbose)

  obs <- as.logical(obs)
  states <- as.logical(states)
  as.data.frame <- as.logical(as.data.frame)
  include.data <- as.logical(include.data)

  if (length(nsim)!=1 || !is.numeric(nsim) || !is.finite(nsim) || nsim < 1)
    stop(ep,sQuote("nsim")," must be a positive integer.",call.=FALSE)
  nsim <- as.integer(nsim)

  ## set the random seed (be very careful about this)
  seed <- as.integer(seed)
  if (length(seed)>0) {
    if (!exists(".Random.seed",envir=.GlobalEnv)) set.seed(NULL)
    save.seed <- get(".Random.seed",envir=.GlobalEnv)
    set.seed(seed)
  }

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)
  if (is.null(params)) params <- numeric(0)

  params <- as.matrix(params)
  storage.mode(params) <- "double"

  if (missing(times))
    times <- time(object,t0=FALSE)
  else
    times <- as.numeric(times)

  if (missing(t0))
    t0 <- timezero(object)
  else
    t0 <- as.numeric(t0)

  if (!obs && !states)
    object <- as(object,"pomp")

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  retval <- tryCatch(
    .Call(
      simulation_computations,
      object,
      params,
      times,
      t0,
      nsim,
      obs,
      states,
      .getnativesymbolinfo
    ),
    error = function (e) {
      stop(ep,conditionMessage(e),call.=FALSE)
    }
  )
  .getnativesymbolinfo <- FALSE

  ## restore the RNG state
  if (length(seed) > 0) {
    assign('.Random.seed',save.seed,envir=.GlobalEnv)
  }

  if (as.data.frame) {
    if (obs && states) {
      dm <- dim(retval$obs)
      nsim <- dm[2L]
      nm <- rownames(retval$obs)
      dim(retval$obs) <- c(dm[1L],prod(dm[-1L]))
      rownames(retval$obs) <- nm
      dm <- dim(retval$states)
      nm <- rownames(retval$states)
      dim(retval$states) <- c(dm[1L],prod(dm[-1L]))
      rownames(retval$states) <- nm
      retval <- cbind(
        as.data.frame(t(retval$obs)),
        as.data.frame(t(retval$states))
      )
      retval$sim <- seq_len(nsim)
      retval$time <- rep(times,each=nsim)
    } else if (obs || states) {
      dm <- dim(retval)
      nsim <- dm[2L]
      nm <- rownames(retval)
      dim(retval) <- c(dm[1L],prod(dm[-1L]))
      rownames(retval) <- nm
      retval <- as.data.frame(t(retval))
      retval$sim <- seq_len(nsim)
      retval$time <- rep(times,each=nsim)
    } else {
      nsim <- length(retval)
      if (nsim > 1) {
        retval <- lapply(
          seq_len(nsim),
          function (k) {
            x <- as.data.frame(retval[[k]])
            x$sim <- as.integer(k)
            x
          }
        )
        retval <- do.call(rbind,retval)
      } else {
        retval <- as.data.frame(retval)
        retval$sim <- 1L
      }
    }

    if (include.data) {
      od <- as.data.frame(object)
      od$sim <- 0L
      tryCatch(
        {
          retval <- merge(od,retval,all=TRUE)
        },
        error = function (e) {
          stop(ep,"error in merging actual and simulated data.\n",
            "Check names of data, covariates, and states for conflicts.\n",
            sQuote("merge")," error message: ",conditionMessage(e),call.=FALSE)
        }
      )
    }

    retval$sim <- ordered(retval$sim)
    if (include.data) levels(retval$sim)[1L] <- "data"
    retval <- retval[order(retval$sim,retval[[object@timename]]),]

  }

  retval
}
