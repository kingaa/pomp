##' Simulations of a partially-observed Markov process
##'
##' \code{simulate} generates simulations of the state and measurement
##' processes.
##'
##' @name simulate
##' @docType methods
##' @rdname simulate
##' @include workhorses.R pomp_class.R pomp.R
##' @author Aaron A. King
##' @family elementary algorithms
##' @inheritSection pomp Note for Windows users
##' @inheritParams pomp
##' @param object optional;
##' if present, it should be a data frame or a \sQuote{pomp} object.
##' @param params a named numeric vector or a matrix with rownames
##' containing the parameters at which the simulations are to be performed.
##' @param nsim The number of simulations to perform.
##' Note that the number of replicates will be \code{nsim} times \code{ncol(params)}.
##' @param seed optional integer;
##' if set, the pseudorandom number generator (RNG) will be initialized with \code{seed}.
##' The RNG will be restored to its original state afterward.
##' @param format the format in which to return the results.
##'
##' \code{format = "pomps"} causes the results to be returned as a single \dQuote{pomp} object, if \code{params} is a vector, or a list of \dQuote{pomp} objects, if \code{params} is a matrix with more than one column.
##' Each of these will be identical to \code{object} except in that the latent states and observations will have been replaced by their simulated values.
##'
##' \code{format = "arrays"} causes the results to be returned as a list of two arrays.
##' The \dQuote{states} element will contain the simulated state trajectories in a rank-3 array with dimensions
##' \code{nvar} x \code{(ncol(params)*nsim)} x \code{ntimes}.
##' Here, \code{nvar} is the number of state variables and \code{ntimes} the length of the argument \code{times}.
##' The \dQuote{obs} element will contain the simulated data, returned as a rank-3 array with dimensions
##' \code{nobs} x \code{(ncol(params)*nsim)} x \code{ntimes}.
##' Here, \code{nobs} is the number of observables.
##'
##' \code{format = "data.frame"} causes the results to be returned as a single data frame containing
##' the time, states, and observations.
##' An ordered factor variable, \sQuote{.id}, distinguishes one simulation from another.
##'
##' @param include.data if \code{TRUE}, the original data and covariates (if any) are included (with \code{.id = "data"}).
##' This option is ignored unless \code{format = "data.frame"}.
##' @param ... additional arguments are passed to \code{\link{pomp}}.
##'
##' @return
##' A single \dQuote{pomp} object,
##' a \dQuote{pompList} object,
##' a named list of two arrays,
##' or a data frame, according to the \code{format} option.
##'
##' If \code{params} is a matrix, each column is treated as a distinct parameter set.
##' In this case, if \code{nsim=1},
##' then \code{simulate} will return one simulation for each parameter set.
##' If \code{nsim>1},
##' then \code{simulate} will yield \code{nsim} simulations for each parameter set.
##' These will be ordered such that
##' the first \code{ncol(params)} simulations represent one simulation
##' from each of the distinct parameter sets,
##' the second \code{ncol(params)} simulations represent a second simulation from each,
##' and so on.
##'
##' Adding column names to \code{params} can be helpful.
##'
NULL

##' @importFrom stats simulate
setGeneric("simulate")

##' @rdname simulate
##' @export
setMethod(
  "simulate",
  signature=signature(object="missing"),
  definition=function (
    object,
    nsim = 1,
    seed = NULL,
    ...,
    times, t0,
    params, rinit, rprocess, rmeasure,
    format = c("pomps", "arrays", "data.frame"),
    include.data = FALSE,
    verbose = getOption("verbose", FALSE)
  ) {

    tryCatch(
      simulate_internal(
        object=NULL,
        ...,
        nsim=nsim,
        seed=seed,
        times=times,
        t0=t0,
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        rmeasure=rmeasure,
        format=match.arg(format),
        include.data=include.data,
        verbose=verbose
      ),
      error = function (e) pStop(who="simulate",conditionMessage(e))
    )

  }
)

##' @rdname simulate
##' @export
setMethod(
  "simulate",
  signature=signature(object="data.frame"),
  definition=function (
    object,
    nsim = 1,
    seed = NULL,
    ...,
    times, t0,
    params, rinit, rprocess, rmeasure,
    format = c("pomps", "arrays", "data.frame"),
    include.data = FALSE,
    verbose = getOption("verbose", FALSE)
  ) {

    tryCatch(
      simulate_internal(
        object,
        ...,
        nsim=nsim,
        seed=seed,
        times=times,
        t0=t0,
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        rmeasure=rmeasure,
        format=match.arg(format),
        include.data=include.data,
        verbose=verbose
      ),
      error = function (e) pStop(who="simulate",conditionMessage(e))
    )

  }
)

##' @rdname simulate
##' @export
setMethod(
  "simulate",
  signature=signature(object="pomp"),
  definition=function (
    object,
    nsim = 1,
    seed = NULL,
    ...,
    format = c("pomps", "arrays", "data.frame"),
    include.data = FALSE,
    verbose = getOption("verbose", FALSE)
  ) {

    tryCatch(
      simulate_internal(
        object,
        ...,
        nsim=nsim,
        seed=seed,
        format=match.arg(format),
        include.data=include.data,
        verbose=verbose
      ),
      error = function (e) pStop(who="simulate",conditionMessage(e))
    )

  }
)

simulate_internal <- function (
  object,
  ...,
  nsim = 1L,
  seed = NULL,
  params,
  format,
  include.data = FALSE,
  .gnsi = TRUE,
  verbose
) {

  object <- pomp(object,...,verbose=verbose)

  if (undefined(object@rprocess)) pStop_(sQuote("rprocess")," is undefined.")

  include.data <- as.logical(include.data)

  if (length(nsim)!=1 || !is.numeric(nsim) || !is.finite(nsim) || nsim < 1)
    pStop_(sQuote("nsim")," must be a positive integer.")
  nsim <- as.integer(nsim)

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)
  if (is.null(params)) params <- numeric(0)
  if (!is.numeric(params))
    pStop_(sQuote("params")," must be named and numeric.")
  params <- as.matrix(params)
  storage.mode(params) <- "double"

  if (ncol(params) == 1L) coef(object) <- params[,1L]

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  return.type <- switch(format,arrays=0L,data.frame=0L,pomps=1L)

  sims <- freeze(
    .Call(P_do_simulate,object,params,nsim,return.type,.gnsi),
    seed=seed
  )

  if (format == "data.frame") {

    nsims <- ncol(sims$states)
    ntimes <- length(time(object))
    simnames <- colnames(sims$states)
    if (is.null(simnames)) simnames <- seq_len(nsims)

    dm <- dim(sims$states)
    nm <- rownames(sims$states)
    dim(sims$states) <- c(dm[1L],prod(dm[-1L]))
    rownames(sims$states) <- nm

    dm <- dim(sims$obs)
    nm <- rownames(sims$obs)
    dim(sims$obs) <- c(dm[1L],prod(dm[-1L]))
    rownames(sims$obs) <- nm

    sims <- cbind(
      time=rep(time(object),each=length(simnames)),
      as.data.frame(t(sims$states)),
      as.data.frame(t(sims$obs)),
      .id=rep(simnames,times=ntimes)
    )

    names(sims)[[1]] <- object@timename

    if (include.data) {
      dat <- as.data.frame(object)
      dat$.id <- "data"
      allnm <- union(names(dat),names(sims))
      for (n in setdiff(allnm,names(dat))) dat[[n]] <- NA
      for (n in setdiff(allnm,names(sims))) sims[[n]] <- rep(dat[[n]],each=nsim)
      sims <- rbind(dat,sims)
      sims$.id <- ordered(sims$.id,levels=c("data",simnames))
    } else {
      sims$.id <- ordered(sims$.id,levels=simnames)
    }
    othernm <- setdiff(names(sims),c(object@timename,".id"))
    sims <- sims[,c(object@timename,".id",othernm)]

  } else if (format == "pomps") {

    if (length(sims) == 1) sims <- sims[[1]]
    else sims <- do.call(concat,sims)

  }

  sims

}
