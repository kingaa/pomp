##' Trajectory of a deterministic model
##'
##' Compute trajectories of the deterministic skeleton of a Markov process.
##'
##' In the case of a discrete-time system, the deterministic skeleton is a map and a trajectory is obtained by iterating the map.
##' In the case of a continuous-time system, the deterministic skeleton is a vector-field;
##' \code{trajectory} uses the numerical solvers in \pkg{\link[deSolve]{deSolve}} to integrate the vectorfield.
##'
##' @name trajectory
##' @rdname trajectory
##' @include workhorses.R pomp_class.R flow.R
##' @aliases trajectory,missing-method trajectory,ANY-method
##' @family elementary algorithms
##' @family deterministic methods
##'
##' @importFrom deSolve ode diagnostics
##' @importFrom stats setNames
##'
##' @inheritParams simulate
##' @inheritParams pomp
##'
##' @param object optional;
##' if present, it should be a data frame or a \sQuote{pomp} object.
##' @param ode_control optional list;
##' the elements of this list will be passed to \code{\link[deSolve]{ode}} if the skeleton is a vectorfield, and ignored if it is a map.
##' @param format the format in which to return the results.
##'
##' \code{format = "pomps"} causes the trajectories to be returned as a single \sQuote{pomp} object (if a single parameter vector have been furnished to \code{trajectory}) or as a \sQuote{pompList} object (if multiple parameters have been furnished).
##' In each of these, the \code{states} slot will have been replaced by the computed trajectory.
##' Use \code{\link{states}} to view these.
##' 
##' \code{format = "array"} causes the trajectories to be returned
##' in a rank-3 array with dimensions
##' \code{nvar} x \code{ncol(params)} x \code{ntimes}.
##' Here, \code{nvar} is the number of state variables and \code{ntimes} the length of the argument \code{times}.
##' Thus if \code{x} is the returned array, \code{x[i,j,k]} is the i-th component of the state vector at time \code{times[k]} given parameters \code{params[,j]}.
##'
##' \code{format = "data.frame"} causes the results to be returned as a single data frame containing the time and states.
##' An ordered factor variable, \sQuote{.id}, distinguishes the trajectories from one another.
##'
##' @return
##' The \code{format} option controls the nature of the return value of \code{trajectory}.
##' See above for details.
##'
##' @example examples/trajectory.R
##' @example examples/ricker-bifdiag.R
##'
NULL

setGeneric(
  "trajectory",
  function (object, ...)
    standardGeneric("trajectory")
)

setMethod(
  "trajectory",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("trajectory","object")
  }
)

setMethod(
  "trajectory",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("trajectory",object)
  }
)

##' @rdname trajectory
##' @export
setMethod(
  "trajectory",
  signature=signature(object="missing"),
  definition=function (
    t0, times, params,
    skeleton, rinit,
    ...,
    ode_control = list(),
    format = c("pomps", "array", "data.frame"),
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      trajectory_internal(
        object=NULL,
        t0=t0,
        times=times,
        params=params,
        skeleton=skeleton,
        rinit=rinit,
        ...,
        ode_control=ode_control,
        format=format,
        verbose=verbose
      ),
      error = function (e) pStop(who="trajectory",conditionMessage(e))
    )

  }
)

##' @rdname trajectory
##' @export
setMethod(
  "trajectory",
  signature=signature(object="data.frame"),
  definition=function (
    object,
    ...,
    t0, times, params,
    skeleton, rinit,
    ode_control = list(),
    format = c("pomps", "array", "data.frame"),
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      trajectory_internal(
        object=object,
        t0=t0,
        times=times,
        params=params,
        skeleton=skeleton,
        rinit=rinit,
        ...,
        ode_control=ode_control,
        format=format,
        verbose=verbose
      ),
      error = function (e) pStop(who="trajectory",conditionMessage(e))
    )
    
  }
)

##' @rdname trajectory
##' @export
setMethod(
  "trajectory",
  signature=signature(object="pomp"),
  definition=function (
    object,
    params,
    ...,
    skeleton, rinit,
    ode_control = list(),
    format = c("pomps", "array", "data.frame"),
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      trajectory_internal(
        object=object,
        params=params,
        skeleton=skeleton,
        rinit=rinit,
        ...,
        ode_control=ode_control,
        format=format,
        verbose=verbose
      ),
      error = function (e) pStop(who="trajectory",conditionMessage(e))
    )

  }
)

trajectory_internal <- function (
  object, params,
  ...,
  format = c("pomps", "array", "data.frame"),
  ode_control = list(),
  .gnsi = TRUE,
  verbose
) {

  format <- match.arg(format)
  verbose <- as.logical(verbose)

  object <- pomp(object,...,verbose=verbose)

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)
  if (is.null(params)) params <- numeric(0)
  if (!is.numeric(params))
    pStop_(sQuote("params")," must be named and numeric.")
  params <- as.matrix(params)
  storage.mode(params) <- "double"

  if (ncol(params) == 1) object@params <- params[,1L]

  pompLoad(object)
  on.exit(pompUnload(object))

  x0 <- rinit(object,params=params,verbose=verbose,.gnsi=.gnsi)

  x <- do.call(
    flow,
    c(
      list(object,x0=x0,params=params),
      ode_control,
      list(.gnsi=.gnsi,verbose=verbose)
    )
  )

  repnames <- colnames(x)

  if (format == "data.frame") {
    x <- lapply(
      seq_len(ncol(x)),
      function (k) {
        nm <- rownames(x)
        y <- x[,k,,drop=FALSE]
        dim(y) <- dim(y)[c(1L,3L)]
        y <- as.data.frame(t(y))
        names(y) <- nm
        y[[object@timename]] <- object@times
        y$.id <- as.integer(k)
        y
      }
    )
    x <- do.call(rbind,x)
    x$.id <- ordered(x$.id)
    if (!is.null(repnames)) {
      levels(x$.id) <- repnames
    }

  } else if (format == "pomps") {
    rv <- rep(list(object),ncol(x))
    dy <- dim(x)[c(1L,3L)]
    ny <- dimnames(x)[c(1L,3L)]
    for (k in seq_len(ncol(x))) {
      y <- x[,k,,drop=FALSE]
      dim(y) <- dy
      dimnames(y) <- ny
      rv[[k]]@states <- y
      rv[[k]]@params <- params[,k]
    }
    if (length(rv)>1) {
      names(rv) <- colnames(x)
      x <- do.call(c,rv)
    } else {
      x <- rv[[1]]
    }
  }

  x
}
