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
##' @inheritParams pomp
##'
##' @param object optional;
##' if present, it should be a data frame or a \sQuote{pomp} object.
##' @param ode_control optional list;
##' the elements of this list will be passed to \code{\link[=deSolve]{ode}} if the skeleton is a vectorfield, and ignored if it is a map.
##' @param format the format in which to return the results.
##'
##' \code{format = "array"} causes the trajectories to be returned
##' in a rank-3 array with dimensions
##' \code{nvar} x \code{ncol(params)} x \code{ntimes}.
##' Here, \code{nvar} is the number of state variables and \code{ntimes} the length of the argument \code{times}.
##'
##' \code{format = "data.frame"} causes the results to be returned as a single data frame containing
##' the time and states.
##' An ordered factor variable, \sQuote{.id}, distinguishes the trajectories from one another.
##'
##' @return
##' \code{trajectory} returns an array of dimensions \code{nvar} x \code{nrep} x \code{ntimes}.
##' If \code{x} is the returned matrix, \code{x[i,j,k]} is the i-th component of the state vector at time \code{times[k]} given parameters \code{params[,j]}.
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
    times, t0, params,
    skeleton, rinit,
    ...,
    ode_control = list(),
    format = c("array", "data.frame"),
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      trajectory.internal(
        object=NULL,
        times=times,
        t0=t0,
        params=params,
        skeleton=skeleton,
        rinit=rinit,
        ...,
        ode_control=ode_control,
        format=format,
        verbose=verbose
      ),
      error = function (e) pStop("trajectory",conditionMessage(e))
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
    times, t0, params,
    skeleton, rinit,
    ode_control = list(),
    format = c("array", "data.frame"),
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      trajectory.internal(
        object=object,
        params=params,
        times=times,
        t0=t0,
        skeleton=skeleton,
        rinit=rinit,
        ...,
        ode_control=ode_control,
        format=format,
        verbose=verbose
      ),
      error = function (e) pStop("trajectory",conditionMessage(e))
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
    format = c("array", "data.frame"),
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      trajectory.internal(
        object=object,
        params=params,
        skeleton=skeleton,
        rinit=rinit,
        ...,
        ode_control=ode_control,
        format=format,
        verbose=verbose
      ),
      error = function (e) pStop("trajectory",conditionMessage(e))
    )

  }
)

trajectory.internal <- function (
  object, params,
  ...,
  format = c("array", "data.frame"),
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
  }

  x
}
