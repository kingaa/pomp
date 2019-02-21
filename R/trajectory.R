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
##' @include workhorses.R pomp_class.R skeleton_spec.R
##' @aliases trajectory trajectory,missing-method trajectory,ANY-method

##' @importFrom deSolve ode diagnostics
##' @importFrom stats setNames
##'
##' @inheritParams dmeasure
##' @inheritParams rinit
##'
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
##' @param times a numeric vector (length \code{ntimes}) containing times.
##' These must be in strictly increasing order.
##'
##' @param verbose logical; if \code{TRUE}, more information will be displayed.
##'
##' @param \dots Additional arguments are passed to the ODE integrator (if the skeleton is a vectorfield) and are ignored if it is a map.
##' See \code{\link[deSolve]{ode}} for a description of the additional arguments accepted by the ODE integrator.
##'
##' Note that this behavior differs from most other functions in \pkg{pomp}.
##' It is not possible to modify the model structure in a call to \code{trajectory}.
##'
##' @return
##' \code{trajectory} returns an array of dimensions \code{nvar} x \code{nrep} x \code{ntimes}.
##' If \code{x} is the returned matrix, \code{x[i,j,k]} is the i-th component of the state vector at time \code{times[k]} given parameters \code{params[,j]}.
##'
##' @seealso \code{\link{skeleton}}
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

##' @name trajectory-pomp
##' @aliases trajectory,pomp-method
##' @rdname trajectory
##' @export
setMethod(
  "trajectory",
  signature=signature(object="pomp"),
  definition=function (object, params, times, t0,
    format = c("array", "data.frame"), ...,
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      trajectory.internal(object=object,params=params,times=times,t0=t0,
        format=format,...,verbose=verbose),
      error = function (e) pStop("trajectory",conditionMessage(e))
    )

  }
)

trajectory.internal <- function (object, params, times, t0,
  format = c("array", "data.frame"), .gnsi = TRUE, ...,
  verbose) {

  if (missing(t0)) t0 <- timezero(object)
  else t0 <- as.numeric(t0)

  if (missing(times)) times <- time(object,t0=FALSE)
  else times <- as.numeric(times)

  x0 <- rinit(object,params=params,t0=t0)

  x <- flow(object, x0, t0, params, times,
          format = c("array", "data.frame"), ...,
          verbose = getOption("verbose", FALSE))

  if (format == "data.frame") {
    x <- lapply(
      seq_len(ncol(x)),
      function (k) {
        nm <- rownames(x)
        y <- x[,k,,drop=FALSE]
        dim(y) <- dim(y)[c(1L,3L)]
        y <- as.data.frame(t(y))
        names(y) <- nm
        y[[object@timename]] <- times
        y$.id <- as.integer(k)
        y
      }
    )
    x <- do.call(rbind,x)
    x$.id <- ordered(x$.id)
  }

  x
}
