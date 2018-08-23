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
##' @include workhorses.R pomp_class.R
##' @aliases trajectory trajectory,missing-method trajectory,ANY-method

##' @importFrom deSolve ode diagnostics
##' @importFrom stats setNames
##'
##' @inheritParams rmeasure
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
    stop("in ",sQuote("trajectory"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "trajectory",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("trajectory")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
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

    trajectory.internal(object=object,params=params,times=times,t0=t0,
      format=format,...,verbose=verbose)

  }
)

trajectory.internal <- function (object, params, times, t0,
  format = c("array", "data.frame"), .getnativesymbolinfo = TRUE, ...,
  verbose) {

  format <- match.arg(format)

  verbose <- as.logical(verbose)

  if (missing(times))
    times <- time(object,t0=FALSE)
  else
    times <- as.numeric(times)

  if (length(times)==0)
    pomp_stop("trajectory",sQuote("times")," is empty, there is no work to do.")

  if (any(diff(times)<=0))
    pomp_stop("trajectory",sQuote("times"),
      " must be an increasing sequence of times.")

  if (missing(t0))
    t0 <- timezero(object)
  else
    t0 <- as.numeric(t0)

  if (t0>times[1L])
    pomp_stop("trajectory","the zero-time ",sQuote("t0"),
      " must occur no later than the first observation.")
  ntimes <- length(times)

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)
  if (is.null(params)) params <- numeric(0)

  if (length(params)==0)
    pomp_stop("trajectory",sQuote("params")," must be supplied.")
  storage.mode(params) <- "double"

  params <- as.matrix(params)
  nrep <- ncol(params)
  paramnames <- rownames(params)
  if (is.null(paramnames))
    pomp_stop("trajectory",sQuote("params")," must have names (or rownames).")

  x0 <- rinit(object,params=params,t0=t0)
  nvar <- nrow(x0)
  statenames <- rownames(x0)
  dim(x0) <- c(nvar,nrep,1)
  dimnames(x0) <- list(statenames,NULL,NULL)

  type <- object@skeleton@type          # map or vectorfield?

  pompLoad(object)
  on.exit(pompUnload(object))

  if (type == 2L) {          ## MAP

    x <- tryCatch(
      .Call(iterate_map,object,times,t0,x0,params,.getnativesymbolinfo),
      error = function (e)
        pomp_stop("trajectory","in map iterator: ",conditionMessage(e))
    )
    .getnativesymbolinfo <- FALSE

  } else if (type == 1L) {   ## VECTORFIELD

    znames <- object@zeronames
    if (length(znames)>0) x0[znames,,] <- 0

    .Call(pomp_desolve_setup,object,x0,params,.getnativesymbolinfo)
    .getnativesymbolinfo <- FALSE

    X <- tryCatch(
      ode(
        y=x0,
        times=c(t0,times),
        func="pomp_vf_eval",
        dllname="pomp",
        initfunc=NULL,
        parms=NULL,
        ...
      ),
      error = function (e) {
        pomp_stop("trajectory","error in ODE integrator: ",conditionMessage(e))
      }
    )

    .Call(pomp_desolve_takedown)

    if (attr(X,"istate")[1L]!=2)
      pomp_warn("trajectory",
        "abnormal exit from ODE integrator, istate = ",attr(X,'istate')[1L])

    if (verbose) diagnostics(X)

    x <- array(data=t(X[-1L,-1L]),dim=c(nvar,nrep,ntimes),
      dimnames=list(statenames,NULL,NULL))

    for (z in znames)
      for (r in seq_len(ncol(x)))
        x[z,r,-1] <- diff(x[z,r,])

  } else {

    x <- array(data=NA_real_,dim=c(dim(x0),length(times)),
      dimnames=list(rownames(x0),NULL,NULL))

  }

  dimnames(x) <- setNames(dimnames(x),c("variable","rep",object@timename))

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
