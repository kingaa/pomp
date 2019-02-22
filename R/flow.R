##' Flow of a deterministic model
##'
##' Get realizations (evolutes) from a dynamic model.
##'
##' In the case of a discrete-time system, the deterministic skeleton is a map and a trajectory is obtained by iterating the map.
##' In the case of a continuous-time system, the deterministic skeleton is a vector-field;
##' \code{flow} uses the numerical solvers in \pkg{\link[deSolve]{deSolve}} to integrate the vectorfield starting with an initial
##' state and for a pre-specified time range.
##'
##' @name flow
##' @rdname flow
##' @include flow.R pomp_class.R skeleton_spec.R
##' @aliases flow flow ,missing-method flow,ANY-method

##' @importFrom deSolve ode diagnostics
##' @importFrom stats setNames
##'
##' @param format the format in which to return the results.
##'
##' \code{format = "array"} causes the evolutes to be returned
##' in a rank-3 array with dimensions
##' \code{nvar} x \code{ncol(xstart)} x \code{ntimes}.
##' Here, \code{nvar} is the number of state variables and \code{ntimes} the length of the argument \code{times}.
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
##' It is not possible to modify the model structure in a call to \code{flow}.
##'
##' @return
##' \code{flow} returns an array of dimensions \code{nvar} x \code{nrep} x \code{ntimes}.
##' If \code{x} is the returned matrix, \code{x[i,j,k]} is the i-th component of the state vector at time \code{times[k]} given parameters \code{params[,j]}.
##'
##' @seealso \code{\link{skeleton}}
NULL

setGeneric(
  "flow",
  function (object, ...)
    standardGeneric("flow")
)

setMethod(
  "flow",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("flow","object")
  }
)

setMethod(
  "flow",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("flow",object)
  }
)

##' @name flow-pomp
##' @aliases flow,pomp-method
##' @rdname flow
##' @export
setMethod(
  "flow",
  signature=signature(object="pomp"),
  definition=function (object, xstart, tstart, params, times,
    format = c("array", "data.frame"), ...,
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      flow.internal(object=object, xstart = xstart, tstart = tstart, params=params,times=times,
        format=format,...,verbose=verbose),
      error = function (e) pStop("flow",conditionMessage(e))
    )

  }
)

flow.internal <- function (object, xstart, tstart, params, times,
  format = c("array", "data.frame"), .gnsi = TRUE, ...,
  verbose) {

  format <- match.arg(format)
  verbose <- as.logical(verbose)

  if(missing(tstart)) reqd_arg("flow","tstart")

  if (missing(times)) times <- reqd_arg("flow","times")
  else times <- as.numeric(times)

  if (length(times)==0)
    pStop_(sQuote("times")," is empty, there is no work to do.")

  if (any(diff(times)<=0))
    pStop_(sQuote("times")," must be a strictly increasing numeric sequence.")

  if (missing(tstart)) tstart <- reqd_arg("flow","tstart")
  else tstart <- as.numeric(tstart)

  if (tstart>times[1L])
    pStop_("the start time ",sQuote("tstart"),
      " must occur no later than the first observation.")
  ntimes <- length(times)

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)
  if (is.null(params)) params <- numeric(0)

  storage.mode(params) <- "double"

  params <- as.matrix(params)
  nrep <- ncol(xstart) # this was initially ncol(params)

  # x0 <- rinit(object,params=params,t0=t0) DELETE LATER
  nvar <- nrow(xstart)
  statenames <- rownames(xstart)
  dim(xstart) <- c(nvar,nrep,1)
  dimnames(xstart) <- list(statenames,NULL,NULL)

  type <- object@skeleton@type          # map or vectorfield?

  pompLoad(object)
  on.exit(pompUnload(object))

  if (type == skeletontype$map) {                  ## MAP

    x <- .Call(P_iterate_map,object,times,tstart,xstart,params,.gnsi)
    .gnsi <- FALSE

  } else if (type == skeletontype$vectorfield) {   ## VECTORFIELD

    znames <- object@accumvars
    if (length(znames)>0) xstart[znames,,] <- 0

    .Call(P_pomp_desolve_setup,object,xstart,params,.gnsi)
    .gnsi <- FALSE

    X <- tryCatch(
      ode(
        y=xstart,
        times=c(tstart,times),
        func="pomp_vf_eval",
        dllname="pomp2",
        initfunc=NULL,
        parms=NULL,
        ...
      ),
      error = function (e) {
        pStop_("error in ODE integrator: ",conditionMessage(e))
      }
    )

    .Call(P_pomp_desolve_takedown)

    if (attr(X,"istate")[1L] != 2)
      pWarn("flow",
        "abnormal exit from ODE integrator, istate = ",attr(X,'istate')[1L])

    if (verbose) diagnostics(X)

    x <- array(data=t(X[-1L,-1L]),dim=c(nvar,nrep,ntimes),
      dimnames=list(statenames,NULL,NULL))

    for (z in znames)
      for (r in seq_len(ncol(x)))
        x[z,r,-1] <- diff(x[z,r,])

  } else {                  ## DEFAULT SKELETON

    x <- array(data=NA_real_,dim=c(nrow(xstart),ncol(xstart),length(times)),
      dimnames=list(rownames(xstart),NULL,NULL))

  }

  dimnames(x) <- setNames(dimnames(x),c("variable","rep",object@timename))
  x
}
