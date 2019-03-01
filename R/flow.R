##' Flow of a deterministic model
##'
##' Compute the flow induced by a deterministic dynamical system, whether vectorfield or map.
##'
##' In the case of a discrete-time system (map), \code{flow} iterates the map to yield trajectories of the system.
##' In the case of a continuous-time system (vectorfield), \code{flow} uses the numerical solvers in \pkg{\link[deSolve]{deSolve}} to integrate the vectorfield starting from given initial conditions.
##'
##' @name flow
##' @rdname flow
##' @include pomp_class.R skeleton_spec.R workhorses.R
##' @aliases flow flow,missing-method flow,ANY-method
##' @family pomp workhorses
##' 
##' @importFrom deSolve ode diagnostics
##' @importFrom stats setNames
##'
##' @inheritParams dmeasure
##' @inheritParams pomp
##'
##' @param xstart an array with dimensions \code{nvar} x \code{nrep} giving the initial conditions of the trajectories to be computed.
##' 
##' @param times a numeric vector (length \code{ntimes}) containing times.
##' These must be in strictly increasing order.
##' The first entry is assumed to be the time of \code{xstart.}
##' 
##' @param offset a non-negative integer scalar.
##' By default, \code{offset = 0}.
##' See below (Value) for the effect of \code{offset}.
##'
##' @param ... Additional arguments are passed to the ODE integrator (if the skeleton is a vectorfield) and are ignored if it is a map.
##' See \code{\link[deSolve]{ode}} for a description of the additional arguments accepted by the ODE integrator.
##'
##' @return
##' \code{flow} returns an array of dimensions \code{nvar} x \code{nrep} x \code{ntimes-offset}.
##' If \code{x} is the returned matrix, \code{x[i,j,k]} is the i-th component of the state vector at time \code{times[k+offset]} given parameters \code{params[,j]}.
##'
##' @seealso \code{\link{skeleton}}, \code{\link{trajectory}}, \code{\link{rprocess}}
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
  definition=function (object, xstart, params, times, offset = 0, ...,
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      flow.internal(object=object,x0=xstart,params=params,times=times,
        offset=offset,...,verbose=verbose),
      error = function (e) pStop("flow",conditionMessage(e))
    )
    
  }
)

flow.internal <- function (object, x0, params, times, offset, ...,
  .gnsi = TRUE, verbose) {

  verbose <- as.logical(verbose)

  if (missing(times))
    reqd_arg("flow","times")
  else
    times <- as.numeric(times)

  if (length(times)==0)
    pStop_(sQuote("times")," is empty, there is no work to do.")

  if (any(diff(times)<0))
    pStop_(sQuote("times")," must be a non-decreasing numeric sequence.")

  offset <- as.integer(offset)
  if (offset < 0L)
    pStop_(sQuote("offset")," must be a non-negative integer.")

  params <- as.matrix(params)
  nrep <- ncol(params)

  x0 <- as.matrix(x0)
  nvar <- nrow(x0)

  storage.mode(params) <- "double"
  storage.mode(x0) <- "double"

  statenames <- rownames(x0)
  dim(x0) <- c(nvar,nrep,1)
  dimnames(x0) <- list(statenames,NULL,NULL)

  t0 <- times[1L]
  if (offset == 1L)
    times <- times[-1L]
  else if (offset > 1L) 
    times <- times[seq.int(offset+1L,length(times),by=1L)]

  ntimes <- length(times)

  type <- object@skeleton@type          # map or vectorfield?

  pompLoad(object)
  on.exit(pompUnload(object))

  if (type == skeletontype$map) {                  ## MAP

    x <- .Call(P_iterate_map,object,times,t0,x0,params,.gnsi)
    .gnsi <- FALSE

  } else if (type == skeletontype$vectorfield) {   ## VECTORFIELD

    znames <- object@accumvars
    if (length(znames)>0) x0[znames,,] <- 0

    .Call(P_pomp_desolve_setup,object,x0,params,.gnsi)
    .gnsi <- FALSE

    X <- tryCatch(
      ode(
        y=x0,
        times=c(t0,times),
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
      pWarn("trajectory",
        "abnormal exit from ODE integrator, istate = ",attr(X,'istate')[1L])

    if (verbose) diagnostics(X)

    x <- array(data=t(X[-1L,-1L]),dim=c(nvar,nrep,ntimes),
      dimnames=list(statenames,NULL,NULL))

    for (z in znames)
      for (r in seq_len(ncol(x)))
        x[z,r,-1] <- diff(x[z,r,])

  } else {                  ## DEFAULT SKELETON

    x <- array(data=NA_real_,dim=c(nrow(x0),ncol(x0),length(times)),
      dimnames=list(rownames(x0),NULL,NULL))

  }

  dimnames(x) <- setNames(dimnames(x),c("variable","rep",object@timename))

  x
}
