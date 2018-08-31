##' Workhorse functions for the \pkg{pomp} algorithms.
##'
##' These functions mediate the interface between the user's model and the package algorithms.
##' They are low-level functions that do the work needed by the package's inference methods.
##'
##' They include \describe{
##' \item{\code{\link{dmeasure}}}{which evaluates the measurement model density,}
##' \item{\code{\link{rmeasure}}}{which samples from the measurement model distribution,}
##' \item{\code{\link{dprocess}}}{which evaluates the process model density,}
##' \item{\code{\link{rprocess}}}{which samples from the process model distribution,}
##' \item{\code{\link{dprior}}}{which evaluates the prior probability density,}
##' \item{\code{\link{rprior}}}{which samples from the prior distribution,}
##' \item{\code{\link{skeleton}}}{which evaluates the model's deterministic skeleton,}
##' \item{\code{\link{partrans}}}{which performs parameter transformations associated with the model.}
##' }
##'
##' @name workhorses
##' @include pomp_class.R pomp_fun.R load.R pstop.R
##' @docType methods
##' @family pomp workhorses
##' @seealso \code{\link[=simulate-pomp]{simulate}}, \code{\link{trajectory}},
##' \code{\link{pfilter}}, \code{\link{probe}}
##'
##' @author Aaron A. King
##'
##' @keywords programming internal
NULL

##' dmeasure
##'
##' \code{dmeasure} evaluates the probability density of observations given states.
##'
##' @name dmeasure
##' @docType methods
##' @aliases dmeasure dmeasure,pomp-method dmeasure-pomp
##' dmeasure,ANY-method dmeasure,missing-method
##' @keywords internal
##' @family pomp workhorses
##'
##' @param object an object of class \sQuote{pomp}, or of a class that extends \sQuote{pomp}.
##' This will typically be the output of \code{pomp}, \code{simulate}, or one of the \pkg{pomp} inference algorithms.
##'
##' @param x an array containing states of the unobserved process.
##' The dimensions of \code{x} are \code{nvars} x \code{nrep} x \code{ntimes},
##' where \code{nvars} is the number of state variables,
##' \code{nrep} is the number of replicates,
##' and \code{ntimes} is the length of \code{times}.
##' One can also pass \code{x} as a named numeric vector, which is equivalent to the \code{nrep=1}, \code{ntimes=1} case.
##'
##' @param y a matrix containing observations.
##' The dimensions of \code{y} are \code{nobs} x \code{ntimes}, where \code{nobs} is the number of observables
##' and \code{ntimes} is the length of \code{times}.
##'
##' @param times a numeric vector (length \code{ntimes}) containing times.
##' These must be in non-decreasing order.
##'
##' @param params a \code{npar} x \code{nrep} matrix of parameters.
##' Each column is treated as an independent parameter set, in correspondence with the corresponding column of \code{x}.
##'
##' @param log if TRUE, log probabilities are returned.
##'
##' @param \dots additional arguments are ignored.
##'
##' @return
##' \code{dmeasure} returns a matrix of dimensions \code{nreps} x \code{ntimes}.
##' If \code{d} is the returned matrix, \code{d[j,k]} is the likelihood (or log likelihood if \code{log = TRUE}) of the observation \code{y[,k]} at time \code{times[k]} given the state \code{x[,j,k]}.
##'
NULL

setGeneric(
  "dmeasure",
  function (object, ...)
    standardGeneric("dmeasure")
)

setMethod(
  "dmeasure",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("dmeasure","object")
  }
)

setMethod(
  "dmeasure",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("dmeasure",object)
  }
)

##' @export
##' @rdname dmeasure
setMethod(
  "dmeasure",
  signature=signature(object="pomp"),
  definition=function (object, y, x, times, params, ..., log = FALSE) {
    tryCatch(
      dmeasure.internal(object=object,y=y,x=x,times=times,
        params=params,log=log,...),
      error = function (e) pStop("dmeasure",conditionMessage(e))
    )
  }
)

dmeasure.internal <- function (object, y, x, times, params, ..., log = FALSE,
  .getnativesymbolinfo = TRUE) {
  storage.mode(y) <- "double"
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_dmeasure,object,y,x,times,params,log,.getnativesymbolinfo)
}

##' dprior
##'
##' Evaluates the prior probability density.
##'
##' @name dprior
##' @docType methods
##' @aliases dprior dprior,pomp-method dprior-pomp
##' dprior,ANY-method dprior,missing-method
##' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##'
##' @return
##' The required density (or log density), as a numeric vector.
##'
NULL

setGeneric(
  "dprior",
  function (object, ...)
    standardGeneric("dprior")
)

setMethod(
  "dprior",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("dprior","object")
  }
)

setMethod(
  "dprior",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("dprior",object)
  }
)

##' @export
##' @rdname dprior
setMethod(
  "dprior",
  signature=signature(object="pomp"),
  definition=function (object, params, ..., log = FALSE) {
    tryCatch(
      dprior.internal(object=object,params=params,log=log,...),
      error = function (e) pStop("dprior",conditionMessage(e))
    )
  }
)

dprior.internal <- function (object, params, log = FALSE,
  .getnativesymbolinfo = TRUE, ...) {
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_dprior,object,params,log,.getnativesymbolinfo)
}

##' dprocess
##'
##' Evaluates the probability density of a sequence of consecutive state transitions.
##'
##' @name dprocess
##' @docType methods
##' @aliases dprocess dprocess,pomp-method dprocess-pomp
##' dprocess,ANY-method dprocess,missing-method
##' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##'
##' @return
##' \code{dprocess} returns a matrix of dimensions \code{nrep} x \code{ntimes-1}.
##' If \code{d} is the returned matrix, \code{d[j,k]} is the likelihood (or the log likelihood if \code{log=TRUE}) of the transition from state \code{x[,j,k-1]} at time \code{times[k-1]} to state \code{x[,j,k]} at time \code{times[k]}.
##'
NULL

setGeneric(
  "dprocess",
  function (object, ...)
    standardGeneric("dprocess")
)

setMethod(
  "dprocess",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("dprocess","object")
  }
)

setMethod(
  "dprocess",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("dprocess",object)
  }
)

##' @export
##' @rdname dprocess
setMethod(
  "dprocess",
  signature=signature(object="pomp"),
  definition = function (object, x, times, params, ..., log = FALSE) {
    tryCatch(
      dprocess.internal(object=object,x=x,times=times,params=params,log=log,...),
      error = function (e) pStop("dprocess",conditionMessage(e))
    )
  }
)

dprocess.internal <- function (object, x, times, params, log = FALSE, .getnativesymbolinfo = TRUE, ...) {
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_dprocess,object,x,times,params,log,.getnativesymbolinfo)
}

##' partrans
##'
##' Performs parameter transformations.
##'
##' @name partrans
##' @docType methods
##' @aliases  partrans partrans,pomp-method partrans-pomp
##' partrans,ANY-method partrans,missing-method
##' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##' @param dir the direction of the transformation to perform.
##'
##' @return
##' If \code{dir=fromEst}, the parameters in \code{params} are assumed to be on the estimation scale and are transformed onto the natural scale.
##' If \code{dir=toEst}, they are transformed onto the estimation scale.
##' In both cases, the parameters are returned as a named numeric vector or an array with rownames, as appropriate.
##'
NULL

setGeneric(
  "partrans",
  function (object, ...)
    standardGeneric("partrans")
)

setMethod(
  "partrans",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("partrans","object")
  }
)

setMethod(
  "partrans",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("partrans",object)
  }
)

##' @export
##' @rdname partrans
setMethod(
  "partrans",
  signature=signature(object="pomp"),
  definition=function (object, params, dir = c("fromEst", "toEst"), ...) {
    dir <- match.arg(dir)
    tryCatch(
      partrans.internal(object=object,params=params,dir=dir,...),
      error = function (e) pStop("partrans",conditionMessage(e))
    )
  }
)

partrans.internal <- function (object, params, dir = c("fromEst", "toEst"),
  .getnativesymbolinfo = TRUE, ...) {
  if (object@partrans@has) {
    dir <- switch(dir,fromEst=-1L,toEst=1L)
    pompLoad(object)
    on.exit(pompUnload(object))
    params <- .Call(P_do_partrans,object,params,dir,.getnativesymbolinfo)
  }
  params
}

##' rinit
##'
##' Samples from the initial-state distribution.
##'
##' @name rinit
##' @docType methods
##' @aliases rinit rinit,pomp-method rinit-pomp
##' rinit,ANY-method rinit,missing-method
##' @keywords internal
##' @family pomp workhorses
##' @seealso Specification of the initial-state distribution: \code{\link{rinit_spec}}
##'
##' @inheritParams dmeasure
##' @param t0 the initial time, i.e., the time corresponding to the initial-state distribution.
##' @param nsim optional integer; the number of initial states to simulate per column of \code{params}.
##'
##' @return
##' \code{rinit} returns an \code{nvar} x \code{nsim*ncol(params)} matrix of state-process initial conditions when given an \code{npar} x \code{nsim} matrix of parameters, \code{params}, and an initial time \code{t0}.
##' By default, \code{t0} is the initial time defined when the \sQuote{pomp} object ws constructed.
##'
NULL

setGeneric(
  "rinit",
  function (object, ...)
    standardGeneric("rinit")
)

setMethod(
  "rinit",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("rinit","object")
  }
)

setMethod(
  "rinit",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("rinit",object)
  }
)

##' @export
##' @rdname rinit
setMethod(
  "rinit",
  signature=signature("pomp"),
  definition=function (object, params, t0, nsim = 1, ...) {
    tryCatch(
      rinit.internal(object=object,params=params,t0=t0,nsim=nsim,...),
      error = function (e) pStop("rinit",conditionMessage(e))
    )
  }
)

rinit.internal <- function (object, params, t0, nsim = 1,
  .getnativesymbolinfo = TRUE, ...) {
  if (missing(t0)) t0 <- object@t0
  if (missing(params)) params <- coef(object)
  else storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_rinit,object,params,t0,nsim,.getnativesymbolinfo)
}

##' rmeasure
##'
##' Sample from the measurement model distribution, given values of the latent states and the parameters.
##'
##' @name rmeasure
##' @docType methods
##' @aliases rmeasure rmeasure,pomp-method rmeasure-pomp
##' rmeasure,ANY-method rmeasure,missing-method
###' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##'
##' @return
##' \code{rmeasure} returns a rank-3 array of dimensions
##' \code{nobs} x \code{nrep} x \code{ntimes},
##' where \code{nobs} is the number of observed variables.
##'
NULL

setGeneric(
  "rmeasure",
  function (object, ...)
    standardGeneric("rmeasure")
)

setMethod(
  "rmeasure",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("rmeasure","object")
  }
)

setMethod(
  "rmeasure",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("rmeasure",object)
  }
)

##' @export
##' @rdname rmeasure
setMethod(
  "rmeasure",
  signature=signature(object="pomp"),
  definition=function (object, x, times, params, ...) {
    tryCatch(
      rmeasure.internal(object=object,x=x,times=times,params=params,...),
      error = function (e) pStop("rmeasure",conditionMessage(e))
    )
  }
)

rmeasure.internal <- function (object, x, times, params,
  .getnativesymbolinfo = TRUE, ...) {
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_rmeasure,object,x,times,params,.getnativesymbolinfo)
}

##' rprior
##'
##' Sample from the prior probability distribution.
##'
##' @name rprior
##' @docType methods
##' @aliases rprior rprior,pomp-method rprior-pomp
##' rprior,ANY-method rprior,missing-method
##' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##'
##' @return
##' A numeric matrix containing the required samples.
##'
NULL

setGeneric(
  "rprior",
  function (object, ...)
    standardGeneric("rprior")
)

setMethod(
  "rprior",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("rprior","object")
  }
)

setMethod(
  "rprior",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("rprior",object)
  }
)

##' @export
##' @rdname rprior
setMethod(
  "rprior",
  signature=signature(object="pomp"),
  definition=function (object, params, ...)
    tryCatch(
      rprior.internal(object=object,params=params,...),
      error = function (e) pStop("rprior",conditionMessage(e))
    )
)

rprior.internal <- function (object, params, .getnativesymbolinfo = TRUE, ...) {
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_rprior,object,params,.getnativesymbolinfo)
}

##' rprocess
##'
##' \code{rprocess} simulates the process-model portion of partially-observed Markov process.
##'
##' When \code{rprocess} is called, the first entry of \code{times} is taken to be the initial time (i.e., that corresponding to \code{xstart}).
##' Subsequent times are the additional times at which the state of the simulated processes are required.
##'
##' @name rprocess
##' @docType methods
##' @aliases rprocess rprocess,pomp-method rprocess-pomp
##' rprocess,ANY-method rprocess,missing-method
##' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##'
##' @param xstart an \code{nvar} x \code{nrep} matrix containing the starting state of the system.
##' Columns of \code{xstart} correspond to states;
##' rows to components of the state vector.
##' One independent simulation will be performed for each column.
##' Note that in this case, \code{params} must also have \code{nrep} columns.
##'
##' @param params a \code{npar} x \code{nrep} matrix of parameters.
##' Each column is treated as an independent parameter set, in correspondence with the corresponding column of \code{xstart}.
##'
##' @param offset integer; the first \code{offset} times in \code{times} will be discarded.
##'
##' @return
##' \code{rprocess} returns a rank-3 array with rownames.
##' Suppose \code{x} is the array returned.
##' Then \preformatted{dim(x)=c(nvars,nrep,ntimes-offset),}
##' where \code{nvars} is the number of state variables (=\code{nrow(xstart)}),
##' \code{nrep} is the number of independent realizations simulated (=\code{ncol(xstart)}), and
##' \code{ntimes} is the length of the vector \code{times}.
##' \code{x[,j,k]} is the value of the state process in the \code{j}-th realization at time \code{times[k+offset]}.
##' The rownames of \code{x} must correspond to those of \code{xstart}.
##'
NULL

setGeneric(
  "rprocess",
  function (object, ...)
    standardGeneric("rprocess")
)

setMethod(
  "rprocess",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("rprocess","object")
  }
)

setMethod(
  "rprocess",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("rprocess",object)
  }
)

##' @export
##' @rdname rprocess
setMethod(
  "rprocess",
  signature=signature(object="pomp"),
  definition=function (object, xstart, times, params, offset = 0L, ...) {
    tryCatch(
      rprocess.internal(object=object,xstart=xstart,times=times,params=params,
        offset=offset,...),
      error = function (e) pStop("rprocess",conditionMessage(e))
    )
  }
)

rprocess.internal <- function (object, xstart, times, params, offset = 0L,
  .getnativesymbolinfo = TRUE, ...) {
  storage.mode(xstart) <- "double"
  storage.mode(params) <- "double"
  if (!is.finite(offset)) pStop_("invalid ",sQuote("offset"),".")
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_rprocess,object,xstart,times,params,offset,.getnativesymbolinfo)
}

##' skeleton
##'
##' Evaluates the deterministic skeleton at a point or points in state space, given parameters.
##' In the case of a discrete-time system, the skeleton is a map.
##' In the case of a continuous-time system, the skeleton is a vectorfield.
##' NB: \code{skeleton} just evaluates the deterministic skeleton;
##' it does not iterate or integrate (see \code{\link{trajectory}} for this).
##'
##' @name skeleton
##' @docType methods
##' @aliases skeleton skeleton,pomp-method skeleton-pomp
##' skeleton,ANY-method skeleton,missing-method
##' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##'
##' @return
##' \code{skeleton} returns an array of dimensions \code{nvar} x \code{nrep} x \code{ntimes}.
##' If \code{f} is the returned matrix, \code{f[i,j,k]} is the i-th component of the deterministic skeleton at time \code{times[k]} given the state \code{x[,j,k]} and parameters \code{params[,j]}.
##'
NULL

setGeneric(
  "skeleton",
  function (object, ...)
    standardGeneric("skeleton")
)

setMethod(
  "skeleton",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("skeleton","object")
  }
)

setMethod(
  "skeleton",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("skeleton",object)
  }
)

##' @export
##' @rdname skeleton
setMethod(
  "skeleton",
  signature=signature("pomp"),
  definition=function (object, x, times, params, ...)
    tryCatch(
      skeleton.internal(object=object,x=x,times=times,params=params,...),
      error = function (e) pStop("skeleton",conditionMessage(e))
    )
)

skeleton.internal <- function (object, x, times, params, .getnativesymbolinfo = TRUE, ...) {
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_skeleton,object,x,times,params,.getnativesymbolinfo)
}
