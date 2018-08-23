##' Workhorse functions for the \pkg{pomp} algorithms.
##'
##' These functions mediate the interface between the user's model and the package algorithms.
##' They are low-level functions that do the work needed by the package's inference methods.
##'
##' They include \describe{
##' \item{dmeasure}{which evaluates the measurement model density,}
##' \item{rmeasure}{which samples from the measurement model distribution,}
##' \item{dprocess}{which evaluates the process model density,}
##' \item{rprocess}{which samples from the process model distribution,}
##' \item{dprior}{which evaluates the prior probability density,}
##' \item{rprior}{which samples from the prior distribution,}
##' \item{skeleton}{which evaluates the model's deterministic skeleton,}
##' \item{partrans}{which performs parameter transformations associated with the model.}
##' }
##'
##' @name workhorses
##' @include pomp_class.R pomp_fun.R load.R
##' @docType methods
##' @family pomp workhorses
##' @seealso \code{\link[=simulate-pomp]{simulate}} \code{\link{trajectory}}
##'
##' @author Aaron A. King
##'
##' @keywords programming internal
NULL

##' Generic dmeasure
##'
##' @name dmeasure-generic
##' @aliases dmeasure,missing-method dmeasure,ANY-method
##' @keywords internal
##'
setGeneric(
  "dmeasure",
  function (object, ...)
    standardGeneric("dmeasure")
)

setMethod(
  "dmeasure",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("dmeasure"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "dmeasure",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("dmeasure")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' dmeasure
##'
##' \code{dmeasure} evaluates the probability density of observations given states.
##'
##' @name dmeasure
##' @aliases dmeasure dmeasure,pomp-method dmeasure-pomp
##' @keywords internal
##' @family pomp workhorses
##'
##' @param object an object of class \sQuote{pomp}, or of a class that extends \sQuote{pomp}.
##' @param x a rank-3 array containing states of the unobserved process.
##' The dimensions of \code{x} are \code{nvars} x \code{nrep} x \code{ntimes},
##' where \code{nvars} is the number of state variables,
##' \code{nrep} is the number of replicates,
##' and \code{ntimes} is the length of \code{times}.
##' @param y a matrix containing observations.
##' The dimensions of \code{y} are \code{nobs} x \code{ntimes}, where \code{nobs} is the number of observables
##' and \code{ntimes} is the length of \code{times}.
##' @param times a numeric vector (length \code{ntimes}) containing times.
##' These must be in non-decreasing order.
##' @param params a \code{npar} x \code{nrep} matrix of parameters.
##' Each column is an independent parameter set and is paired with the corresponding column of \code{x} or \code{xstart}.
##' @param log if TRUE, log probabilities are returned.
##' @param \dots additional arguments are ignored.
##'
##' @return
##' \code{dmeasure} returns a matrix of dimensions \code{nreps} x \code{ntimes}.
##' If \code{d} is the returned matrix, \code{d[j,k]} is the likelihood (or log likelihood if \code{log = TRUE}) of the observation \code{y[,k]} at time \code{times[k]} given the state \code{x[,j,k]}.
##'
##' @export
setMethod(
  "dmeasure",
  signature=signature(object="pomp"),
  definition=function (object, y, x, times, params, log = FALSE, ...) {
    tryCatch(
      dmeasure.internal(object=object,y=y,x=x,times=times,
        params=params,log=log,...),
      error = function (e) pomp_stop("dmeasure",conditionMessage(e))
    )
  }
)

dmeasure.internal <- function (object, y, x, times, params, log = FALSE,
  .getnativesymbolinfo = TRUE, ...) {
  storage.mode(y) <- "double"
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(do_dmeasure,object,y,x,times,params,log,.getnativesymbolinfo)
}

##' Generic dprior
##'
##' @name dprior-generic
##' @aliases dprior,missing-method dprior,ANY-method
##' @keywords internal
setGeneric(
  "dprior",
  function (object, ...)
    standardGeneric("dprior")
)

setMethod(
  "dprior",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("dprior"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "dprior",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("dprior")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' dprior
##'
##' Evaluates the prior probability density.
##'
##' @name dprior
##' @aliases dprior dprior,pomp-method dprior-pomp
##' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##'
##' @return
##' The required density (or log density), as a numeric vector.
##'
##' @export
setMethod(
  "dprior",
  signature=signature(object="pomp"),
  definition=function (object, params, log = FALSE, ...) {
    tryCatch(
      dprior.internal(object=object,params=params,log=log,...),
      error = function (e) pomp_stop("dprior",conditionMessage(e))
    )
  }
)

dprior.internal <- function (object, params, log = FALSE,
  .getnativesymbolinfo = TRUE, ...) {
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(do_dprior,object,params,log,.getnativesymbolinfo)
}

##' Generic dprocess
##'
##' @name dprocess-generic
##' @aliases dprocess,missing-method dprocess,ANY-method
##' @keywords internal
##' @inheritParams dmeasure

setGeneric(
  "dprocess",
  function (object, ...)
    standardGeneric("dprocess")
)

setMethod(
  "dprocess",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("dprocess"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "dprocess",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("dprocess")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' dprocess
##'
##' Evaluates the probability density of a sequence of consecutive state transitions.
##'
##' @name dprocess
##' @aliases dprocess dprocess,pomp-method dprocess-pomp
##' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##'
##' @return
##' \code{dprocess} returns a matrix of dimensions \code{nrep} x \code{ntimes-1}.
##' If \code{d} is the returned matrix, \code{d[j,k]} is the likelihood (or the log likelihood if \code{log=TRUE}) of the transition from state \code{x[,j,k-1]} at time \code{times[k-1]} to state \code{x[,j,k]} at time \code{times[k]}.
##'
##' @export
setMethod(
  "dprocess",
  signature=signature(object="pomp"),
  definition = function (object, x, times, params, log = FALSE, ...) {
    tryCatch(
      dprocess.internal(object=object,x=x,times=times,params=params,log=log,...),
      error = function (e) pomp_stop("dprocess",conditionMessage(e))
    )
  }
)

dprocess.internal <- function (object, x, times, params, log = FALSE, .getnativesymbolinfo = TRUE, ...) {
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(do_dprocess,object,x,times,params,log,.getnativesymbolinfo)
}

##' Generic partrans
##'
##' @name partrans-generic
##' @aliases partrans,missing-method partrans,ANY-method
##' @keywords internal
##'
setGeneric(
  "partrans",
  function (object, ...)
    standardGeneric("partrans")
)

setMethod(
  "partrans",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("partrans"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "partrans",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("partrans")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' partrans
##'
##' Performs parameter transformations.
##'
##' @name partrans
##' @aliases partrans partrans,pomp-method partrans-pomp
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
##' @export
setMethod(
  "partrans",
  signature=signature(object="pomp"),
  definition=function (object, params, dir = c("fromEst", "toEst"), ...) {
    dir <- match.arg(dir)
    tryCatch(
      partrans.internal(object=object,params=params,dir=dir,...),
      error = function (e) pomp_stop("partrans",conditionMessage(e))
    )
  }
)

partrans.internal <- function (object, params, dir = c("fromEst", "toEst"),
  .getnativesymbolinfo = TRUE, ...) {
  if (object@partrans@has) {
    dir <- switch(dir,fromEst=-1L,toEst=1L)
    pompLoad(object)
    on.exit(pompUnload(object))
    params <- .Call(do_partrans,object,params,dir,.getnativesymbolinfo)
  }
  params
}

##' Generic rinit
##'
##' @name rinit-generic
##' @aliases rinit,missing-method rinit,ANY-method
##' @keywords internal
setGeneric(
  "rinit",
  function (object, ...)
    standardGeneric("rinit")
)

setMethod(
  "rinit",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("rinit"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "rinit",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("rinit")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' rinit
##'
##' Samples from the initial-state distribution.
##'
##' @name rinit
##' @aliases rinit rinit,pomp-method rinit-pomp
##' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##' @param t0 the initial time, i.e., the time corresponding to the initial-state distribution.
##' @param nsim optional integer; the number of initial states to simulate per column of \code{params}.
##'
##' @return
##' \code{rinit} returns an \code{nvar} x \code{nsim*ncol(params)} matrix of state-process initial conditions when given an \code{npar} x \code{nsim} matrix of parameters, \code{params}, and an initial time \code{t0}.
##' By default, \code{t0} is the initial time defined when the \sQuote{pomp} object ws constructed.
##'
##' @export
setMethod(
  "rinit",
  signature=signature("pomp"),
  definition=function (object, params, t0, nsim = 1, ...) {
    tryCatch(
      rinit.internal(object=object,params=params,t0=t0,nsim=nsim,...),
      error = function (e) pomp_stop("rinit",conditionMessage(e))
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
  .Call(do_rinit,object,params,t0,nsim,.getnativesymbolinfo)
}

##' Generic rmeasure
##'
##' @name rmeasure-generic
##' @aliases rmeasure,missing-method rmeasure,ANY-method
##' @keywords internal
setGeneric(
  "rmeasure",
  function (object, ...)
    standardGeneric("rmeasure")
)

setMethod(
  "rmeasure",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("rmeasure"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "rmeasure",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("rmeasure")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' rmeasure
##'
##' Sample from the measurement model distribution, given values of the latent states and the parameters.
##'
##' @name rmeasure
##' @aliases rmeasure rmeasure,pomp-method rmeasure-pomp
##' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##'
##' @return
##' \code{rmeasure} returns a rank-3 array of dimensions \code{nobs} x \code{nrep} x \code{ntimes}, where \code{nobs} is the number of observed variables.
##'
##' @export
setMethod(
  "rmeasure",
  signature=signature(object="pomp"),
  definition=function (object, x, times, params, ...) {
    tryCatch(
      rmeasure.internal(object=object,x=x,times=times,params=params,...),
      error = function (e) pomp_stop("rmeasure",conditionMessage(e))
    )
  }
)

rmeasure.internal <- function (object, x, times, params,
  .getnativesymbolinfo = TRUE, ...) {
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(do_rmeasure,object,x,times,params,.getnativesymbolinfo)
}

##' Generic rprior
##'
##' @name rprior-generic
##' @aliases rprior,missing-method rprior,ANY-method
##' @keywords internal
setGeneric(
  "rprior",
  function (object, ...)
    standardGeneric("rprior")
)

setMethod(
  "rprior",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("rprior"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "rprior",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("rprior")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' rprior
##'
##' Sample from the prior probability distribution.
##'
##' @name rprior
##' @aliases rprior rprior,pomp-method rprior-pomp
##' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##'
##' @return
##' A numeric matrix containing the required samples.
##'
##' @export
setMethod(
  "rprior",
  signature=signature(object="pomp"),
  definition=function (object, params, ...)
    tryCatch(
      rprior.internal(object=object,params=params,...),
      error = function (e) pomp_stop("rprior",conditionMessage(e))
    )
)

rprior.internal <- function (object, params, .getnativesymbolinfo = TRUE, ...) {
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(do_rprior,object,params,.getnativesymbolinfo)
}

##' Generic rprocess
##'
##' @name rprocess-generic
##' @aliases rprocess,missing-method rprocess,ANY-method
##' @keywords internal
setGeneric(
  "rprocess",
  function (object, ...)
    standardGeneric("rprocess")
)

setMethod(
  "rprocess",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("rprocess"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "rprocess",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("rprocess")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' rprocess
##'
##' \code{rprocess} simulates the process-model portion of partially-observed Markov process.
##'
##' When \code{rprocess} is called, the first entry of \code{times} is taken to be the initial time (i.e., that corresponding to \code{xstart}).
##' Subsequent times are the additional times at which the state of the simulated processes are required.
##'
##' @name rprocess
##' @aliases rprocess rprocess,pomp-method rprocess-pomp
##' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##' @param xstart an \code{nvar} x \code{nrep} matrix containing the starting state of the system.
##' Columns of \code{xstart} correspond to states;
##' rows to components of the state vector.
##' One independent simulation will be performed for each column.
##' Note that in this case, \code{params} must also have \code{nrep} columns.
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
##' @export
setMethod(
  "rprocess",
  signature=signature(object="pomp"),
  definition=function (object, xstart, times, params, offset = 0L, ...) {
    tryCatch(
      rprocess.internal(object=object,xstart=xstart,times=times,params=params,offset=offset,...),
      error = function (e) pomp_stop("rprocess",conditionMessage(e))
    )
  }
)

rprocess.internal <- function (object, xstart, times, params, offset = 0L, .getnativesymbolinfo = TRUE, ...) {
  storage.mode(xstart) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(do_rprocess,object,xstart,times,params,offset,.getnativesymbolinfo)
}

##' Generic skeleton
##'
##' @name skeleton-generic
##' @aliases skeleton,missing-method skeleton,ANY-method
##' @keywords internal
setGeneric(
  "skeleton",
  function (object, ...)
    standardGeneric("skeleton")
)

setMethod(
  "skeleton",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("skeleton"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "skeleton",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("skeleton")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' skeleton
##'
##' Evaluates the deterministic skeleton at a point or points in state space, given parameters.
##' In the case of a discrete-time system, the skeleton is a map.
##' In the case of a continuous-time system, the skeleton is a vectorfield.
##' NB: \code{skeleton} just evaluates the deterministic skeleton;
##' it does not iterate or integrate (see \code{\link{trajectory}} for this).
##'
##' @name skeleton
##' @aliases skeleton skeleton,pomp-method skeleton-pomp
##' @keywords internal
##' @family pomp workhorses
##'
##' @inheritParams dmeasure
##'
##' @return
##' \code{skeleton} returns an array of dimensions \code{nvar} x \code{nrep} x \code{ntimes}.
##' If \code{f} is the returned matrix, \code{f[i,j,k]} is the i-th component of the deterministic skeleton at time \code{times[k]} given the state \code{x[,j,k]} and parameters \code{params[,j]}.
##'
##' @export
setMethod(
  "skeleton",
  signature=signature("pomp"),
  definition=function (object, x, times, params, ...)
    tryCatch(
      skeleton.internal(object=object,x=x,times=times,params=params,...),
      error = function (e) pomp_stop("skeleton",conditionMessage(e))
    )
)

skeleton.internal <- function (object, x, times, params, .getnativesymbolinfo = TRUE, ...) {
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(do_skeleton,object,x,times,params,.getnativesymbolinfo)
}
