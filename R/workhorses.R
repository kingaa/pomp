##' Workhorse functions for the \pkg{pomp} algorithms.
##'
##' These functions mediate the interface between the user's model and the package algorithms.
##' They are low-level functions that do the work needed by the package's inference methods.
##'
##' They include \describe{
##' \item{\code{\link{rinit}}}{which samples from the initial-state distribution,}
##' \item{\code{\link{dinit}}}{which evaluates the initial-state density,}
##' \item{\code{\link{dmeasure}}}{which evaluates the measurement model density,}
##' \item{\code{\link{rmeasure}}}{which samples from the measurement model distribution,}
##' \item{\code{\link{emeasure}}}{which computes the expectation of the observed variables conditional on the latent state,}
##' \item{\code{\link{vmeasure}}}{which computes the covariance matrix of the observed variables conditional on the latent state,}
##' \item{\code{\link{dprocess}}}{which evaluates the process model density,}
##' \item{\code{\link{rprocess}}}{which samples from the process model distribution,}
##' \item{\code{\link{dprior}}}{which evaluates the prior probability density,}
##' \item{\code{\link{rprior}}}{which samples from the prior distribution,}
##' \item{\code{\link{skeleton}}}{which evaluates the model's deterministic skeleton,}
##' \item{\code{\link{flow}}}{which iterates or integrates the deterministic skeleton to yield trajectories,}
##' \item{\code{\link{partrans}}}{which performs parameter transformations associated with the model.}
##' }
##'
##' @name workhorses
##' @include pomp_class.R pomp_fun.R load.R pstop.R
##' @docType methods
##' @family pomp workhorses
##' @concept extending the pomp package
##' @concept low-level interface
##' @seealso \link[=basic components]{basic model components},
##' \link[=elementary algorithms]{elementary algorithms},
##' \link[=estimation algorithms]{estimation algorithms}
##' @author Aaron A. King
##' @concept extending the pomp package
##' @concept low-level interface
NULL

##' dmeasure
##'
##' \code{dmeasure} evaluates the probability density of observations given states.
##'
##' @name dmeasure
##' @docType methods
##' @aliases dmeasure,ANY-method dmeasure,missing-method
##' @family pomp workhorses
##' @concept extending the pomp package
##' @concept low-level interface
##' @seealso Specification of the measurement density evaluator: \link{dmeasure specification}
##' @param object an object of class \sQuote{pomp}, or of a class that extends \sQuote{pomp}.
##' This will typically be the output of \code{pomp}, \code{simulate}, or one of the \pkg{pomp} inference algorithms.
##' @param x an array containing states of the unobserved process.
##' The dimensions of \code{x} are \code{nvars} x \code{nrep} x \code{ntimes},
##' where \code{nvars} is the number of state variables,
##' \code{nrep} is the number of replicates,
##' and \code{ntimes} is the length of \code{times}.
##' One can also pass \code{x} as a named numeric vector, which is equivalent to the \code{nrep=1}, \code{ntimes=1} case.
##' @param y a matrix containing observations.
##' The dimensions of \code{y} are \code{nobs} x \code{ntimes}, where \code{nobs} is the number of observables
##' and \code{ntimes} is the length of \code{times}.
##' @param times a numeric vector (length \code{ntimes}) containing times.
##' These must be in non-decreasing order.
##' @param params a \code{npar} x \code{nrep} matrix of parameters.
##' Each column is treated as an independent parameter set, in correspondence with the corresponding column of \code{x}.
##' @param log if TRUE, log probabilities are returned.
##' @param \dots additional arguments are ignored.
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
  definition=function (
    object,
    y = obs(object),
    x = states(object),
    times = time(object),
    params = coef(object),
    ...,
    log = FALSE
  ) {
    tryCatch(
      dmeasure_internal(object=object,y=y,x=x,times=times,
        params=params,log=log,...),
      error = function (e) pStop(who="dmeasure",conditionMessage(e))
    )
  }
)

dmeasure_internal <- function (object, y, x, times, params, ..., log = FALSE,
  .gnsi = TRUE) {
  storage.mode(y) <- "double"
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_dmeasure,object,y,x,times,params,log,.gnsi)  
}

##' dprior
##'
##' Evaluates the prior probability density.
##'
##' @name dprior
##' @docType methods
##' @aliases dprior,ANY-method dprior,missing-method
##' @family pomp workhorses
##' @family Bayesian methods
##' @concept extending the pomp package
##' @concept low-level interface
##' @seealso Specification of the prior density evaluator: \link{prior specification}
##' @inheritParams dmeasure
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
  definition=function (
    object,
    params = coef(object),
    ...,
    log = FALSE
  ) {
    tryCatch(
      dprior_internal(object=object,params=params,log=log,...),
      error = function (e) pStop(who="dprior",conditionMessage(e))
    )
  }
)

dprior_internal <- function (object, params, log = FALSE,
  .gnsi = TRUE, ...) {
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_dprior,object,params,log,.gnsi)
}

##' dprocess
##'
##' Evaluates the probability density of a sequence of consecutive state transitions.
##'
##' @name dprocess
##' @docType methods
##' @aliases dprocess,ANY-method dprocess,missing-method
##' @family pomp workhorses
##' @concept extending the pomp package
##' @concept low-level interface
##' @seealso Specification of the process-model density evaluator: \link{dprocess specification}
##' @inheritParams dmeasure
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
  definition = function (
    object,
    x = states(object),
    times = time(object),
    params = coef(object),
    ...,
    log = FALSE
  ) {
    tryCatch(
      dprocess_internal(object=object,x=x,times=times,params=params,log=log,...),
      error = function (e) pStop(who="dprocess",conditionMessage(e))
    )
  }
)

dprocess_internal <- function (object, x, times, params, log = FALSE, .gnsi = TRUE, ...) {
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_dprocess,object,x,times,params,log,.gnsi)
}

##' partrans
##'
##' Performs parameter transformations.
##'
##' @name partrans
##' @docType methods
##' @aliases partrans,ANY-method partrans,missing-method
##' @family pomp workhorses
##' @concept extending the pomp package
##' @concept low-level interface
##' @seealso Specification of parameter transformations: \code{\link{parameter_trans}}
##' @inheritParams dmeasure
##' @param dir the direction of the transformation to perform.
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
  definition=function (
    object, params,
    dir = c("fromEst", "toEst"),
    ...
  ) {
    dir <- match.arg(dir)
    tryCatch(
      partrans_internal(object=object,params=params,dir=dir,...),
      error = function (e) pStop(who="partrans",conditionMessage(e))
    )
  }
)

partrans_internal <- function (object, params, dir = c("fromEst", "toEst"),
  .gnsi = TRUE, ...) {
  if (object@partrans@has) {
    dir <- switch(dir,fromEst=-1L,toEst=1L)
    pompLoad(object)
    on.exit(pompUnload(object))
    params <- .Call(P_do_partrans,object,params,dir,.gnsi)
  }
  params
}

##' rinit
##'
##' Samples from the initial-state distribution.
##'
##' @name rinit
##' @docType methods
##' @aliases rinit,ANY-method rinit,missing-method
##' @family pomp workhorses
##' @concept extending the pomp package
##' @concept low-level interface
##' @seealso Specification of the initial-state distribution: \link{rinit specification}
##' @inheritParams dmeasure
##' @param t0 the initial time, i.e., the time corresponding to the initial-state distribution.
##' @param nsim optional integer; the number of initial states to simulate per column of \code{params}.
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
  definition=function (
    object,
    params = coef(object),
    t0 = timezero(object),
    nsim = 1,
    ...
  ) {
    tryCatch(
      rinit_internal(object=object,params=params,t0=t0,nsim=nsim,...),
      error = function (e) pStop(who="rinit",conditionMessage(e))
    )
  }
)

rinit_internal <- function (object, params, t0, nsim = 1,
  .gnsi = TRUE, ...) {
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_rinit,object,params,t0,nsim,.gnsi)
}

##' dinit
##'
##' Evaluates the initial-state density.
##'
##' @name dinit
##' @docType methods
##' @aliases dinit,ANY-method dinit,missing-method
##' @family pomp workhorses
##' @concept extending the pomp package
##' @concept low-level interface
##' @seealso Specification of the initial-state distribution: \link{dinit specification}
##' @inheritParams dmeasure
##' @param t0 the initial time, i.e., the time corresponding to the initial-state distribution.
##' @return
##' \code{dinit} returns a 1-D numerical array containing the likelihoods (or log likelihoods if \code{log=TRUE}).
##' By default, \code{t0} is the initial time defined when the \sQuote{pomp} object ws constructed.
##'
NULL

setGeneric(
  "dinit",
  function (object, ...)
    standardGeneric("dinit")
)

setMethod(
  "dinit",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("dinit","object")
  }
)

setMethod(
  "dinit",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("dinit",object)
  }
)

##' @export
##' @rdname dinit
setMethod(
  "dinit",
  signature=signature("pomp"),
  definition=function (
    object,
    params = coef(object),
    t0 = timezero(object),
    x,
    log = FALSE,
    ...
  ) {
    tryCatch(
      dinit_internal(object=object,x=x,params=params,t0=t0,log=log,...),
      error = function (e) pStop(who="dinit",conditionMessage(e))
    )
  }
)

dinit_internal <- function (object, x, params, t0, log, .gnsi = TRUE, ...) {
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_dinit,object,t0,x,params,log,.gnsi)
}

##' rmeasure
##'
##' Sample from the measurement model distribution, given values of the latent states and the parameters.
##'
##' @name rmeasure
##' @docType methods
##' @aliases rmeasure,ANY-method rmeasure,missing-method
##' @family pomp workhorses
##' @concept extending the pomp package
##' @concept low-level interface
##' @seealso Specification of the measurement-model simulator: \link{rmeasure specification}
##' @inheritParams dmeasure
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
  definition=function (
    object,
    x = states(object),
    times = time(object),
    params = coef(object),
    ...
  ) {
    tryCatch(
      rmeasure_internal(object=object,x=x,times=times,params=params,...),
      error = function (e) pStop(who="rmeasure",conditionMessage(e))
    )
  }
)

rmeasure_internal <- function (object, x, times, params,
  .gnsi = TRUE, ...) {
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_rmeasure,object,x,times,params,.gnsi)
}

##' emeasure
##'
##' Return the expected value of the observed variables, given values of the latent states and the parameters.
##'
##' @name emeasure
##' @docType methods
##' @aliases emeasure,ANY-method emeasure,missing-method
##' @family pomp workhorses
##' @seealso Specification of the measurement-model expectation: \link{emeasure specification}
##' @inheritParams dmeasure
##' @return
##' \code{emeasure} returns a rank-3 array of dimensions
##' \code{nobs} x \code{nrep} x \code{ntimes},
##' where \code{nobs} is the number of observed variables.
##'
NULL

setGeneric(
  "emeasure",
  function (object, ...)
    standardGeneric("emeasure")
)

setMethod(
  "emeasure",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("emeasure","object")
  }
)

setMethod(
  "emeasure",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("emeasure",object)
  }
)

##' @export
##' @rdname emeasure
setMethod(
  "emeasure",
  signature=signature(object="pomp"),
  definition=function (
    object,
    x = states(object),
    times = time(object),
    params = coef(object),
    ...
  ) {
    tryCatch(
      emeasure_internal(object=object,x=x,times=times,params=params,...),
      error = function (e) pStop(who="emeasure",conditionMessage(e))
    )
  }
)

emeasure_internal <- function (object, x, times, params,
  .gnsi = TRUE, ...) {
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_emeasure,object,x,times,params,.gnsi)
}

##' vmeasure
##'
##' Return the covariance matrix of the observed variables, given values of the latent states and the parameters.
##'
##' @name vmeasure
##' @docType methods
##' @aliases vmeasure,ANY-method vmeasure,missing-method
##' @family pomp workhorses
##' @seealso Specification of the measurement-model covariance matrix: \link{vmeasure specification}
##' @inheritParams dmeasure
##' @return
##' \code{vmeasure} returns a rank-4 array of dimensions
##' \code{nobs} x \code{nobs} x \code{nrep} x \code{ntimes},
##' where \code{nobs} is the number of observed variables.
##' If \code{v} is the returned array, \code{v[,,j,k]} contains the
##' covariance matrix at time \code{times[k]} given the state \code{x[,j,k]}.
##'
NULL

setGeneric(
  "vmeasure",
  function (object, ...)
    standardGeneric("vmeasure")
)

setMethod(
  "vmeasure",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("vmeasure","object")
  }
)

setMethod(
  "vmeasure",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("vmeasure",object)
  }
)

##' @export
##' @rdname vmeasure
setMethod(
  "vmeasure",
  signature=signature(object="pomp"),
  definition=function (
    object,
    x = states(object),
    times = time(object),
    params = coef(object),
    ...
  ) {
    tryCatch(
      vmeasure_internal(object=object,x=x,times=times,params=params,...),
      error = function (e) pStop(who="vmeasure",conditionMessage(e))
    )
  }
)

vmeasure_internal <- function (object, x, times, params,
  .gnsi = TRUE, ...) {
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_vmeasure,object,x,times,params,.gnsi)
}

##' rprior
##'
##' Sample from the prior probability distribution.
##'
##' @name rprior
##' @docType methods
##' @aliases rprior,ANY-method rprior,missing-method
##' @family pomp workhorses
##' @family Bayesian methods
##' @concept extending the pomp package
##' @concept low-level interface
##' @seealso Specification of the prior distribution simulator: \link{prior specification}
##' @inheritParams dmeasure
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
  definition=function (
    object,
    params = coef(object),
    ...
  )
    tryCatch(
      rprior_internal(object=object,params=params,...),
      error = function (e) pStop(who="rprior",conditionMessage(e))
    )
)

rprior_internal <- function (object, params, .gnsi = TRUE, ...) {
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_rprior,object,params,.gnsi)
}

##' rprocess
##'
##' \code{rprocess} simulates the process-model portion of partially-observed Markov process.
##'
##' When \code{rprocess} is called, \code{t0} is taken to be the initial time (i.e., that corresponding to \code{x0}).
##' The values in \code{times} are the times at which the state of the simulated processes are required.
##'
##' @name rprocess
##' @docType methods
##' @aliases rprocess,ANY-method rprocess,missing-method
##' @family pomp workhorses
##' @concept extending the pomp package
##' @concept low-level interface
##' @seealso Specification of the process-model simulator: \link{rprocess specification}
##' @inheritParams dmeasure
##' @param x0 an \code{nvar} x \code{nrep} matrix containing the starting state of the system.
##' Columns of \code{x0} correspond to states;
##' rows to components of the state vector.
##' One independent simulation will be performed for each column.
##' Note that in this case, \code{params} must also have \code{nrep} columns.
##' @param t0 the initial time, i.e., the time corresponding to the state in \code{x0}.
##' @param params a \code{npar} x \code{nrep} matrix of parameters.
##' Each column is treated as an independent parameter set, in correspondence with the corresponding column of \code{x0}.
##' @return
##' \code{rprocess} returns a rank-3 array with rownames.
##' Suppose \code{x} is the array returned.
##' Then \preformatted{dim(x)=c(nvars,nrep,ntimes),}
##' where \code{nvars} is the number of state variables (=\code{nrow(x0)}),
##' \code{nrep} is the number of independent realizations simulated (=\code{ncol(x0)}), and
##' \code{ntimes} is the length of the vector \code{times}.
##' \code{x[,j,k]} is the value of the state process in the \code{j}-th realization at time \code{times[k]}.
##' The rownames of \code{x} will correspond to those of \code{x0}.
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
  definition=function (
    object,
    x0 = rinit(object),
    t0 = timezero(object),
    times = time(object),
    params = coef(object),
    ...
  ) {
    tryCatch(
      rprocess_internal(object=object,x0=x0,t0=t0,times=times,params=params,...),
      error = function (e) pStop(who="rprocess",conditionMessage(e))
    )
  }
)

rprocess_internal <- function (object, x0, t0, times, params, ...,
  .gnsi = TRUE) {
  storage.mode(x0) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_rprocess,object,x0,t0,times,params,.gnsi)
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
##' @aliases skeleton,ANY-method skeleton,missing-method
##' @family pomp workhorses
##' @family deterministic methods
##' @concept extending the pomp package
##' @concept low-level interface
##' @seealso Specification of the deterministic skeleton: \link{skeleton specification}
##' @inheritParams dmeasure
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
  definition=function (
    object,
    x = states(object),
    times = time(object),
    params = coef(object),
    ...
  )
    tryCatch(
      skeleton_internal(object=object,x=x,times=times,params=params,...),
      error = function (e) pStop(who="skeleton",conditionMessage(e))
    )
)

skeleton_internal <- function (object, x, times, params, .gnsi = TRUE, ...) {
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  pompLoad(object)
  on.exit(pompUnload(object))
  .Call(P_do_skeleton,object,x,times,params,.gnsi)
}
