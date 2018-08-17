##' Spectrum matching
##'
##' Estimation of parameters by matching power spectra
##'
##' In spectrum matching, one attempts to minimize the discrepancy between a \acronym{POMP} model's predictions and data, as measured in the frequency domain by the power spectrum.
##'
##' \code{spect.match.objfun} constructs an objective function that measures the discrepancy.
##' It can be passed to any one of a variety of numerical optimization routines, which will adjust model parameters to minimize the discrepancies between the power spectrum of model simulations and that of the data.
##'
##' @name spect.match
##' @docType methods
##' @rdname spect_match
##' @include spect.R probe_match.R loglik.R summary.R
##' @aliases spect.match.objfun,missing-method spect.match.objfun,ANY-method
##'
##' @return
##' \code{spect.match.objfun} construct a stateful objective function for spectrum matching.
##' Specfically, \code{spect.match.objfun} returns an object of class \sQuote{spect_match_objfun}, which is a function suitable for use in an \code{\link{optim}}-like optimizer.
##' In particular, this function takes a single numeric-vector argument that is assumed to contain the parameters named in \code{est}, in that order.
##' When called, it will return the (optionally weighted) \eqn{L^2}{L2} distance between the data spectrum and simulated spectra.
##' It is a stateful function:
##' Each time it is called, it will remember the values of the parameters and its estimate of the synthetic likelihood.
##'
##' @section Important Note:
##' Since \pkg{pomp} cannot guarantee that the \emph{final} call an optimizer makes to the function is a call \emph{at} the optimum, it cannot guarantee that the parameters stored in the function are the optimal ones.
##' For this reason, \code{\link[=logLik,spect_match_objfun-method]{logLik}} and \code{\link[=summary,spect_match_objfun-method]{summary}} both call \code{spect} on the estimated parameters.
##' One should check that the parameters agree with those that are returned by the optimizer.
##' The best practice is to call \code{\link[=spect,spect_match_objfun-method]{spect}} on the objective function after the optimization has been performed, thus obtaining a \sQuote{spectd_pomp} object containing the (putative) optimal parameters.
##' 
##' @seealso \code{\link{trajectory}}, \code{\link{optim}},
##' \code{\link[subplex]{subplex}}, \code{\link[nloptr]{nloptr}}
##' 
NULL

setClass(
  "spect_match_objfun",
  contains="function",
  slots=c(
    env="environment"
  )
)

setGeneric(
  "spect.match.objfun",
  function (object, ...)
    standardGeneric("spect.match.objfun")
)

setMethod(
  "spect.match.objfun",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("spect.match.objfun"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "spect.match.objfun",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("spect.match.objfun")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' @name spect.match.objfun-pomp
##' @aliases spect.match.objfun spect.match.objfun,pomp-method
##' @rdname spect_match
##'
##' @inheritParams probe.match
##' @inheritParams spect
##' @param weights optional numeric or function.
##' The mismatch between model and data is measured by a weighted average of mismatch at each frequency.
##' By default, all frequencies are weighted equally.
##' \code{weights} can be specified either as a vector (which must have length equal to the number of frequencies) or as a function of frequency.
##' If the latter, \code{weights(freq)} must return a nonnegative weight for each frequency.
##'
setMethod(
  "spect.match.objfun",
  signature=signature(object="pomp"),
  definition=function(object, params, est = character(0),
    vars, nsim, seed = NULL, kernel.width, transform.data,
    detrend = c("none","mean","linear","quadratic"),
    weights = 1, fail.value = NA, transform = FALSE, ...) {

    transform.data <- match.fun(transform.data)
    detrend <- match.arg(detrend)

    smof.internal(object,params=params,est=est,vars=vars,nsim=nsim,
      seed=seed,kernel.width=kernel.width,transform.data=transform.data,
      detrend=detrend,weights=weights,fail.value=fail.value,
      transform=transform,...)

  }
)

##' @name spect.match.objfun-spectd_pomp
##' @aliases spect.match.objfun,spectd_pomp-method
##' @rdname spect_match
setMethod(
  "spect.match.objfun",
  signature=signature(object="spectd_pomp"),
  definition=function(object, params, est, vars, nsim, seed = NULL,
    kernel.width, transform.data, detrend, ...) {

    if (missing(params)) params <- object@params
    if (missing(vars)) vars <- object@vars
    if (missing(nsim)) nsim <- object@nsim
    if (missing(kernel.width)) kernel.width <- object@kernel.width
    if (missing(transform.data)) transform.data <- object@transform.data
    if (missing(detrend)) detrend <- object@detrend

    spect.match.objfun(object,params=params,est=est,vars=vars,nsim=nsim,
      seed=seed,kernel.width=kernel.width, transform.data=transform.data,
      detrend=detrend,...)

  }
)

##' @name spect.match.objfun-spect_match_objfun
##' @aliases spect.match.objfun,spect_match_objfun-method
##' @rdname spect_match
setMethod(
  "spect.match.objfun",
  signature=signature(object="spect_match_objfun"),
  definition=function(object, ...) {

    spect.match.objfun(object@env$object,...)

  }
)

smof.internal <- function (object, params, est, vars, nsim, seed,
  kernel.width, transform.data, detrend, weights, fail.value, transform, ...) {

  ep <- paste0("in ",sQuote("spect.match.objfun"),": ")

  transform <- as.logical(transform)
  fail.value <- as.numeric(fail.value)

  if (missing(params)) params <- coef(object)

  if (missing(est)) est <- character(0)
  est <- as.character(est)
  est <- est[nzchar(est)]

  if (missing(vars)) vars <- rownames(object@data)
  if (missing(nsim)) stop(ep,sQuote("nsim")," must be supplied",call.=FALSE)
  if (missing(kernel.width))
    stop(ep,sQuote("kernel.width")," must be specified",call.=FALSE)
  if (missing(transform.data)) transform.data <- identity

  object <- spect(object,params=params,vars=vars,nsim=nsim,seed=seed,
    kernel.width=kernel.width,transform.data=transform.data,detrend=detrend,
    ...)

  if (is.numeric(weights)) {
    if (length(weights)==1) {
      weights <- rep(weights,length(object@freq))
    } else if ((length(weights) != length(object@freq)))
      stop(ep,"if ",sQuote("weights"),
        " is provided as a vector, it must have length ",
        length(object@freq),call.=FALSE)
  } else if (is.function(weights)) {
    weights <- tryCatch(
      vapply(object@freq,weights,numeric(1)),
      error = function (e) {
        stop(ep,"problem with ",sQuote("weights")," function: ",
          conditionMessage(e),call.=FALSE)
      }
    )
  } else {
    stop(ep,sQuote("weights"),
      " must be specified as a vector or as a function",call.=FALSE)
  }

  if (any(!is.finite(weights) | weights<0))
    stop(ep,sQuote("weights")," should be nonnegative and finite",call.=FALSE)
  weights <- weights/mean(weights)

  params <- coef(object,transform=transform)

  idx <- match(est,names(params))
  if (any(is.na(idx))) {
    missing <- est[is.na(idx)]
    stop(ep,ngettext(length(missing),"parameter","parameters")," ",
      paste(sQuote(missing),collapse=","),
      " not found in ",sQuote("params"),call.=FALSE)
  }

  pompLoad(object)

  ker <- reuman.kernel(kernel.width)
  discrep <- spect.discrep(object,ker=ker,weights=weights)

  ofun <- function (par) {
    params[idx] <- par
    coef(object,transform=transform) <<- params
    discrep <<- spect.discrep(object,ker=ker,weights=weights)
    if (is.finite(discrep) || is.na(fail.value)) discrep else fail.value
  }

  environment(ofun) <- list2env(
    list(object=object,transform=transform,fail.value=fail.value,
      params=params,idx=idx,discrep=discrep,seed=seed,ker=ker,weights=weights),
    parent=parent.frame(2)
  )

  new("spect_match_objfun",ofun,env=environment(ofun))
}

## compute a measure of the discrepancies between simulations and data
spect.discrep <- function (object, ker, weights) {

  ## estimate power spectra of simulations
  simvals <- compute.spect.sim(
    object,
    vars=object@vars,
    params=object@params,
    nsim=object@nsim,
    seed=object@seed,
    transform.data=object@transform.data,
    detrend=object@detrend,
    ker=ker
  )

  discrep <- array(dim=c(length(object@freq),length(object@vars)))
  sim.means <- colMeans(simvals)
  for (j in seq_along(object@freq)) {
    for (k in seq_along(object@vars)) {
      discrep[j,k] <- ((object@datspec[j,k]-sim.means[j,k])^2)/
        mean((simvals[,j,k]-sim.means[j,k])^2)
    }
    discrep[j,] <- weights[j]*discrep[j,]
  }

  sum(discrep)

}


##' @name spect-spect_match_objfun
##' @rdname spect
##' @aliases spect,spect_match_objfun-method
setMethod(
  "spect",
  signature=signature(object="spect_match_objfun"),
  definition=function (object, ...) {
    spect(object@env$object,...)
  }
)

##' @name summary-spect_match_objfun
##' @rdname summary
##' @aliases summary,spect_match_objfun-method
setMethod(
  "summary",
  signature=signature(object="spect_match_objfun"),
  definition=function (object) {
    summary(spect(object@env$object))
  }
)

##' @name logLik-spect_match_objfun
##' @rdname loglik
##' @aliases logLik,spect_match_objfun-method
setMethod(
  "logLik",
  signature=signature(object="spect_match_objfun"),
  definition=function (object) {
    discrep <- spect.discrep(object,ker=object@env$ker,
      weights=object@env$weights)
    -discrep
  }
)
