## trajectory matching

setClass(
  "traj_matched_pomp",
  contains="pomp",
  slots=c(
    transform="logical",
    est="character",
    evals="integer",
    convergence="integer",
    msg="character",
    value="numeric"
  )
)

setGeneric("traj.match.objfun",
  function(object,...)standardGeneric("traj.match.objfun"))
setGeneric("traj.match",
  function(object,...)standardGeneric("traj.match"))

setMethod(
  "logLik",
  signature=signature(object="traj_matched_pomp"),
  definition=function(object, ...)
    object@value
)

setMethod(
  "traj.match.objfun",
  signature=signature(object="pomp"),
  function (object, params, est, transform = FALSE, ...) {

    tmof.internal(
      object=object,
      params=params,
      est=est,
      transform=transform,
      ...
    )

  }
)

setMethod(
  "traj.match",
  signature=signature(object="pomp"),
  function (object, start, est = character(0),
    method = c("Nelder-Mead","subplex","SANN","BFGS",
      "sannbox","nloptr"),
    transform = FALSE, ...)
  {

    if (missing(start)) start <- coef(object)
    if (is.list(start)) start <- unlist(start)

    method <- match.arg(method)
    est <- as.character(est)
    transform <- as.logical(transform)

    m <- minim.internal(
      objfun=traj.match.objfun(
        object=object,
        params=start,
        est=est,
        transform=transform
      ),
      start=start,
      est=est,
      object=object,
      method=method,
      transform=transform,
      ...
    )

    ## fill params slot appropriately
    coef(object) <- m$params

    ## fill states slot appropriately
    x <- trajectory(object)
    object@states <- array(data=x,dim=dim(x)[c(1L,3L)])
    rownames(object@states) <- rownames(x)

    new(
      "traj_matched_pomp",
      object,
      transform=transform,
      est=est,
      value=-m$value,
      evals=m$evals,
      convergence=m$convergence,
      msg=m$msg
    )
  }
)

setMethod(
  "traj.match",
  signature=signature(object="traj_matched_pomp"),
  function (object, est, transform, ...)
  {
    if (missing(est)) est <- object@est
    if (missing(transform)) transform <- object@transform

    traj.match(as(object,"pomp"),est=est,transform=transform,...)
  }
)

setMethod(
  "traj.match",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("traj.match"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "traj.match",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("traj.match")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "traj.match.objfun",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("traj.match.objfun"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "traj.match.objfun",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("traj.match.objfun")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

tmof.internal <- function (object, params, est, transform, ...) {

  ep <- paste0("in ",sQuote("traj.match.objfun"),": ")
  object <- as(object,"pomp")
  if (missing(est)) est <- character(0)
  est <- as.character(est)
  transform <- as.logical(transform)

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)
  if ((!is.numeric(params))||(is.null(names(params))))
    stop(ep,sQuote("params")," must be a named numeric vector",call.=FALSE)
  if (transform)
    params <- partrans(object,params,dir="toEst")
  par.est.idx <- match(est,names(params))
  if (any(is.na(par.est.idx)))
    stop(ep,"parameter(s): ",
      paste(sapply(est[is.na(par.est.idx)],sQuote),collapse=","),
      " not found in ",sQuote("params"),call.=FALSE)

  function (par) {
    pompLoad(object)
    on.exit(pompUnload(object))
    params[par.est.idx] <- par
    if (transform)
      tparams <- partrans(object,params,dir="fromEst")
    d <- dmeasure(
      object,
      y=object@data,
      x=trajectory(
        object,
        params=if (transform) tparams else params,
        ...
      ),
      times=time(object),
      params=if (transform) tparams else params,
      log=TRUE
    )
    -sum(d)
  }
}
