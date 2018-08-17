##' Extract, set, or alter coefficients
##'
##' Extract, set, or modify the estimated parameters from a fitted model.
##'
##' @name coef
##' @rdname coef
##' @docType methods
##' @aliases coef,missing-method coef,ANY-method
##' coef<-,missing-method coef<-,ANY-method
##' @include pomp_class.R listies.R
NULL

setGeneric(
  "coef",
  function (object, ...)
    standardGeneric("coef")
)

setGeneric(
  "coef<-",
  function (object, ..., value)
    standardGeneric("coef<-")
)

setMethod(
  "coef",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("coef"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "coef",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("coef")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "coef<-",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("coef<-"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "coef<-",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("coef<-")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' @name coef-listies
##' @aliases coef,listies-method
##' @rdname coef
setMethod(
  "coef",
  signature=signature(object="listies"),
  definition=function(object, ...) {
    do.call(cbind,lapply(object,coef))
  }
)

##' @name coef-pomp
##' @aliases coef coef,pomp-method
##' @rdname coef
##'
##' @param object an object of class \sQuote{pomp}, or of a class extending \sQuote{pomp}
##' @param pars optional character; names of parameters to be retrieved or set.
##' @param transform logical; perform parameter transformation?
##' @param \dots ignored
##'
##' @details
##' \code{coef} allows one to extract the parameters from a fitted model.
##'
##' \code{coef(object,transform=TRUE)} returns the parameters transformed onto
##' the estimation scale.
##'
setMethod(
  "coef",
  signature=signature(object="pomp"),
  definition=function (object, pars, transform = FALSE, ...) {
    if (length(object@params)>0) {
      if (transform)
        params <- partrans(object,params=object@params,dir="toEst")
      else
        params <- object@params
      if (missing(pars))
        pars <- names(params)
      else {
        excl <- setdiff(pars,names(params))
        if (length(excl)>0) {
          stop("in ",sQuote("coef"),": name(s) ",
            paste(sQuote(excl),collapse=","),
            " correspond to no parameter(s).",
            call.=FALSE)
        }
      }
      params[pars]
    } else {
      numeric(0)
    }
  }
)

##' @name coef<--pomp
##' @aliases coef<- coef<-,pomp-method
##' @rdname coef
##'
##' @param value numeric vector or list; values to be assigned.
##' If \code{value = NULL}, the parameters are unset.
##'
##' @details
##' \code{coef(object) <- value} sets or alters the coefficients of a
##' \sQuote{pomp} object.
##'
##' \code{coef(object,transform=TRUE) <- value} assumes that \code{value} is on
##' the estimation scale, and applies the \dQuote{from estimation scale}
##' parameter transformation from \code{object} before altering the
##' coefficients.
##'
setMethod(
  "coef<-",
  signature=signature(object="pomp"),
  definition=function (object, pars, transform = FALSE, ..., value) {
    ep <- paste0("in ",sQuote("coef<-"),": ")
    if (is.null(value)) value <- numeric(0)
    if (is.list(value)) value <- unlist(value)
    if (missing(pars)) {          ## replace the whole params slot with 'value'
      if (length(value)>0) {
        if (transform)
          value <- partrans(object,params=value,dir="fromEst")
        pars <- names(value)
        if (is.null(pars)) {
          stop(ep,sQuote("value")," must be a named vector",call.=FALSE)
        }
      }
      object@params <- value
    } else { ## replace or append only the parameters named in 'pars'
      if (!is.null(names(value))) ## we ignore the names of 'value'
        warning(ep," names of ",sQuote("value")," are being discarded",call.=FALSE)
      if (length(object@params)==0) { ## no pre-existing 'params' slot
        val <- numeric(length(pars))
        names(val) <- pars
        val[] <- value
        if (transform) val <- partrans(object,params=val,dir="fromEst")
        object@params <- val
      } else { ## pre-existing params slot
        params <- coef(object,transform=transform)
        val <- numeric(length(pars))
        names(val) <- pars
        val[] <- value
        excl <- !(pars%in%names(params)) ## new parameter names
        if (any(excl)) { ## append parameters
          warning(ep,"name(s) ",
            paste(sQuote(pars[excl]),collapse=","),
            " do not refer to existing parameter(s);",
            " they are being concatenated",
            call.=FALSE)
          params <- c(params,val[excl])
        }
        params[pars] <- val
        if (transform)
          params <- partrans(object,params=params,dir="fromEst")
        object@params <- params
      }
    }
    storage.mode(object@params) <- "double"
    object
  }
)
