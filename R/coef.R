##' Extract, set, or alter coefficients
##'
##' Extract, set, or modify the estimated parameters from a fitted model.
##'
##' @name coef
##' @rdname coef
##' @docType methods
##' @aliases coef,missing-method
##' coef<-,missing-method
##' @include pomp_class.R listie.R
##' @importFrom stats coef
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
  "coef<-",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("coef<-","object")
  }
)

##' @name coef-listie
##' @aliases coef,listie-method
##' @rdname coef
##' @export
setMethod(
  "coef",
  signature=signature(object="listie"),
  definition=function(object, ...) {
    x <- do.call(cbind,lapply(object,coef))
    names(dimnames(x)) <- c("parameter",".id")
    x
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
##' @export
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
        nexcl <- length(excl)
        if (nexcl > 0) {
          pStop("coef","name",ngettext(nexcl," ","s "),
            paste(sQuote(excl),collapse=","),
            " correspond",ngettext(nexcl,"s",""),
            " to no parameter",ngettext(nexcl,".","s."))
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
##' @export
setMethod(
  "coef<-",
  signature=signature(object="pomp"),
  definition=function (object, pars, transform = FALSE, ..., value) {
    if (is.null(value)) value <- numeric(0)
    if (is.list(value)) value <- unlist(value)
    if (missing(pars)) {          ## replace the whole params slot with 'value'
      if (length(value)>0) {
        if (transform)
          value <- partrans(object,params=value,dir="fromEst")
        pars <- names(value)
        if (is.null(pars)) {
          pStop("coef<-",sQuote("value")," must be a named vector.")
        }
      }
      object@params <- value
    } else { ## replace or append only the parameters named in 'pars'
      if (!is.null(names(value))) ## we ignore the names of 'value'
        pWarn("coef<-","names of ",sQuote("value")," are being discarded.")
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
          nexcl <- sum(excl)
          pWarn("coef<-","name",ngettext(nexcl,"","s")," ",
            paste(sQuote(pars[excl]),collapse=","),
            " refer",ngettext(nexcl,"s","")," to no existing parameter",
            ngettext(nexcl,"","s"),"; ",
            ngettext(nexcl,"it is","they are")," being concatenated.")
          params <- c(params,val[excl])
        }
        params[pars] <- val
        if (transform) params <- partrans(object,params=params,dir="fromEst")
        object@params <- params
      }
    }
    storage.mode(object@params) <- "double"
    object
  }
)
