## extract coefficients (estimated parameters)

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
  "coef",
  signature=signature(object="listies"),
  definition=function(object, ...) {
    do.call(cbind,lapply(object,coef))
  }
)

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

## modify the coefficients
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
