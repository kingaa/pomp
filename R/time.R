##' Methods to manipulate the obseration times
##'
##' Get and set the vector of observation times.
##'
##' @name time
##' @rdname time
##' @aliases time time<- time,ANY-method time,missing-method
NULL

setGeneric(
  "time",
  function (x, ...)
    standardGeneric("time")
)

setGeneric(
  "time<-",
  function (object, ..., value)
    standardGeneric("time<-")
)

##' @name time-pomp
##' @aliases time,pomp-method
##' @rdname time
##' @param x  a \sQuote{pomp} object
##' @param t0 logical; should the zero time be included?
##' @param \dots ignored
##'
##' @details
##' \code{time(object)} returns the vector of observation times.
##' \code{time(object,t0=TRUE)} returns the vector of observation
##' times with the zero-time \code{t0} prepended.

setMethod(
  "time",
  signature=signature(x="pomp"),
  definition=function (x, t0 = FALSE, ...) {
    if (t0) c(x@t0,x@times) else x@times
  }
)

##' @name time<--pomp
##' @aliases time<-,pomp-method
##' @rdname time
##' @param object  a \sQuote{pomp} object
##' @param value numeric vector; the new vector of times
##'
##' @details
##' \code{time(object) <- value} replaces the observation times slot (\code{times}) of \code{object} with \code{value}.
##' \code{time(object,t0=TRUE) <- value} has the same effect, but the first element in \code{value} is taken to be the initial time.
##' The second and subsequent elements of \code{value} are taken to be the observation times.
##' Those data and states (if they exist) corresponding to the new times are retained.
setMethod(
  "time<-",
  signature=signature(object="pomp"),
  definition=function (object, t0 = FALSE, ..., value) {
    ep <- paste0("in ",sQuote("time<-"),": ")
    if (!is.numeric(value))
      stop(ep,sQuote("value")," must be a numeric vector.",call.=FALSE)
    storage.mode(value) <- "double"
    tt <- object@times
    dd <- object@data
    ss <- object@states
    if (t0) {
      object@t0 <- value[1]
      object@times <- value[-1]
    } else {
      object@times <- value
    }
    if (!all(diff(object@times)>0))
      stop(ep,"the times specified must be an increasing sequence.",call.=FALSE)
    if (object@t0>object@times[1])
      stop(ep,"the zero-time ",sQuote("t0")," must occur no later than the first observation.",call.=FALSE)
    object@data <- array(
      data=NA,
      dim=c(nrow(dd),length(object@times)),
      dimnames=list(rownames(dd),NULL)
    )
    object@data[,object@times%in%tt] <- dd[,tt%in%object@times]
    if (length(ss)>0) {
      object@states <- array(
        data=NA,
        dim=c(nrow(ss),length(object@times)),
        dimnames=list(rownames(ss),NULL)
      )
      for (kt in seq_along(object@times)) {
        wr <- which(object@times[kt]==tt)
        if (length(wr)>0)
          object@states[,kt] <- ss[,wr[1]]
      }
    }
    object
  }
)
