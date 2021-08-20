##' Methods to manipulate the obseration times
##'
##' Get and set the vector of observation times.
##'
##' @name time
##' @rdname time
##' @aliases time<- time,missing-method
##'
##' @importFrom stats time
NULL

setGeneric("time")

setGeneric(
  "time<-",
  function (object, ..., value)
    standardGeneric("time<-")
)

##' @rdname time
##' @param x  a \sQuote{pomp} object
##' @param t0 logical; should the zero time be included?
##' @param \dots ignored or passed to the more primitive function
##' @details
##' \code{time(object)} returns the vector of observation times.
##' \code{time(object,t0=TRUE)} returns the vector of observation
##' times with the zero-time \code{t0} prepended.
##' @export
setMethod(
  "time",
  signature=signature(x="pomp"),
  definition=function (x, t0 = FALSE, ...) {
    if (t0) c(x@t0,x@times) else x@times
  }
)

##' @rdname time
##' @param object  a \sQuote{pomp} object
##' @param value numeric vector; the new vector of times
##' @details
##' \code{time(object) <- value} replaces the observation times slot (\code{times}) of \code{object} with \code{value}.
##' \code{time(object,t0=TRUE) <- value} has the same effect, but the first element in \code{value} is taken to be the initial time.
##' The second and subsequent elements of \code{value} are taken to be the observation times.
##' Those data and states (if they exist) corresponding to the new times are retained.
##' @export
setMethod(
  "time<-",
  signature=signature(object="pomp"),
  definition=function (object, t0 = FALSE, ..., value) {
    tryCatch(
      time.repl.internal(object,t0=t0,...,value=value),
      error = function (e) pStop("time<-",conditionMessage(e))
    )
  }
)

time.repl.internal <- function (object, t0 = FALSE, ..., value) {
  if (!is.numeric(value)) pStop_(sQuote("value")," must be a numeric vector.")
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
  if (!all(diff(object@times)>=0))
    pStop_("times must be a non-decreasing numeric sequence.")
  if (object@t0>object@times[1])
    pStop_("the zero-time ",sQuote("t0"),
      " must occur no later than the first observation.")
  object@data <- array(
    data=NA,
    dim=c(nrow(dd),length(object@times)),
    dimnames=list(rownames(dd),NULL)
  )
  idx2 <- match(object@times,tt,nomatch=0)
  idx1 <- idx2 > 0
  object@data[,idx1] <- dd[,idx2]
  if (length(ss)>0) {
    object@states <- array(
      data=NA,
      dim=c(nrow(ss),length(object@times)),
      dimnames=list(rownames(ss),NULL)
    )
    object@states[,idx1] <- ss[,idx2]
  }
  object
}
