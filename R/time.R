## methods to manipulate the vector of times

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

setMethod(
  "time",
  signature=signature(x="pomp"),
  definition=function (x, t0 = FALSE, ...) {
    if (t0) c(x@t0,x@times) else x@times
  }
)

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
