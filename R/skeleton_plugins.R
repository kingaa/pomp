setClass(
  "skelPlugin",
  slots=c(
    csnippet='logical',
    slotname='character',
    type='integer',
    skel.fn="ANY"
  ),
  prototype=prototype(
    csnippet=FALSE,
    slotname=character(0),
    type=0L,
    skel.fn=NULL
  )
)

setClass(
  "vectorfieldPlugin",
  contains="skelPlugin"
)

setClass(
  "mapPlugin",
  contains="skelPlugin",
  slots=c(
    delta.t="numeric"
  )
)

skel_plugin <- function (object, skel.fn) {
  if (missing(object)) {
    new("skelPlugin")
  } else {
    if (!missing(skel.fn)) object@skel.fn <- skel.fn
    object
  }
}

vectorfield <- function (f) {
  new("vectorfieldPlugin",skel.fn=f,type=1L)
}

map <- function (f, delta.t = 1) {
  if (!isTRUE(delta.t > 0 && length(delta.t)==1))
    stop("in ",sQuote("map"),": ",sQuote("delta.t"),
      " must be a positive number.",call.=FALSE)
  new("mapPlugin",skel.fn=f,delta.t=delta.t,type=2L)
}
