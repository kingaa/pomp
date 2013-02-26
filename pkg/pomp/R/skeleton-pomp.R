## evaluate the deterministic skeleton
setGeneric("skeleton",function(object,...)standardGeneric("skeleton"))

skeleton.internal <- function (object, x, t, params, .getnativesymbolinfo = TRUE, ...) {
  .Call(do_skeleton,object,x,t,params,.getnativesymbolinfo)
}

setMethod("skeleton","pomp",
          function (object, x, t, params, ...)
          skeleton.internal(object=object,x=x,t=t,params=params,...)
          )
