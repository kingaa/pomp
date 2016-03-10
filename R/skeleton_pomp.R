## evaluate the deterministic skeleton

skeleton.internal <- function (object, x, t, params, .getnativesymbolinfo = TRUE, ...) {
  pompLoad(object)
  rv <- .Call(do_skeleton,object,x,t,params,.getnativesymbolinfo)
  pompUnload(object)
  rv
}

setMethod("skeleton",
          signature=signature("pomp"),
          definition=function (object, x, t, params, ...)
          skeleton.internal(object=object,x=x,t=t,params=params,...)
          )
