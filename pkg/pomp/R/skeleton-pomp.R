## evaluate the deterministic skeleton
setGeneric("skeleton",function(object,...)standardGeneric("skeleton"))

skeleton.internal <- function (object, x, t, params, ...) {
  skel <- .Call(get_pomp_fun,object@skeleton)
  .Call(do_skeleton,object,x,t,params,skel)
}

setMethod("skeleton","pomp",skeleton.internal)
