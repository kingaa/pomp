## simulate from the prior

setGeneric("rprior",function(object,...)standardGeneric("rprior"))

rprior.internal <- function (object, params, .getnativesymbolinfo = TRUE, ...) {
  .Call(do_rprior,object,params,.getnativesymbolinfo)
}

setMethod("rprior","pomp",
          function (object, params, ...)
          rprior.internal(object=object,params=params,...)
          )
