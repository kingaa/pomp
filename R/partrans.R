## parameter transformations

setMethod(
  "partrans",
  signature=signature(object="pomp"),
  definition=function (object, params,
    dir = c("fromEst", "toEst"), ...) {
    dir <- match.arg(dir)
    partrans.internal(object=object,params=params,dir=dir,...)
  }
)

partrans.internal <- function (object, params, dir = c("fromEst", "toEst"),
  .getnativesymbolinfo = TRUE, ...) {

  if (object@partrans@has) {
    dir <- switch(dir,fromEst=-1L,toEst=1L)
    pompLoad(object)
    on.exit(pompUnload(object))
    params <- .Call(do_partrans,object,params,dir,.getnativesymbolinfo)
  }
  params
}
