## simulate from the prior

rprior.internal <- function (object, params, .getnativesymbolinfo = TRUE, ...) {
    storage.mode(params) <- "double"
    pompLoad(object)
    on.exit(pompUnload(object))
    rv <- .Call(do_rprior,object,params,.getnativesymbolinfo)
    rv
}

setMethod(
    "rprior",
    signature=signature(object="pomp"),
    definition=function (object, params, ...)
        rprior.internal(object=object,params=params,...)
)
