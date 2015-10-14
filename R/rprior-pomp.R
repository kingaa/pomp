## simulate from the prior

rprior.internal <- function (object, params, .getnativesymbolinfo = TRUE, ...) {
    pompLoad(object)
    rv <- .Call(do_rprior,object,params,.getnativesymbolinfo)
    pompUnload(object)
    rv
}

setMethod(
    "rprior",
    signature=signature(object="pomp"),
    definition=function (object, params, ...)
        rprior.internal(object=object,params=params,...)
)
