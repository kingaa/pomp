## evaluate the prior probability density

dprior.internal <- function (object, params, log = FALSE,
                             .getnativesymbolinfo = TRUE, ...) {
    storage.mode(params) <- "double"
    pompLoad(object)
    on.exit(pompUnload(object))
    rv <- .Call(do_dprior,object,params,log,.getnativesymbolinfo)
    rv
}

setMethod(
    "dprior",
    signature=signature(object="pomp"),
    definition=function (object, params, log = FALSE, ...)
        dprior.internal(object=object,params=params,log=log,...)
)
