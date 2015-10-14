## evaluate the prior probability density

dprior.internal <- function (object, params, log = FALSE,
                             .getnativesymbolinfo = TRUE, ...) {
    pompLoad(object)
    rv <- .Call(do_dprior,object,params,log,.getnativesymbolinfo)
    pompUnload(object)
    rv
}

setMethod(
    "dprior",
    signature=signature(object="pomp"),
    definition=function (object, params, log = FALSE, ...)
        dprior.internal(object=object,params=params,log=log,...)
)
