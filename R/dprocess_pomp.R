## evaluate the process model density function

dprocess.internal <- function (object, x, times, params, log = FALSE, .getnativesymbolinfo = TRUE, ...) {
    storage.mode(x) <- "double"
    storage.mode(params) <- "double"
    pompLoad(object)
    rv <- .Call(do_dprocess,object,x,times,params,log,.getnativesymbolinfo)
    pompUnload(object)
    rv
}

setMethod(
    "dprocess",
    signature=signature(object="pomp"),
    definition = function (object, x, times, params, log = FALSE, ...)
        dprocess.internal(object=object,x=x,times=times,params=params,log=log,...)
)
