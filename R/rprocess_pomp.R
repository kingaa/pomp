## simulate the process model

rprocess.internal <- function (object, xstart, times, params, offset = 0L, .getnativesymbolinfo = TRUE, ...) {
    storage.mode(xstart) <- "double"
    storage.mode(params) <- "double"
    pompLoad(object)
    rv <- .Call(do_rprocess,object,xstart,times,params,offset,.getnativesymbolinfo)
    pompUnload(object)
    rv
}

setMethod(
    "rprocess",
    signature=signature(object="pomp"),
    definition=function (object, xstart, times, params, offset = 0L, ...) {
        rprocess.internal(object=object,xstart=xstart,times=times,params=params,offset=offset,...)
    }
)
