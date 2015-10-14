## simulate the process model

rprocess.internal <- function (object, xstart, times, params, offset = 0, .getnativesymbolinfo = TRUE, ...) {
    pompLoad(object)
    rv <- .Call(do_rprocess,object,xstart,times,params,offset,.getnativesymbolinfo)
    pompUnload(object)
    rv
}

setMethod(
    "rprocess",
    signature=signature(object="pomp"),
    definition=function (object, xstart, times, params, offset = 0, ...) {
        rprocess.internal(object=object,xstart=xstart,times=times,params=params,offset=offset,...)
    }
)
