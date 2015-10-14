## simulate the measurement model

rmeasure.internal <- function (object, x, times, params,
                               .getnativesymbolinfo = TRUE, ...) {
    pompLoad(object)
    rv <- .Call(do_rmeasure,object,x,times,params,.getnativesymbolinfo)
    pompUnload(object)
    rv
}

setMethod(
    "rmeasure",
    signature=signature(object="pomp"),
    definition=function (object, x, times, params, ...)
        rmeasure.internal(object=object,x=x,times=times,params=params,...)
)
