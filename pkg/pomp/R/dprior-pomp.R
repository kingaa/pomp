## evaluate the prior probability density

dprior.internal <- function (object, params, log = FALSE,
                             .getnativesymbolinfo = TRUE, ...) {
  .Call(do_dprior,object,params,log,.getnativesymbolinfo)
}

setMethod("dprior","pomp",
          function (object, params, log = FALSE, ...)
          dprior.internal(object=object,params=params,log=log,...)
          )
