## evaluate the prior probability density

dprior.internal <- function (object, params, log = FALSE,
                             .getnativesymbolinfo = TRUE, ...) {
  pompLoad(object)
  rv <- .Call(do_dprior,object,params,log,.getnativesymbolinfo)
  pompUnload(object)
  rv
}

setMethod("dprior",
          signature=signature("pomp"),
          function (object, params, log = FALSE, ...)
          dprior.internal(object=object,params=params,log=log,...)
          )
