## initialize the state variables of the process model

init.state.internal <- function (object, params, t0,
                                 .getnativesymbolinfo = TRUE, ...) {
  if (missing(t0)) t0 <- object@t0
  if (missing(params)) params <- coef(object)
  pompLoad(object)
  x <- .Call(do_init_state,object,params,t0,.getnativesymbolinfo)
  pompUnload(object)
  x
}

setMethod('init.state',
          signature=signature('pomp'),
          definition=function (object, params, t0, ...) {
            init.state.internal(object=object,params=params,t0=t0,...)
          }
          )
