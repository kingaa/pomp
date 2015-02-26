## initialize the state variables of the process model

init.state.internal <- function (object, params, t0, ...) {
  if (missing(t0)) t0 <- object@t0
  if (missing(params)) params <- coef(object)
  pompLoad(object)
  rv <- .Call(do_init_state,object,params,t0)
  pompUnload(object)
  rv
}

setMethod('init.state',
          signature=signature('pomp'),
          definition=function (object, params, t0, ...) {
            init.state.internal(object=object,params=params,t0=t0,...)
          }
          )
