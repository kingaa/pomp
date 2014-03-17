## initialize the state variables of the process model

init.state.internal <- function (object, params, t0, ...) {
  if (missing(t0)) t0 <- object@t0
  if (missing(params)) params <- coef(object)
  .Call(do_init_state,object,params,t0)
}

setMethod(
          'init.state',
          'pomp',
          function (object, params, t0, ...) {
            init.state.internal(object=object,params=params,t0=t0,...)
          }
          )
